# JFEM - Julia Nastran-compatible FEM Solver (SOL 101 Linear Static)
#
# Usage:  julia main.jl <model.bdf> [output_dir]
#
# Supported elements: CQUAD4, CTRIA3, CBAR (PBARL ROD), RBE3
# Supported loads:    FORCE, MOMENT, PLOAD4, LOAD combos
# Supported BCs:      SPC1, SPCADD
# Output:             VTK (.vtu) + JSON results

using LinearAlgebra
using SparseArrays
using Printf
using WriteVTK
using Statistics
using JSON

include("FEMKernels.jl"); using .FEM
include("NastranMath.jl"); using .NastranMath
include("NastranParser.jl"); using .NastranParser
include("Solver.jl"); using .Solver

function transform_geometry!(model)
    grids = model["GRIDs"]
    cords = model["CORDs"]
    for (sid, g) in grids
        if g["CP"] != 0 && haskey(cords, string(g["CP"]))
            c = cords[string(g["CP"])]
            g["X"] = NastranMath.CORDR_transform(c["Origin"], c["U"], c["V"], c["W"], g["X"])
        end
    end
end

function main(filename::String, output_dir::String="")
    if !isfile(filename)
        println("ERROR: File not found: $filename")
        return
    end

    if isempty(output_dir)
        output_dir = joinpath(dirname(filename), "FEM_output")
    end
    if !isdir(output_dir)
        mkpath(output_dir)
    end

    println(">>> Reading BDF file: $filename")
    lines = readlines(filename)

    println(">>> Parsing Bulk Data...")
    cc, bulk = NastranParser.read_bulk_and_case(lines)
    cards = NastranParser.process_cards(bulk)

    println(">>> Constructing Model Data...")
    model = Dict(
        "CASE_CONTROL" => cc,
        "GRIDs"       => haskey(cards,"GRID")   ? NastranParser.extract_grid(cards["GRID"]) : Dict(),
        "CORDs"       => merge(haskey(cards,"CORD2R") ? NastranParser.extract_coords(cards["CORD2R"]) : Dict(),
                               haskey(cards,"CORD1R") ? NastranParser.extract_coords(cards["CORD1R"]) : Dict()),
        "CSHELLs"     => merge(haskey(cards,"CTRIA3") ? NastranParser.extract_shells(cards["CTRIA3"]) : Dict(),
                               haskey(cards,"CQUAD4") ? NastranParser.extract_shells(cards["CQUAD4"]) : Dict()),
        "CBARs"       => haskey(cards,"CBAR")   ? NastranParser.extract_cbar(cards["CBAR"]) : Dict(),
        "RBE3s"       => haskey(cards,"RBE3")   ? NastranParser.extract_rbe3(cards["RBE3"]) : Dict(),
        "PSHELLs"     => haskey(cards,"PSHELL") ? NastranParser.extract_props_shell(cards["PSHELL"]) : Dict(),
        "PBARLs"      => haskey(cards,"PBARL")  ? NastranParser.extract_pbarl(cards["PBARL"]) : Dict(),
        "MATs"        => haskey(cards,"MAT1")   ? NastranParser.extract_mats(cards["MAT1"]) : Dict(),
        "FORCEs"      => haskey(cards,"FORCE")  ? NastranParser.extract_loads(cards["FORCE"]) : [],
        "MOMENTs"     => haskey(cards,"MOMENT") ? NastranParser.extract_moments(cards["MOMENT"]) : [],
        "PLOAD4s"     => haskey(cards,"PLOAD4") ? NastranParser.extract_pload4(cards["PLOAD4"]) : [],
        "LOAD_COMBOS" => haskey(cards,"LOAD")   ? NastranParser.extract_load_combos(cards["LOAD"]) : [],
        "SPC1s"       => haskey(cards,"SPC1")   ? NastranParser.extract_spc1(cards["SPC1"]) : [],
        "SPCADDs"     => haskey(cards,"SPCADD") ? NastranParser.extract_spcadd(cards["SPCADD"]) : Dict()
    )

    if haskey(cards, "PARAM")
        for c in cards["PARAM"]
            pname = uppercase(strip(string(NastranParser.safe_get(c, 3))))
            pval = NastranParser.parse_nastran_number(NastranParser.safe_get(c, 4), 0.0)
            model["PARAM_$pname"] = pval
        end
    end

    transform_geometry!(model)
    K, id_map, X, ndof, node_R = Solver.assemble_stiffness(model)

    global_results = Dict(
        "displacements" => [],
        "spc_forces" => [],
        "forces" => Dict("cbar" => [], "quad4" => [], "tria3" => []),
        "stresses" => Dict("cbar" => [], "quad4" => [], "tria3" => []),
        "strains" => Dict("cbar" => [], "quad4" => [], "tria3" => [])
    )

    sorted_sids = sort(collect(keys(cc["SUBCASES"])))

    for sid in sorted_sids
        sub = cc["SUBCASES"][sid]
        println("\n>>> Solving Subcase $sid...")
        load_id = get(sub, "LOAD", nothing)
        spc_id = get(sub, "SPC", nothing)
        u, stresses, sub_res = Solver.solve_case(K, ndof, model, id_map, X, load_id, spc_id, node_R)

        append!(global_results["displacements"], sub_res["displacements"])
        append!(global_results["spc_forces"], sub_res["spc_forces"])
        for k in keys(sub_res["forces"]); append!(global_results["forces"][k], sub_res["forces"][k]); end
        for k in keys(sub_res["stresses"]); append!(global_results["stresses"][k], sub_res["stresses"][k]); end
        for k in keys(sub_res["strains"]); append!(global_results["strains"][k], sub_res["strains"][k]); end

        # VTK export
        base_name = replace(basename(filename), ".bdf" => "")
        vtk_base = base_name * "_Subcase_$sid"
        vtk_path = joinpath(output_dir, vtk_base)
        points = zeros(3, length(id_map))
        disp = zeros(3, length(id_map))
        for (nid, idx) in id_map
            points[:, idx] = X[idx, :]
            disp[:, idx] = u[(idx-1)*6+1:(idx-1)*6+3]
        end
        cells = MeshCell[]
        data_vonmises = Float64[]
        for (id, el) in model["CSHELLs"]
            if !haskey(el, "NODES"); continue; end
            eid = parse(Int, id)
            nids = [get(id_map, n, 0) for n in el["NODES"]]
            if 0 in nids; continue; end
            if length(nids) == 3
                push!(cells, MeshCell(VTKCellTypes.VTK_TRIANGLE, nids))
            elseif length(nids) == 4
                push!(cells, MeshCell(VTKCellTypes.VTK_QUAD, nids))
            end
            push!(data_vonmises, get(stresses, eid, 0.0))
        end
        for (id, bar) in model["CBARs"]
            if !haskey(bar, "GA"); continue; end
            eid = parse(Int, id)
            nids = [get(id_map, bar["GA"], 0), get(id_map, bar["GB"], 0)]
            if 0 in nids; continue; end
            push!(cells, MeshCell(VTKCellTypes.VTK_LINE, nids))
            push!(data_vonmises, get(stresses, eid, 0.0))
        end
        if !isempty(cells)
            vtk = vtk_grid(vtk_path, points, cells)
            vtk["Displacement", VTKPointData()] = disp
            vtk["VonMises_Stress", VTKCellData()] = data_vonmises
            vtk_save(vtk)
            println("  VTK saved: $vtk_path.vtu")
        end
    end

    # JSON export
    base_name = replace(basename(filename), ".bdf" => "")
    json_name = base_name * ".JU.JSON"
    json_path = joinpath(output_dir, json_name)
    println("\n>>> Exporting JSON: $json_path")
    open(json_path, "w") do f; JSON.print(f, global_results, 4); end
    println(">>> Done.")
end

# Entry point
if length(ARGS) >= 1
    bdf_file = ARGS[1]
    out_dir = length(ARGS) >= 2 ? ARGS[2] : ""
    main(bdf_file, out_dir)
else
    println("Usage: julia main.jl <model.bdf> [output_dir]")
    println("  output_dir defaults to FEM_output/ next to the BDF file")
end
