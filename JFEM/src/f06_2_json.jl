# F06 to JSON Converter for Nastran results
# Usage: julia f06_2_json.jl <model.f06> [output_dir]
# Parses Nastran F06 output and writes a JSON file with structured results

using JSON

function parse_nastran_float(s::AbstractString)
    s = strip(s)
    if isempty(s); return 0.0; end
    val = tryparse(Float64, s)
    if val !== nothing; return val; end
    try
        for i in length(s):-1:2
            c = s[i]
            prev = s[i-1]
            if (c == '+' || c == '-') && (prev != 'E' && prev != 'e')
                new_s = s[1:i-1] * "E" * s[i:end]
                return parse(Float64, new_s)
            end
        end
    catch; end
    return 0.0
end

function is_header(clean_line::AbstractString, keyword::String)
    normalized = replace(clean_line, " " => "")
    return occursin(keyword, normalized)
end

function parse_f06_to_json(file_content::String)
    results = Dict(
        "displacements" => [],
        "spc_forces" => [],
        "forces" => Dict("cbar" => [], "quad4" => [], "tria3" => []),
        "stresses" => Dict("cbar" => [], "quad4" => [], "tria3" => []),
        "strains" => Dict("cbar" => [], "quad4" => [], "tria3" => [])
    )

    current_section = :none
    cbar_buffer = Dict()

    lines = split(file_content, "\n")
    println("... Scanning $(length(lines)) lines...")

    for line in lines
        if startswith(line, "1"); continue; end
        raw_line = string(line)
        line_stripped = strip(raw_line)
        if isempty(line_stripped); continue; end

        if is_header(line_stripped, "DISPLACEMENTVECTOR")
            current_section = :displacement; continue
        elseif is_header(line_stripped, "FORCESOFSINGLE-POINTCONSTRAINT")
            current_section = :spc_force; continue
        elseif is_header(line_stripped, "FORCESINBARELEMENTS")
            current_section = :force_cbar; continue
        elseif is_header(line_stripped, "FORCESINQUADRILATERAL")
            current_section = :force_quad4; continue
        elseif is_header(line_stripped, "FORCESINTRIANGULAR")
            current_section = :force_tria3; continue
        elseif is_header(line_stripped, "STRESSESINBARELEMENTS")
            current_section = :stress_cbar; continue
        elseif is_header(line_stripped, "STRAINSINBARELEMENTS")
            current_section = :strain_cbar; continue
        elseif is_header(line_stripped, "STRESSESINQUADRILATERAL")
            current_section = :stress_quad4; continue
        elseif is_header(line_stripped, "STRAINSINQUADRILATERAL")
            current_section = :strain_quad4; continue
        elseif is_header(line_stripped, "STRESSESINTRIANGULAR")
            current_section = :stress_tria3; continue
        elseif is_header(line_stripped, "STRAINSINTRIANGULAR")
            current_section = :strain_tria3; continue
        end

        if is_header(line_stripped, "POINTID") || is_header(line_stripped, "ELEMENTID"); continue; end
        if is_header(line_stripped, "SUBCASE") || is_header(line_stripped, "DAREA"); continue; end
        if is_header(line_stripped, "BEND-MOMENT") || occursin("MSC Nastran", line_stripped); continue; end

        try
            parts = split(line_stripped)
            if isempty(parts); continue; end

            if parts[1] == "0"
                popfirst!(parts)
            elseif startswith(parts[1], "0") && length(parts[1]) > 1 && isdigit(parts[1][2])
                parts[1] = parts[1][2:end]
            end
            if isempty(parts); continue; end

            if current_section == :displacement
                id = tryparse(Int, parts[1])
                if id !== nothing
                    idx = (length(parts) == 8) ? 3 : 2
                    if length(parts) >= idx+5
                        push!(results["displacements"], Dict(
                            "grid_id" => id,
                            "t1" => parse_nastran_float(parts[idx]),
                            "t2" => parse_nastran_float(parts[idx+1]),
                            "t3" => parse_nastran_float(parts[idx+2]),
                            "r1" => parse_nastran_float(parts[idx+3]),
                            "r2" => parse_nastran_float(parts[idx+4]),
                            "r3" => parse_nastran_float(parts[idx+5])
                        ))
                    end
                end

            elseif current_section == :spc_force
                id = tryparse(Int, parts[1])
                if id !== nothing
                    idx = (length(parts) == 8) ? 3 : 2
                    if length(parts) >= idx+5
                        push!(results["spc_forces"], Dict(
                            "grid_id" => id,
                            "t1" => parse_nastran_float(parts[idx]),
                            "t2" => parse_nastran_float(parts[idx+1]),
                            "t3" => parse_nastran_float(parts[idx+2]),
                            "r1" => parse_nastran_float(parts[idx+3]),
                            "r2" => parse_nastran_float(parts[idx+4]),
                            "r3" => parse_nastran_float(parts[idx+5])
                        ))
                    end
                end

            elseif current_section == :force_cbar
                id = tryparse(Int, parts[1])
                if id !== nothing && length(parts) >= 9
                    push!(results["forces"]["cbar"], Dict(
                        "eid" => id,
                        "moment_a1" => parse_nastran_float(parts[2]),
                        "moment_a2" => parse_nastran_float(parts[3]),
                        "moment_b1" => parse_nastran_float(parts[4]),
                        "moment_b2" => parse_nastran_float(parts[5]),
                        "shear_1"   => parse_nastran_float(parts[6]),
                        "shear_2"   => parse_nastran_float(parts[7]),
                        "axial"     => parse_nastran_float(parts[8]),
                        "torque"    => parse_nastran_float(parts[9])
                    ))
                end

            elseif (current_section == :force_quad4 || current_section == :force_tria3)
                key = (current_section == :force_quad4) ? "quad4" : "tria3"
                id = tryparse(Int, parts[1])
                if id !== nothing && length(parts) >= 9
                    push!(results["forces"][key], Dict(
                        "eid" => id,
                        "fx" => parse_nastran_float(parts[2]), "fy" => parse_nastran_float(parts[3]),
                        "fxy" => parse_nastran_float(parts[4]), "mx" => parse_nastran_float(parts[5]),
                        "my" => parse_nastran_float(parts[6]), "mxy" => parse_nastran_float(parts[7]),
                        "qx" => parse_nastran_float(parts[8]), "qy" => parse_nastran_float(parts[9])
                    ))
                end

            elseif current_section in [:stress_cbar, :strain_cbar]
                cat_key = (current_section == :stress_cbar) ? "stresses" : "strains"
                id = tryparse(Int, parts[1])
                if id !== nothing && length(parts) >= 6
                    cbar_buffer = Dict(
                        "eid" => id,
                        "end_a" => Dict(
                            "p1" => parse_nastran_float(parts[2]), "p2" => parse_nastran_float(parts[3]),
                            "p3" => parse_nastran_float(parts[4]), "p4" => parse_nastran_float(parts[5])
                        ),
                        "axial" => parse_nastran_float(parts[6])
                    )
                elseif !isempty(cbar_buffer) && length(parts) >= 4
                    cbar_buffer["end_b"] = Dict(
                        "p1" => parse_nastran_float(parts[1]), "p2" => parse_nastran_float(parts[2]),
                        "p3" => parse_nastran_float(parts[3]), "p4" => parse_nastran_float(parts[4])
                    )
                    push!(results[cat_key]["cbar"], deepcopy(cbar_buffer))
                    empty!(cbar_buffer)
                end

            elseif current_section in [:stress_quad4, :stress_tria3, :strain_quad4, :strain_tria3]
                is_stress = occursin("stress", string(current_section))
                elem_type = occursin("quad", string(current_section)) ? "quad4" : "tria3"
                root_cat = is_stress ? "stresses" : "strains"

                id_val = tryparse(Int, parts[1])
                is_id = (id_val !== nothing) && !occursin(".", parts[1])

                if is_id && length(parts) >= 9
                    entry = Dict(
                        "fiber_dist" => parse_nastran_float(parts[2]),
                        "normal_x"   => parse_nastran_float(parts[3]),
                        "normal_y"   => parse_nastran_float(parts[4]),
                        "shear_xy"   => parse_nastran_float(parts[5]),
                        "angle"      => parse_nastran_float(parts[6]),
                        "major"      => parse_nastran_float(parts[7]),
                        "minor"      => parse_nastran_float(parts[8]),
                        "von_mises"  => parse_nastran_float(parts[9])
                    )
                    push!(results[root_cat][elem_type], Dict("eid" => id_val, "z1" => entry))

                elseif !is_id && !isempty(results[root_cat][elem_type]) && length(parts) >= 8
                    entry = Dict(
                        "fiber_dist" => parse_nastran_float(parts[1]),
                        "normal_x"   => parse_nastran_float(parts[2]),
                        "normal_y"   => parse_nastran_float(parts[3]),
                        "shear_xy"   => parse_nastran_float(parts[4]),
                        "angle"      => parse_nastran_float(parts[5]),
                        "major"      => parse_nastran_float(parts[6]),
                        "minor"      => parse_nastran_float(parts[7]),
                        "von_mises"  => parse_nastran_float(parts[8])
                    )
                    results[root_cat][elem_type][end]["z2"] = entry
                end
            end
        catch e
        end
    end

    return results
end

function convert_f06(f06_path::String, output_dir::String="")
    if !isfile(f06_path)
        println("ERROR: File not found: $f06_path")
        return
    end

    if isempty(output_dir)
        output_dir = joinpath(dirname(f06_path), "FEM_output")
    end
    if !isdir(output_dir)
        mkpath(output_dir)
    end

    model_name = replace(basename(f06_path), ".f06" => "", ".F06" => "")
    json_filename = "$(model_name).NAST.json"
    output_path = joinpath(output_dir, json_filename)

    println(">>> Reading F06: $f06_path")
    file_content = read(f06_path, String)

    println("... Parsing content")
    json_data = parse_f06_to_json(file_content)

    println("... Extracted $(length(json_data["displacements"])) displacement entries")
    println("... Extracted $(length(json_data["spc_forces"])) spc force entries")
    println("... Extracted $(length(json_data["forces"]["cbar"])) cbar force entries")
    println("... Extracted $(length(json_data["stresses"]["quad4"])) quad4 stress entries")

    println("... Writing JSON to: $output_path")
    open(output_path, "w") do f
        JSON.print(f, json_data, 4)
    end
    println(">>> Success!")
end

# Entry point
if length(ARGS) >= 1
    f06_file = ARGS[1]
    out_dir = length(ARGS) >= 2 ? ARGS[2] : ""
    convert_f06(f06_file, out_dir)
else
    println("Usage: julia f06_2_json.jl <model.f06> [output_dir]")
    println("  Converts Nastran F06 output to JSON format")
    println("  output_dir defaults to FEM_output/ next to the F06 file")
end
