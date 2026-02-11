module Solver

using LinearAlgebra
using SparseArrays
using Statistics
using ..FEM
using ..NastranMath

function shell_element_frame(Xc, n)
    if n == 4
        # Nastran CQUAD4 element coordinate system:
        # x_e bisects the angle between diagonals G1-G3 and G2-G4
        # z_e = normal from right-hand rule (G1->G3) x (G2->G4)
        d13 = Xc[3,:] - Xc[1,:]
        d24 = Xc[4,:] - Xc[2,:]
        v3 = normalize(cross(d13, d24))
        # Bisector of the two diagonals: difference of unit vectors
        # (d13n - d24n gives the xi-like bisector direction)
        d13n = normalize(d13)
        d24n = normalize(d24)
        x_raw = d13n - d24n
        if norm(x_raw) < 1e-10
            x_raw = d13n + d24n
        end
        x_raw = normalize(x_raw)
        v1 = normalize(x_raw - dot(x_raw, v3)*v3)
        v2 = cross(v3, v1)
    else
        # CTRIA3: x_e = edge 1->2
        v1 = normalize(Xc[2,:] - Xc[1,:])
        v3 = normalize(cross(v1, Xc[3,:] - Xc[1,:]))
        v2 = cross(v3, v1)
    end
    return v1, v2, v3
end

function get_nodal_transform(model, cid)
    if cid == 0; return Matrix(1.0I, 3, 3); end
    if !haskey(model["CORDs"], string(cid)); return Matrix(1.0I, 3, 3); end
    cord = model["CORDs"][string(cid)]
    return hcat(cord["U"], cord["V"], cord["W"])
end

function get_coord_transform(model, cid, vec)
    if cid == 0; return vec; end
    if !haskey(model["CORDs"], string(cid)); return vec; end
    cord = model["CORDs"][string(cid)]
    return NastranMath.CORDR_transform_vector(cord["U"], cord["V"], cord["W"], vec)
end

function rotation_from_normal(n_avg::Vector{Float64})
    z = normalize(n_avg)
    ref = abs(z[3]) < 0.9 ? [0.0, 0.0, 1.0] : [1.0, 0.0, 0.0]
    x = normalize(cross(ref, z))
    y = cross(z, x)
    return hcat(x, y, z)
end

function compute_snorm(model, id_map, node_coords, node_R)
    snorm_angle = get(model, "PARAM_SNORM", 20.0)
    if snorm_angle <= 0.0
        return 0
    end

    node_elems = Dict{Int, Vector{Tuple{String, Vector{Float64}}}}()
    for (eid, el) in model["CSHELLs"]
        if !haskey(el, "NODES"); continue; end
        nids_mapped = [get(id_map, n, 0) for n in el["NODES"]]
        if any(x->x==0, nids_mapped); continue; end
        Xc = node_coords[nids_mapped, :]
        n = length(nids_mapped)
        if n == 4
            d13 = Xc[3,:] - Xc[1,:]; d24 = Xc[4,:] - Xc[2,:]
            nrm = cross(d13, d24)
        else
            v1e = Xc[2,:] - Xc[1,:]; v2e = Xc[3,:] - Xc[1,:]
            nrm = cross(v1e, v2e)
        end
        nn = norm(nrm)
        if nn < 1e-12; continue; end
        nrm = nrm / nn
        for nidx in nids_mapped
            if !haskey(node_elems, nidx); node_elems[nidx] = []; end
            push!(node_elems[nidx], (eid, nrm))
        end
    end

    modified_avg = 0
    for (nidx, elem_list) in node_elems
        if length(elem_list) == 1
            eid, nrm = elem_list[1]
            node_R[nidx] = rotation_from_normal(nrm)
            modified_avg += 1
            continue
        end
        avg_n = copy(elem_list[1][2])
        for i in 2:length(elem_list)
            nrm = elem_list[i][2]
            if dot(nrm, avg_n) >= 0
                avg_n .+= nrm
            else
                avg_n .-= nrm
            end
        end
        nn = norm(avg_n)
        if nn < 1e-12; continue; end
        avg_n ./= nn
        all_within = true
        for (_, nrm) in elem_list
            angle = acosd(clamp(abs(dot(nrm, avg_n)), 0, 1))
            if angle > snorm_angle
                all_within = false
                break
            end
        end
        if all_within
            node_R[nidx] = rotation_from_normal(avg_n)
            modified_avg += 1
        end
    end

    println(">>> [SOLVER] SNORM: $modified_avg averaged nodes")
    return modified_avg
end

function resolve_loads(model, sid, scale, id_map, elem_map, node_coords, F_acc)
    raw_forces = Dict{Int, Vector{Float64}}()
    add_force = (gid, vec) -> begin
        if !haskey(raw_forces, gid); raw_forces[gid] = zeros(6); end
        raw_forces[gid] .+= vec
    end

    for frc in model["FORCEs"]
        if Int(frc["SID"]) == sid
            global_dir = get_coord_transform(model, Int(frc["CID"]), frc["Dir"])
            v = normalize(global_dir) * frc["Mag"] * scale
            f_vec = zeros(6); f_vec[1:3] = v
            add_force(frc["GID"], f_vec)
        end
    end

    for mom in model["MOMENTs"]
        if Int(mom["SID"]) == sid
            global_dir = get_coord_transform(model, Int(mom["CID"]), mom["Dir"])
            v = normalize(global_dir) * mom["Mag"] * scale
            f_vec = zeros(6); f_vec[4:6] = v
            add_force(mom["GID"], f_vec)
        end
    end

    for pload in model["PLOAD4s"]
        if Int(pload["SID"]) == sid
            eid = pload["EID"]
            el_def = nothing
            if haskey(model["CSHELLs"], string(eid)); el_def = model["CSHELLs"][string(eid)]; end
            if !isnothing(el_def)
                nids = [get(id_map, n, 0) for n in el_def["NODES"]]
                if !any(x->x==0, nids)
                    Xc = node_coords[nids, :]
                    v1 = Xc[2,:] - Xc[1,:]
                    v2 = Xc[3,:] - Xc[1,:]
                    if length(nids) == 4
                        d1 = Xc[3,:] - Xc[1,:]
                        d2 = Xc[4,:] - Xc[2,:]
                        normal_vec = cross(d1, d2); area = 0.5 * norm(normal_vec)
                    else
                        normal_vec = cross(v1, v2); area = 0.5 * norm(normal_vec)
                    end
                    unit_normal = normalize(normal_vec)
                    total_force = area * pload["P"] * scale
                    f_node = unit_normal * (total_force / length(nids))
                    for idx in nids
                        dof = (idx-1)*6; F_acc[dof+1 : dof+3] .+= f_node
                    end
                end
            end
        end
    end

    # RBE3 load distribution
    rbe_map = Dict{Int, Any}()
    for (id, rbe) in model["RBE3s"]
        rbe_map[rbe["REFGRID"]] = rbe
    end

    for (gid, f_vec) in raw_forces
        if haskey(rbe_map, gid)
            rbe = rbe_map[gid]
            deps = rbe["DEP_GRIDS"]
            N = length(deps)
            if N > 0
                ref_idx = get(id_map, gid, 0)
                X_ref = (ref_idx > 0) ? node_coords[ref_idx, :] : zeros(3)
                dep_indices = Int[]; dep_coords = Vector{Float64}[]
                for dep_id in deps
                    idx = get(id_map, dep_id, 0)
                    if idx > 0
                        push!(dep_indices, idx)
                        push!(dep_coords, node_coords[idx, :])
                    end
                end
                Nd = length(dep_indices)
                if Nd == 0; continue; end
                centroid = sum(dep_coords) ./ Nd
                F_ref = f_vec[1:3]
                M_ref = f_vec[4:6]
                M_total = M_ref .+ cross(X_ref .- centroid, F_ref)
                F_per = F_ref ./ Nd
                A = zeros(3, 3*Nd)
                for i in 1:Nd
                    r = dep_coords[i] .- centroid
                    j = (i-1)*3
                    A[1, j+2] = -r[3]; A[1, j+3] = r[2]
                    A[2, j+1] = r[3];  A[2, j+3] = -r[1]
                    A[3, j+1] = -r[2]; A[3, j+2] = r[1]
                end
                AAT = A * A'
                if abs(det(AAT)) > 1e-20
                    f_all = A' * (AAT \ M_total)
                else
                    f_all = zeros(3*Nd)
                end
                for i in 1:Nd
                    idx = dep_indices[i]
                    dof = (idx-1)*6
                    F_acc[dof+1:dof+3] .+= F_per .+ f_all[(i-1)*3+1:(i-1)*3+3]
                end
            end
        else
            if haskey(id_map, gid)
                idx = id_map[gid]
                dof = (idx-1)*6
                F_acc[dof+1:dof+6] .+= f_vec
            end
        end
    end

    for c in model["LOAD_COMBOS"]
        if Int(c["SID"]) == sid
            for sub in c["COMPS"]
                resolve_loads(model, Int(sub["LID"]), scale * c["S"] * sub["S"], id_map, elem_map, node_coords, F_acc)
            end
        end
    end
end

function assemble_stiffness(model)
    println(">>> [SOLVER] Indexing...")
    ids = sort(collect(keys(model["GRIDs"])), by=x->parse(Int,x))
    n_nodes = length(ids)
    id_map = Dict(parse(Int, k)=>i for (i,k) in enumerate(ids))
    ndof = n_nodes * 6

    node_R = Vector{Matrix{Float64}}(undef, n_nodes)
    node_coords = zeros(n_nodes, 3)
    for (sid, g) in model["GRIDs"]
        idx = id_map[g["ID"]]
        node_coords[idx, :] = g["X"]
        node_R[idx] = get_nodal_transform(model, g["CD"])
    end

    compute_snorm(model, id_map, node_coords, node_R)

    println(">>> [SOLVER] Assembling K (Size: $ndof)...")
    I_idx = Int[]; J_idx = Int[]; V_val = Float64[]

    for (id, el) in model["CSHELLs"]
        if !haskey(model["PSHELLs"], string(el["PID"])); continue; end
        prop = model["PSHELLs"][string(el["PID"])]
        if !haskey(model["MATs"], string(prop["MID"])); continue; end
        mat = model["MATs"][string(prop["MID"])]
        nids = [get(id_map, n, 0) for n in el["NODES"]]
        if any(x->x==0, nids); continue; end
        Xc = node_coords[nids, :]
        n = length(nids)
        v1, v2, v3 = shell_element_frame(Xc, n)
        R_el = hcat(v1, v2, v3)
        T = zeros(6*n, 6*n)
        for k in 1:n
            R_node = node_R[nids[k]]
            TR = R_el' * R_node
            base = (k-1)*6; T[base+1:base+3, base+1:base+3] = TR; T[base+4:base+6, base+4:base+6] = TR
        end
        lc = zeros(n, 2); c = mean(Xc, dims=1)[:]
        for k=1:n; d = Xc[k,:] - c; lc[k,1] = dot(d, v1); lc[k,2] = dot(d, v2); end
        br = get(prop, "BEND_RATIO", 1.0); tst = get(prop, "TS_T", 5.0/6.0)
        Ke_loc = (n == 4) ? FEM.stiffness_quad4(lc, mat["E"], mat["NU"], prop["T"]; bend_ratio=br, ts_t=tst) : FEM.stiffness_tria3(lc, mat["E"], mat["NU"], prop["T"])
        Ke = T' * Ke_loc * T
        dofs = reduce(vcat, [(idx-1)*6 .+ (1:6) for idx in nids])
        append!(I_idx, repeat(dofs, inner=6*n)); append!(J_idx, repeat(dofs, outer=6*n)); append!(V_val, vec(Ke))
    end

    for (id, bar) in model["CBARs"]
        if !haskey(model["PBARLs"], string(bar["PID"])); continue; end
        prop = model["PBARLs"][string(bar["PID"])]
        if !haskey(model["MATs"], string(prop["MID"])); continue; end
        mat = model["MATs"][string(prop["MID"])]
        ia, ib = get(id_map, bar["GA"], 0), get(id_map, bar["GB"], 0)
        if ia==0 || ib==0; continue; end
        Pa, Pb = node_coords[ia,:], node_coords[ib,:]
        L = norm(Pb - Pa)
        vx = normalize(Pb - Pa); v_ref = bar["V"]
        if norm(v_ref) < 1e-6; v_ref=[0.0,0.0,1.0]; if abs(dot(vx, v_ref)) > 0.9; v_ref=[0.0,1.0,0.0]; end; end
        vz = normalize(cross(vx, v_ref)); vy = cross(vz, vx); R_el = hcat(vx, vy, vz)
        T = zeros(12,12)
        Ra = node_R[ia]; Tra = R_el' * Ra; T[1:3, 1:3] = Tra; T[4:6, 4:6] = Tra
        Rb = node_R[ib]; Trb = R_el' * Rb; T[7:9, 7:9] = Trb; T[10:12, 10:12] = Trb
        As_bar = Inf
        if get(prop, "K_SHEAR", 0.0) > 0.0
            nu_bar = mat["NU"]
            kappa_s = 6*(1+nu_bar)/(7+6*nu_bar)
            As_bar = kappa_s * prop["A"]
        end
        Ke = T' * FEM.stiffness_frame3d(L, prop["A"], prop["I"], prop["I"], prop["J"], mat["E"], mat["G"]; As_y=As_bar, As_z=As_bar) * T
        dofs = vcat((ia-1)*6 .+ (1:6), (ib-1)*6 .+ (1:6))
        append!(I_idx, repeat(dofs, inner=12)); append!(J_idx, repeat(dofs, outer=12)); append!(V_val, vec(Ke))
    end

    K = sparse(I_idx, J_idx, V_val, ndof, ndof)
    return K, id_map, node_coords, ndof, node_R
end

function solve_case(K, ndof, model, id_map, X, load_id, spc_id, node_R)
    F_applied = zeros(ndof)
    if !isnothing(load_id)
        elem_map = Dict{Int, Any}()
        F_global_accum = zeros(ndof)
        resolve_loads(model, Int(load_id), 1.0, id_map, elem_map, X, F_global_accum)
        for i in 1:length(id_map)
            idx = (i-1)*6
            R = node_R[i]
            F_applied[idx+1:idx+3] = R' * F_global_accum[idx+1:idx+3]
            F_applied[idx+4:idx+6] = R' * F_global_accum[idx+4:idx+6]
        end
    end
    println("      Applied Force Norm: $(norm(F_applied))")

    fixed = Set{Int}()
    sets = Set{Int}()
    if !isnothing(spc_id)
        sid = Int(spc_id)
        if haskey(model["SPCADDs"], sid); union!(sets, model["SPCADDs"][sid]); else; push!(sets, sid); end
    end
    for spc in model["SPC1s"]
        if Int(spc["SID"]) in sets
            for n in spc["NODES"]
                idx = get(id_map, n, 0)
                if idx > 0
                    for c in spc["C"]; push!(fixed, (idx-1)*6 + parse(Int, c)); end
                end
            end
        end
    end

    free = setdiff(1:ndof, fixed)
    u = zeros(ndof)

    K_diag = diag(K)
    non_zeros = filter(x -> x > 1e-6, K_diag)
    avg_k = isempty(non_zeros) ? 1.0 : mean(non_zeros)
    for i in free
        if abs(K[i, i]) < 1e-9 * avg_k; K[i, i] += avg_k * 1e-6; end
        K[i, i] += avg_k * 1e-12
    end

    u[free] = K[free, free] \ F_applied[free]

    # RBE3 displacement interpolation
    for (id, rbe) in model["RBE3s"]
        ref_idx = get(id_map, rbe["REFGRID"], 0)
        if ref_idx == 0; continue; end
        dep_grids = rbe["DEP_GRIDS"]
        X_ref = X[ref_idx, :]
        n_deps = 0; u_avg = zeros(6); centroid = zeros(3)
        dep_list = Int[]
        for gid in dep_grids
            dep_idx = get(id_map, gid, 0)
            if dep_idx == 0; continue; end
            base = (dep_idx - 1) * 6
            u_avg .+= u[base+1:base+6]
            centroid .+= X[dep_idx, :]
            push!(dep_list, dep_idx)
            n_deps += 1
        end
        if n_deps > 0
            u_avg ./= n_deps
            centroid ./= n_deps
            offset = X_ref .- centroid
            omega_avg = u_avg[4:6]
            u_ref_trans = u_avg[1:3] .+ cross(omega_avg, offset)
            ref_base = (ref_idx - 1) * 6
            u[ref_base+1:ref_base+3] = u_ref_trans
            u[ref_base+4:ref_base+6] = u_avg[4:6]
        end
    end

    R = K * u - F_applied

    u_global = zeros(ndof)
    for i in 1:length(id_map)
        idx = (i-1)*6
        R_mat = node_R[i]
        u_local_node = u[idx+1:idx+6]
        u_global[idx+1:idx+3] = R_mat * u_local_node[1:3]
        u_global[idx+4:idx+6] = R_mat * u_local_node[4:6]
    end

    results_json = Dict(
        "displacements" => [],
        "spc_forces" => [],
        "forces" => Dict("cbar" => [], "quad4" => [], "tria3" => []),
        "stresses" => Dict("cbar" => [], "quad4" => [], "tria3" => []),
        "strains" => Dict("cbar" => [], "quad4" => [], "tria3" => [])
    )

    sorted_nodes = sort(collect(keys(id_map)))
    for nid in sorted_nodes
        idx = id_map[nid]; base = (idx-1)*6
        u_json = copy(u_global[base+1:base+6])
        push!(results_json["displacements"], Dict("grid_id" => nid, "t1" => u_json[1], "t2" => u_json[2], "t3" => u_json[3], "r1" => u_json[4], "r2" => u_json[5], "r3" => u_json[6]))
        R_spc_local = R[base+1:base+6]
        if norm(R_spc_local) > 1e-6
            R_mat = node_R[idx]
            spc_glob_t = R_mat * R_spc_local[1:3]
            spc_glob_r = R_mat * R_spc_local[4:6]
            push!(results_json["spc_forces"], Dict("grid_id" => nid, "t1" => spc_glob_t[1], "t2" => spc_glob_t[2], "t3" => spc_glob_t[3], "r1" => spc_glob_r[1], "r2" => spc_glob_r[2], "r3" => spc_glob_r[3]))
        end
    end

    stresses = Dict{Int, Float64}()
    for (id, el) in model["CSHELLs"]
        eid = parse(Int, id)
        if !haskey(model["PSHELLs"], string(el["PID"])); continue; end
        prop = model["PSHELLs"][string(el["PID"])]
        if !haskey(model["MATs"], string(prop["MID"])); continue; end
        mat = model["MATs"][string(prop["MID"])]
        nids = [get(id_map, n, 0) for n in el["NODES"]]; if any(x->x==0, nids); continue; end
        Xc = X[nids, :]; n = length(nids)
        v1, v2, v3 = shell_element_frame(Xc, n)
        R_el = hcat(v1, v2, v3)
        T = zeros(6*n, 6*n)
        for k=0:(n-1); R_node = node_R[nids[k+1]]; TR = R_el' * R_node; T[k*6+1:k*6+3, k*6+1:k*6+3] = TR; T[k*6+4:k*6+6, k*6+4:k*6+6] = TR; end
        lc = zeros(n, 2); c = mean(Xc,dims=1)[:]
        for k=1:n; d=Xc[k,:]-c; lc[k,1]=dot(d,v1); lc[k,2]=dot(d,v2); end
        dofs = reduce(vcat, [(idx-1)*6 .+ (1:6) for idx in nids]); u_loc = T * u[dofs]
        br = get(prop, "BEND_RATIO", 1.0)
        if n == 4
            N, M, Q, s_z1, s_z2, e_z1, e_z2 = FEM.stress_strain_quad4(lc, u_loc, mat["E"], mat["NU"], prop["T"], prop["T"]; bend_ratio=br)
        else
            N, M, Q, s_z1, s_z2, e_z1, e_z2 = FEM.stress_strain_tria3(lc, u_loc, mat["E"], mat["NU"], prop["T"])
        end

        vm = FEM.compute_principal_2d(s_z1[1], s_z1[2], s_z1[3])[1]; stresses[eid] = vm
        elem_key = (n==4) ? "quad4" : "tria3"
        push!(results_json["forces"][elem_key], Dict("eid" => eid, "fx" => N[1], "fy" => N[2], "fxy" => N[3], "mx" => M[1], "my" => M[2], "mxy" => M[3], "qx" => Q[1], "qy" => Q[2]))
        make_stress_entry(s, t) = Dict("fiber_dist" => t, "normal_x" => s[1], "normal_y" => s[2], "shear_xy" => s[3], "von_mises" => sqrt(s[1]^2-s[1]*s[2]+s[2]^2+3*s[3]^2), "major" => 0.0, "minor" => 0.0)
        make_strain_entry(e, t) = Dict("fiber_dist" => t, "normal_x" => e[1], "normal_y" => e[2], "shear_xy" => e[3], "major" => 0.0, "minor" => 0.0)
        push!(results_json["stresses"][elem_key], Dict("eid" => eid, "z1" => make_stress_entry(s_z1, -prop["T"]/2), "z2" => make_stress_entry(s_z2, prop["T"]/2)))
        push!(results_json["strains"][elem_key], Dict("eid" => eid, "z1" => make_strain_entry(e_z1, -prop["T"]/2), "z2" => make_strain_entry(e_z2, prop["T"]/2)))
    end

    for (id, bar) in model["CBARs"]
        eid = parse(Int, id)
        if !haskey(model["PBARLs"], string(bar["PID"])); continue; end
        prop = model["PBARLs"][string(bar["PID"])]; mat = model["MATs"][string(prop["MID"])]
        ia, ib = get(id_map, bar["GA"], 0), get(id_map, bar["GB"], 0); if ia==0 || ib==0; continue; end
        Pa, Pb = X[ia,:], X[ib,:]; L = norm(Pb - Pa); vx = normalize(Pb - Pa); v_ref = bar["V"]
        if norm(v_ref) < 1e-6; v_ref=[0.0,0.0,1.0]; if abs(dot(vx, v_ref)) > 0.9; v_ref=[0.0,1.0,0.0]; end; end
        vz = normalize(cross(vx, v_ref)); vy = cross(vz, vx); R_el = hcat(vx, vy, vz)
        T = zeros(12,12); Ra = node_R[ia]; Tra = R_el' * Ra; T[1:3, 1:3] = Tra; T[4:6, 4:6] = Tra
        Rb = node_R[ib]; Trb = R_el' * Rb; T[7:9, 7:9] = Trb; T[10:12, 10:12] = Trb
        dofs = vcat((ia-1)*6 .+ (1:6), (ib-1)*6 .+ (1:6)); u_loc = T * u[dofs]
        As_bar_f = Inf
        if get(prop, "K_SHEAR", 0.0) > 0.0
            nu_bar_f = mat["NU"]; kappa_sf = 6*(1+nu_bar_f)/(7+6*nu_bar_f); As_bar_f = kappa_sf * prop["A"]
        end
        forces = FEM.forces_frame3d(u_loc, L, prop["A"], prop["I"], prop["I"], prop["J"], mat["E"], mat["G"]; As_y=As_bar_f, As_z=As_bar_f)
        stresses[eid] = abs(forces["axial"] / prop["A"])
        push!(results_json["forces"]["cbar"], Dict("eid" => eid, "axial" => forces["axial"], "shear_1" => forces["shear_1"], "shear_2" => forces["shear_2"], "torque" => forces["torque"], "moment_a1" => forces["moment_a1"], "moment_a2" => forces["moment_a2"], "moment_b1" => forces["moment_b1"], "moment_b2" => forces["moment_b2"]))
    end

    return u_global, stresses, results_json
end

end
