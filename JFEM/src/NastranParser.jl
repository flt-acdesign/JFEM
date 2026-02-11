module NastranParser

using LinearAlgebra

function safe_get(arr::Vector{Any}, idx::Int, default_val=nothing)
    if idx > length(arr); return default_val; end
    return arr[idx]
end

function parse_nastran_number(field::Any, default_val=nothing)
    if isa(field, Number); return field; end
    if isnothing(field) || strip(string(field)) == ""; return default_val; end
    field_str = String(strip(string(field)))
    clean_field = replace(field_str, r"(\d)([+-])" => s"\1e\2")
    clean_field = replace(clean_field, "ee" => "e")
    try
        val = parse(Float64, clean_field)
        if abs(val - round(val)) < 1e-8; return Int(round(val)); end
        return val
    catch; return default_val; end
end

function to_id(val)
    if isa(val, Integer); return val; end
    if isa(val, AbstractFloat); return Int(round(val)); end
    return 0
end

function expand_nastran_list(raw_fields)
    result = Int[]
    i = 1
    while i <= length(raw_fields)
        val = raw_fields[i]
        if isa(val, AbstractString) && uppercase(strip(val)) == "THRU"
            if isempty(result) || i == length(raw_fields)
                i += 1; continue
            end
            start_id = result[end]
            end_id = to_id(parse_nastran_number(raw_fields[i+1]))
            if end_id > start_id; append!(result, (start_id+1):end_id); end
            i += 2
        else
            parsed = to_id(parse_nastran_number(val, 0))
            if parsed > 0; push!(result, parsed); end
            i += 1
        end
    end
    return result
end

function get_nastran_card_name(line::AbstractString)
    if occursin(",", line)
        parts = split(line, ",")
        return uppercase(strip(parts[1]))
    end
    return uppercase(strip(line[1:min(8, length(line))]))
end

function get_nastran_fields_from_line(line::AbstractString)
    fields = []
    if occursin(",", line)
        parts = split(line, ",")
        if length(parts) > 1
            for p in parts[2:end]; push!(fields, strip(string(p))); end
        end
    else
        padded = rpad(line, 80, ' ')
        for i in 0:7
            s = 9 + i*8; e = s + 7
            if s > length(padded); break; end
            push!(fields, strip(padded[s:min(e, end)]))
        end
    end
    return fields
end

function process_cards(lines)
    processed = Dict{String, Vector{Any}}()
    i = 1
    while i <= length(lines)
        line = lines[i]
        clean_line = strip(line)
        if startswith(clean_line, '$') || isempty(clean_line)
            i += 1; continue
        end
        name = get_nastran_card_name(line)
        if startswith(name, "+") || startswith(name, "*")
            i += 1; continue
        end
        if !isnothing(name) && !isempty(name)
            fields = Any["SMALL"; String(name)]
            append!(fields, get_nastran_fields_from_line(line))
            steps = 1
            while i + steps <= length(lines)
                next_line = lines[i+steps]
                next_clean = strip(next_line)
                if startswith(next_clean, '$')
                    steps += 1; continue
                end
                is_cont = false
                if startswith(next_line, " ") || startswith(next_clean, "+") || startswith(next_clean, "*")
                    is_cont = true
                elseif occursin(",", next_line) && startswith(next_clean, ",")
                    is_cont = true
                end
                if is_cont
                    append!(fields, get_nastran_fields_from_line(next_line))
                    steps += 1
                else
                    break
                end
            end
            if !haskey(processed, name); processed[name] = []; end
            push!(processed[name], fields)
            i += steps
        else
            i += 1
        end
    end
    return processed
end

function extract_grid(cards)
    d = Dict()
    for c in cards
        id = to_id(parse_nastran_number(safe_get(c, 3)))
        cp = to_id(parse_nastran_number(safe_get(c, 4), 0))
        x = [parse_nastran_number(safe_get(c, 5), 0.0),
             parse_nastran_number(safe_get(c, 6), 0.0),
             parse_nastran_number(safe_get(c, 7), 0.0)]
        cd = to_id(parse_nastran_number(safe_get(c, 8), 0))
        if id > 0; d[string(id)] = Dict("ID"=>id, "CP"=>cp, "CD"=>cd, "X"=>x); end
    end
    return d
end

function extract_coords(cards)
    d = Dict()
    for c in cards
        id = to_id(parse_nastran_number(safe_get(c, 3)))
        A = [parse_nastran_number(safe_get(c, 5),0.0), parse_nastran_number(safe_get(c, 6),0.0), parse_nastran_number(safe_get(c, 7),0.0)]
        B = [parse_nastran_number(safe_get(c, 8),0.0), parse_nastran_number(safe_get(c, 9),0.0), parse_nastran_number(safe_get(c, 10),0.0)]
        C = [parse_nastran_number(safe_get(c, 11),0.0), parse_nastran_number(safe_get(c, 12),0.0), parse_nastran_number(safe_get(c, 13),0.0)]
        w = B - A
        if norm(w) < 1e-9; w=[0.0,0.0,1.0]; else; w=normalize(w); end
        v_t = C - A
        v = cross(w, v_t)
        if norm(v) < 1e-9; v=[0.0,1.0,0.0]; else; v=normalize(v); end
        u = normalize(cross(v, w))
        d[string(id)] = Dict("Origin"=>A, "U"=>u, "V"=>v, "W"=>w)
    end
    return d
end

function extract_shells(cards)
    d = Dict()
    for c in cards
        id = to_id(parse_nastran_number(safe_get(c, 3)))
        pid = to_id(parse_nastran_number(safe_get(c, 4)))
        nodes = []
        for k in 5:8
            val = to_id(parse_nastran_number(safe_get(c, k), 0))
            if val > 0; push!(nodes, val); end
        end
        if id > 0; d[string(id)] = Dict("ID"=>id, "PID"=>pid, "NODES"=>nodes); end
    end
    return d
end

function extract_cbar(cards)
    d = Dict()
    for c in cards
        id = to_id(parse_nastran_number(safe_get(c, 3)))
        pid = to_id(parse_nastran_number(safe_get(c, 4)))
        ga = to_id(parse_nastran_number(safe_get(c, 5)))
        gb = to_id(parse_nastran_number(safe_get(c, 6)))
        v = [parse_nastran_number(safe_get(c, 7),0.0), parse_nastran_number(safe_get(c, 8),0.0), parse_nastran_number(safe_get(c, 9),0.0)]
        if norm(v) < 1e-6; v = [0.0, 0.0, 1.0]; end
        d[string(id)] = Dict("ID"=>id, "PID"=>pid, "GA"=>ga, "GB"=>gb, "V"=>v, "TYPE"=>"CBAR")
    end
    return d
end

function extract_rbe3(cards)
    d = Dict()
    for c in cards
        eid = to_id(parse_nastran_number(safe_get(c, 3)))
        refgrid = to_id(parse_nastran_number(safe_get(c, 5)))
        refc = to_id(parse_nastran_number(safe_get(c, 6)))
        wt = parse_nastran_number(safe_get(c, 7), 1.0)
        comps = to_id(parse_nastran_number(safe_get(c, 8)))
        raw_grids = c[9:end]
        dep_grids = expand_nastran_list(raw_grids)
        d[string(eid)] = Dict("ID"=>eid, "REFGRID"=>refgrid, "REFC"=>refc, "WT"=>wt, "COMPS"=>comps, "DEP_GRIDS"=>dep_grids)
    end
    return d
end

function extract_props_shell(cards)
    d = Dict()
    for c in cards
        pid = to_id(parse_nastran_number(safe_get(c, 3)))
        mid = to_id(parse_nastran_number(safe_get(c, 4)))
        t = parse_nastran_number(safe_get(c, 5), 0.0)
        bend_ratio = parse_nastran_number(safe_get(c, 7), 1.0)   # 12I/T^3, default=1.0
        ts_t = parse_nastran_number(safe_get(c, 9), 5.0/6.0)     # TS/T, default=5/6
        d[string(pid)] = Dict("PID"=>pid, "MID"=>mid, "T"=>t, "BEND_RATIO"=>bend_ratio, "TS_T"=>ts_t)
    end
    return d
end

function extract_pbarl(cards)
    d = Dict()
    for c in cards
        pid = to_id(parse_nastran_number(safe_get(c, 3)))
        mid = to_id(parse_nastran_number(safe_get(c, 4)))
        type = "ROD"
        dim_start_idx = 7
        for k in 6:min(12, length(c))
            val = c[k]
            if isa(val, AbstractString) && length(val) > 1
                type = strip(val); dim_start_idx = k + 1; break
            end
        end
        dim1 = 0.0
        for k in dim_start_idx:length(c)
            val = parse_nastran_number(safe_get(c, k), nothing)
            if isa(val, Number) && val > 0
                dim1 = Float64(val); break
            end
        end
        A, I, J = 1.0, 1.0, 1.0
        K_shear = 0.0
        if type == "ROD"
            R = dim1; A = pi*R^2; I = pi*R^4/4; J = pi*R^4/2
            K_shear = 6.0/7.0  # placeholder; actual kappa depends on nu, computed in solver
        end
        d[string(pid)] = Dict("PID"=>pid, "MID"=>mid, "A"=>A, "I"=>I, "J"=>J, "TYPE"=>type, "K_SHEAR"=>K_shear)
    end
    return d
end

function extract_mats(cards)
    d = Dict()
    for c in cards
        mid = to_id(parse_nastran_number(safe_get(c, 3)))
        E = parse_nastran_number(safe_get(c, 4), 0.0)
        nu = parse_nastran_number(safe_get(c, 6), 0.3)
        if E > 0; G = E/(2*(1+nu)); else; G=0.0; end
        d[string(mid)] = Dict("MID"=>mid, "E"=>E, "G"=>G, "NU"=>nu)
    end
    return d
end

function extract_loads(cards)
    f = []
    for c in cards
        sid = to_id(parse_nastran_number(safe_get(c, 3)))
        gid = to_id(parse_nastran_number(safe_get(c, 4)))
        cid = to_id(parse_nastran_number(safe_get(c, 5), 0))
        mag = parse_nastran_number(safe_get(c, 6), 0.0)
        dir = [parse_nastran_number(safe_get(c, 7),0.0), parse_nastran_number(safe_get(c, 8),0.0), parse_nastran_number(safe_get(c, 9),0.0)]
        push!(f, Dict("TYPE"=>"FORCE", "SID"=>sid, "GID"=>gid, "CID"=>cid, "Mag"=>mag, "Dir"=>dir))
    end
    return f
end

function extract_moments(cards)
    m = []
    for c in cards
        sid = to_id(parse_nastran_number(safe_get(c, 3)))
        gid = to_id(parse_nastran_number(safe_get(c, 4)))
        cid = to_id(parse_nastran_number(safe_get(c, 5), 0))
        mag = parse_nastran_number(safe_get(c, 6), 0.0)
        dir = [parse_nastran_number(safe_get(c, 7),0.0), parse_nastran_number(safe_get(c, 8),0.0), parse_nastran_number(safe_get(c, 9),0.0)]
        push!(m, Dict("TYPE"=>"MOMENT", "SID"=>sid, "GID"=>gid, "CID"=>cid, "Mag"=>mag, "Dir"=>dir))
    end
    return m
end

function extract_pload4(cards)
    p = []
    for c in cards
        sid = to_id(parse_nastran_number(safe_get(c, 3)))
        eid = to_id(parse_nastran_number(safe_get(c, 4)))
        press = [parse_nastran_number(safe_get(c, 5), 0.0)]
        push!(p, Dict("TYPE"=>"PLOAD4", "SID"=>sid, "EID"=>eid, "P"=>press[1]))
    end
    return p
end

function extract_load_combos(cards)
    combos = []
    for c in cards
        sid = to_id(parse_nastran_number(safe_get(c, 3)))
        s = parse_nastran_number(safe_get(c, 4), 1.0)
        comps = []
        for i in 5:2:length(c)-1
            s_i = parse_nastran_number(safe_get(c, i), nothing)
            l_i = to_id(parse_nastran_number(safe_get(c, i+1), nothing))
            if !isnothing(s_i) && l_i > 0
                push!(comps, Dict("S"=>s_i, "LID"=>l_i))
            end
        end
        push!(combos, Dict("SID"=>sid, "S"=>s, "COMPS"=>comps))
    end
    return combos
end

function extract_spc1(cards)
    spcs = []
    for c in cards
        sid = to_id(parse_nastran_number(safe_get(c, 3)))
        comp = string(parse_nastran_number(safe_get(c, 4), ""))
        raw_nodes = c[5:end]
        nodes = expand_nastran_list(raw_nodes)
        push!(spcs, Dict("SID"=>sid, "C"=>comp, "NODES"=>nodes))
    end
    return spcs
end

function extract_spcadd(cards)
    d = Dict()
    for c in cards
        sid = to_id(parse_nastran_number(safe_get(c, 3)))
        raw_sets = c[4:end]
        sets = expand_nastran_list(raw_sets)
        d[sid] = sets
    end
    return d
end

function read_bulk_and_case(lines::Vector{String})
    case_control = Dict("SUBCASES" => Dict{Int, Dict{String, Any}}())
    bulk_lines = String[]
    in_bulk = false
    global_load, global_spc = nothing, nothing
    current_sub = 0
    has_grav = false
    for line in lines
        cl = uppercase(split(line, '$')[1])
        if occursin("BEGIN BULK", cl); in_bulk=true; continue; end
        if occursin("ENDDATA", cl); break; end
        if !in_bulk
            if startswith(strip(cl), "SUBCASE")
                parts = split(cl)
                if length(parts) >= 2
                    val = try parse(Int, parts[2]) catch; 0 end
                    if val > 0
                        current_sub = val
                        case_control["SUBCASES"][current_sub] = Dict{String, Any}("LOAD"=>global_load, "SPC"=>global_spc)
                    end
                end
            elseif occursin("=", cl)
                k, v = strip.(split(cl, "="))
                val = try parse(Int, v) catch; v end
                if current_sub > 0
                    case_control["SUBCASES"][current_sub][k] = val
                else
                    if k == "LOAD"; global_load = val; end
                    if k == "SPC"; global_spc = val; end
                end
            end
        else
            if length(strip(cl)) > 1
                if startswith(strip(cl), "GRAV"); has_grav = true; end
                push!(bulk_lines, String(rstrip(cl)))
            end
        end
    end
    if has_grav; println("[WARNING] GRAV detected but not supported."); end
    return case_control, bulk_lines
end

end
