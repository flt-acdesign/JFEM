module NastranMath

using LinearAlgebra

function vec_normalize(v::Vector{Float64})
    return normalize(v)
end

function makevector_from_points(A, B)
    return B - A
end

function vector_product(v1, v2)
    return cross(v1, v2)
end

"""
    CORDR_transform(Origin, U, V, W, P_local)

Transforms a point from a local coordinate system (defined by origin and basis vectors U,V,W)
to the global coordinate system.
"""
function CORDR_transform(Origin, U, V, W, P_local)
    R = hcat(U, V, W)
    return Origin + R * P_local
end

function CORDR_transform_vector(U, V, W, Vec_local)
    R = hcat(U, V, W)
    return R * Vec_local
end

end
