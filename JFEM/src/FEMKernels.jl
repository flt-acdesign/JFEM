module FEM

using LinearAlgebra
using Statistics

function stiffness_frame3d(L, A, Iy, Iz, J, E, G; As_y=Inf, As_z=Inf)
    k = zeros(12, 12)
    if L < 1e-9; return k; end

    X = E * A / L
    k[1,1] = X;  k[1,7] = -X; k[7,1] = -X; k[7,7] = X

    T = G * J / L
    k[4,4] = T;  k[4,10] = -T; k[10,4] = -T; k[10,10] = T

    # Timoshenko shear parameters (Phi=0 reduces to Euler-Bernoulli)
    Phi_y = isinf(As_y) ? 0.0 : 12*E*Iy/(G*As_y*L^2)
    Phi_z = isinf(As_z) ? 0.0 : 12*E*Iz/(G*As_z*L^2)

    # Bending in xz-plane (uses Iy, shear via As_y)
    a_y = 12*E*Iy / (L^3*(1+Phi_y))
    b_y = 6*E*Iy / (L^2*(1+Phi_y))
    c_y = (4+Phi_y)*E*Iy / (L*(1+Phi_y))
    d_y = (2-Phi_y)*E*Iy / (L*(1+Phi_y))
    k[3,3] = a_y;  k[3,9] = -a_y; k[9,3] = -a_y; k[9,9] = a_y
    k[3,5] = -b_y; k[3,11] = -b_y; k[5,3] = -b_y; k[11,3] = -b_y
    k[9,5] = b_y;  k[9,11] = b_y;  k[5,9] = b_y;  k[11,9] = b_y
    k[5,5] = c_y;  k[5,11] = d_y;  k[11,5] = d_y;  k[11,11] = c_y

    # Bending in xy-plane (uses Iz, shear via As_z)
    a_z = 12*E*Iz / (L^3*(1+Phi_z))
    b_z = 6*E*Iz / (L^2*(1+Phi_z))
    c_z = (4+Phi_z)*E*Iz / (L*(1+Phi_z))
    d_z = (2-Phi_z)*E*Iz / (L*(1+Phi_z))
    k[2,2] = a_z;  k[2,8] = -a_z; k[8,2] = -a_z; k[8,8] = a_z
    k[2,6] = b_z;  k[2,12] = b_z;  k[6,2] = b_z;  k[12,2] = b_z
    k[8,6] = -b_z; k[8,12] = -b_z; k[6,8] = -b_z; k[12,8] = -b_z
    k[6,6] = c_z;  k[6,12] = d_z;  k[12,6] = d_z;  k[12,12] = c_z

    return k
end

function forces_frame3d(u_elem, L, A, Iy, Iz, J, E, G; As_y=Inf, As_z=Inf)
    k = stiffness_frame3d(L, A, Iy, Iz, J, E, G; As_y=As_y, As_z=As_z)
    f_local = k * u_elem
    return Dict(
        "axial"     => f_local[7],
        "shear_1"   => -f_local[2],
        "shear_2"   => -f_local[3],
        "torque"    => -f_local[4],
        "moment_a1" => -f_local[6],
        "moment_a2" => f_local[5],
        "moment_b1" => f_local[12],
        "moment_b2" => -f_local[11]
    )
end

function stress_frame3d(u_elem, L, A, E)
    return E * (u_elem[7] - u_elem[1]) / L
end

function shape_derivs_quad(xi, eta)
    dN_dxi  = 0.25 .* [-(1-eta),  (1-eta),  (1+eta), -(1+eta)]
    dN_deta = 0.25 .* [-(1-xi),  -(1+xi),   (1+xi),   (1-xi)]
    return dN_dxi, dN_deta
end

function stiffness_quad4(coords, E, nu, h; bend_ratio=1.0, ts_t=5.0/6.0)
    const_mem = E * h / (1 - nu^2)
    Cm = const_mem .* [1 nu 0; nu 1 0; 0 0 (1-nu)/2]

    const_bend = bend_ratio * (E * h^3) / (12 * (1 - nu^2))
    Cb = const_bend .* [1 nu 0; nu 1 0; 0 0 (1-nu)/2]

    G = E / (2*(1+nu))
    k_shear = ts_t * G * h
    Cs = k_shear .* [1 0; 0 1]

    elem_area = abs(0.5 * ((coords[1,1]*(coords[2,2]-coords[4,2]) + coords[2,1]*(coords[3,2]-coords[1,2]) + coords[3,1]*(coords[4,2]-coords[2,2]) + coords[4,1]*(coords[1,2]-coords[3,2]))))

    Ke = zeros(24, 24)
    K_ab = zeros(24, 4)  # coupling standard-incompatible (membrane only)
    K_bb = zeros(4, 4)   # incompatible-incompatible

    pt = 1.0/sqrt(3.0)
    gauss_2x2 = [-pt -pt; pt -pt; pt pt; -pt pt]

    for i in 1:4
        r, s = gauss_2x2[i,1], gauss_2x2[i,2]
        dNr, dNs = shape_derivs_quad(r, s)
        J = [dNr'; dNs'] * coords
        detJ = det(J); abs_detJ = abs(detJ)
        if abs_detJ < 1e-12; abs_detJ = 1e-12; end
        invJ = inv(J); dN_dxy = invJ * [dNr'; dNs']

        Bm = zeros(3, 24); Bb = zeros(3, 24)
        for k in 1:4
            idx = (k-1)*6
            Bm[1, idx+1] = dN_dxy[1,k]; Bm[2, idx+2] = dN_dxy[2,k]
            Bm[3, idx+1] = dN_dxy[2,k]; Bm[3, idx+2] = dN_dxy[1,k]
            Bb[1, idx+5] = dN_dxy[1,k]
            Bb[2, idx+4] = -dN_dxy[2,k]
            Bb[3, idx+5] = dN_dxy[2,k]
            Bb[3, idx+4] = -dN_dxy[1,k]
        end

        # Incompatible mode B-matrix: derivatives of bubble functions
        # phi1 = 1-xi^2, phi2 = 1-eta^2
        dphi1_dxi = -2.0*r; dphi1_deta = 0.0
        dphi2_dxi = 0.0;    dphi2_deta = -2.0*s
        iJ = inv(J)
        dphi1_dx = iJ[1,1]*dphi1_dxi + iJ[1,2]*dphi1_deta
        dphi1_dy = iJ[2,1]*dphi1_dxi + iJ[2,2]*dphi1_deta
        dphi2_dx = iJ[1,1]*dphi2_dxi + iJ[1,2]*dphi2_deta
        dphi2_dy = iJ[2,1]*dphi2_dxi + iJ[2,2]*dphi2_deta

        Bi = zeros(3, 4)
        Bi[1, 1] = dphi1_dx;  Bi[3, 1] = dphi1_dy
        Bi[2, 2] = dphi1_dy;  Bi[3, 2] = dphi1_dx
        Bi[1, 3] = dphi2_dx;  Bi[3, 3] = dphi2_dy
        Bi[2, 4] = dphi2_dy;  Bi[3, 4] = dphi2_dx

        Ke .+= (Bm' * Cm * Bm .+ Bb' * Cb * Bb) .* abs_detJ
        K_ab .+= (Bm' * Cm * Bi) .* abs_detJ
        K_bb .+= (Bi' * Cm * Bi) .* abs_detJ
    end

    # Static condensation: K_eff = K_aa - K_ab * K_bb^{-1} * K_ab'
    Ke .-= K_ab * (K_bb \ K_ab')

    # MITC4 transverse shear in covariant (natural) coordinates
    # Tying points: edge midpoints
    # e_xiz sampled at A(0,-1), C(0,+1) -> interpolated linearly in eta
    # e_etaz sampled at B(-1,0), D(+1,0) -> interpolated linearly in xi

    # Tying point A: (0,-1)
    dNr_A, dNs_A = shape_derivs_quad(0.0, -1.0)
    J_A = [dNr_A'; dNs_A'] * coords
    N_A = 0.25 .* [(1.0-0.0)*(1.0-(-1.0)), (1.0+0.0)*(1.0-(-1.0)), (1.0+0.0)*(1.0+(-1.0)), (1.0-0.0)*(1.0+(-1.0))]
    B_exi_A = zeros(1, 24)
    for k in 1:4
        idx = (k-1)*6
        B_exi_A[1, idx+3] = dNr_A[k]
        B_exi_A[1, idx+5] = N_A[k] * J_A[1,1]
        B_exi_A[1, idx+4] = -N_A[k] * J_A[1,2]
    end

    # Tying point C: (0,+1)
    dNr_C, dNs_C = shape_derivs_quad(0.0, 1.0)
    J_C = [dNr_C'; dNs_C'] * coords
    N_C = 0.25 .* [(1.0-0.0)*(1.0-1.0), (1.0+0.0)*(1.0-1.0), (1.0+0.0)*(1.0+1.0), (1.0-0.0)*(1.0+1.0)]
    B_exi_C = zeros(1, 24)
    for k in 1:4
        idx = (k-1)*6
        B_exi_C[1, idx+3] = dNr_C[k]
        B_exi_C[1, idx+5] = N_C[k] * J_C[1,1]
        B_exi_C[1, idx+4] = -N_C[k] * J_C[1,2]
    end

    # Tying point B: (-1,0)
    dNr_B, dNs_B = shape_derivs_quad(-1.0, 0.0)
    J_B = [dNr_B'; dNs_B'] * coords
    N_B = 0.25 .* [(1.0-(-1.0))*(1.0-0.0), (1.0+(-1.0))*(1.0-0.0), (1.0+(-1.0))*(1.0+0.0), (1.0-(-1.0))*(1.0+0.0)]
    B_eta_B = zeros(1, 24)
    for k in 1:4
        idx = (k-1)*6
        B_eta_B[1, idx+3] = dNs_B[k]
        B_eta_B[1, idx+5] = N_B[k] * J_B[2,1]
        B_eta_B[1, idx+4] = -N_B[k] * J_B[2,2]
    end

    # Tying point D: (+1,0)
    dNr_D, dNs_D = shape_derivs_quad(1.0, 0.0)
    J_D = [dNr_D'; dNs_D'] * coords
    N_D = 0.25 .* [(1.0-1.0)*(1.0-0.0), (1.0+1.0)*(1.0-0.0), (1.0+1.0)*(1.0+0.0), (1.0-1.0)*(1.0+0.0)]
    B_eta_D = zeros(1, 24)
    for k in 1:4
        idx = (k-1)*6
        B_eta_D[1, idx+3] = dNs_D[k]
        B_eta_D[1, idx+5] = N_D[k] * J_D[2,1]
        B_eta_D[1, idx+4] = -N_D[k] * J_D[2,2]
    end

    # Gauss integration with interpolated covariant shear strains
    for i in 1:4
        r, s = gauss_2x2[i,1], gauss_2x2[i,2]
        dNr_gp, dNs_gp = shape_derivs_quad(r, s)
        J_gp = [dNr_gp'; dNs_gp'] * coords
        detJ_gp = det(J_gp); abs_detJ_gp = abs(detJ_gp)
        if abs_detJ_gp < 1e-12; abs_detJ_gp = 1e-12; end
        invJ_gp = inv(J_gp)

        # Interpolated covariant shear B-matrices
        B_exi = 0.5*(1-s) .* B_exi_A .+ 0.5*(1+s) .* B_exi_C
        B_eta = 0.5*(1-r) .* B_eta_B .+ 0.5*(1+r) .* B_eta_D

        # Transform covariant to physical shear
        B_cov = vcat(B_exi, B_eta)  # 2x24
        Cs_cov = invJ_gp * Cs * invJ_gp'

        Ke .+= (B_cov' * Cs_cov * B_cov) .* abs_detJ_gp
    end

    # Drill rotation stabilization (Hughes-Brezzi approach)
    k_drill = 1e-3 * G * h * elem_area / 4.0
    for k in 1:4; idx=(k-1)*6+6; Ke[idx,idx]+=k_drill; end

    return Ke
end

function compute_principal_2d(s11, s22, s12)
    s_avg = (s11 + s22) / 2.0
    radius = sqrt(((s11 - s22) / 2.0)^2 + s12^2)
    return s_avg + radius, s_avg - radius
end

function stress_strain_quad4(coords, u_elem, E, nu, h, t_shell; bend_ratio=1.0)
    const_mem = E / (1 - nu^2)
    D_mem = const_mem .* [1 nu 0; nu 1 0; 0 0 (1-nu)/2]
    Cm = D_mem * h

    dNr, dNs = shape_derivs_quad(0.0, 0.0)
    J = [dNr'; dNs'] * coords
    invJ = inv(J); dN_dxy = invJ * [dNr'; dNs']

    Bm = zeros(3, 24); Bb = zeros(3, 24)
    for k in 1:4
        idx = (k-1)*6
        Bm[1, idx+1]=dN_dxy[1,k]; Bm[2, idx+2]=dN_dxy[2,k]
        Bm[3, idx+1]=dN_dxy[2,k]; Bm[3, idx+2]=dN_dxy[1,k]
        Bb[1, idx+5] = dN_dxy[1,k]
        Bb[2, idx+4] = -dN_dxy[2,k]
        Bb[3, idx+5] = dN_dxy[2,k]
        Bb[3, idx+4] = -dN_dxy[1,k]
    end

    # Recover incompatible mode amplitudes via static condensation
    K_ab_sr = zeros(24, 4); K_bb_sr = zeros(4, 4)
    pt = 1.0/sqrt(3.0)
    gauss_pts = [-pt -pt; pt -pt; pt pt; -pt pt]
    for i in 1:4
        r, s = gauss_pts[i,1], gauss_pts[i,2]
        dNr_g, dNs_g = shape_derivs_quad(r, s)
        J_g = [dNr_g'; dNs_g'] * coords
        detJ_g = abs(det(J_g))
        if detJ_g < 1e-12; detJ_g = 1e-12; end
        iJ = inv(J_g)
        dN_dxy_g = iJ * [dNr_g'; dNs_g']

        Bm_g = zeros(3, 24)
        for k in 1:4
            idx = (k-1)*6
            Bm_g[1, idx+1] = dN_dxy_g[1,k]; Bm_g[2, idx+2] = dN_dxy_g[2,k]
            Bm_g[3, idx+1] = dN_dxy_g[2,k]; Bm_g[3, idx+2] = dN_dxy_g[1,k]
        end

        dphi1_dx = iJ[1,1]*(-2.0*r); dphi1_dy = iJ[2,1]*(-2.0*r)
        dphi2_dx = iJ[1,2]*(-2.0*s); dphi2_dy = iJ[2,2]*(-2.0*s)

        Bi = zeros(3, 4)
        Bi[1, 1] = dphi1_dx;  Bi[3, 1] = dphi1_dy
        Bi[2, 2] = dphi1_dy;  Bi[3, 2] = dphi1_dx
        Bi[1, 3] = dphi2_dx;  Bi[3, 3] = dphi2_dy
        Bi[2, 4] = dphi2_dy;  Bi[3, 4] = dphi2_dx

        K_ab_sr .+= (Bm_g' * Cm * Bi) .* detJ_g
        K_bb_sr .+= (Bi' * Cm * Bi) .* detJ_g
    end

    alpha = -(K_bb_sr \ (K_ab_sr' * u_elem))

    # GP-averaged strains including incompatible mode contribution
    eps_mem_avg = zeros(3)
    kappa_avg = zeros(3)
    total_area = 0.0
    for i in 1:4
        r, s = gauss_pts[i,1], gauss_pts[i,2]
        dNr_g, dNs_g = shape_derivs_quad(r, s)
        J_g = [dNr_g'; dNs_g'] * coords
        detJ_g = abs(det(J_g))
        if detJ_g < 1e-12; detJ_g = 1e-12; end
        iJ = inv(J_g)
        dN_dxy_g = iJ * [dNr_g'; dNs_g']

        Bm_g = zeros(3, 24)
        for k in 1:4
            idx = (k-1)*6
            Bm_g[1, idx+1] = dN_dxy_g[1,k]; Bm_g[2, idx+2] = dN_dxy_g[2,k]
            Bm_g[3, idx+1] = dN_dxy_g[2,k]; Bm_g[3, idx+2] = dN_dxy_g[1,k]
        end

        Bb_g = zeros(3, 24)
        for k in 1:4
            idx = (k-1)*6
            Bb_g[1, idx+5] = dN_dxy_g[1,k]
            Bb_g[2, idx+4] = -dN_dxy_g[2,k]
            Bb_g[3, idx+5] = dN_dxy_g[2,k]
            Bb_g[3, idx+4] = -dN_dxy_g[1,k]
        end

        dphi1_dx = iJ[1,1]*(-2.0*r); dphi1_dy = iJ[2,1]*(-2.0*r)
        dphi2_dx = iJ[1,2]*(-2.0*s); dphi2_dy = iJ[2,2]*(-2.0*s)

        Bi = zeros(3, 4)
        Bi[1, 1] = dphi1_dx;  Bi[3, 1] = dphi1_dy
        Bi[2, 2] = dphi1_dy;  Bi[3, 2] = dphi1_dx
        Bi[1, 3] = dphi2_dx;  Bi[3, 3] = dphi2_dy
        Bi[2, 4] = dphi2_dy;  Bi[3, 4] = dphi2_dx

        eps_gp = Bm_g * u_elem .+ Bi * alpha
        kappa_gp = Bb_g * u_elem

        eps_mem_avg .+= eps_gp .* detJ_g
        kappa_avg .+= kappa_gp .* detJ_g
        total_area += detJ_g
    end
    eps_mem_avg ./= total_area
    kappa_avg ./= total_area

    eps_mem = eps_mem_avg
    kappa = kappa_avg

    N = (D_mem * eps_mem) * h
    M = -bend_ratio * (D_mem * kappa) * (h^3/12.0)
    Q = [0.0, 0.0]

    z1 = -h/2.0; z2 = h/2.0

    strain_z1 = eps_mem .+ z1 .* kappa
    stress_z1 = D_mem * strain_z1
    strain_z2 = eps_mem .+ z2 .* kappa
    stress_z2 = D_mem * strain_z2

    return N, M, Q, stress_z1, stress_z2, strain_z1, strain_z2
end

function stiffness_tria3(coords, E, nu, h)
    x, y = coords[:,1], coords[:,2]
    A = 0.5 * abs(x[1]*(y[2]-y[3]) + x[2]*(y[3]-y[1]) + x[3]*(y[1]-y[2]))
    if A < 1e-12; return zeros(18,18); end
    b = [y[2]-y[3], y[3]-y[1], y[1]-y[2]] ./ (2*A)
    c = [x[3]-x[2], x[1]-x[3], x[2]-x[1]] ./ (2*A)
    B = zeros(3, 6)
    for i in 1:3; B[1, i*2-1]=b[i]; B[2, i*2]=c[i]; B[3, i*2-1]=c[i]; B[3, i*2]=b[i]; end
    D = (E * h / (1 - nu^2)) .* [1 nu 0; nu 1 0; 0 0 (1-nu)/2]
    K_mem = B' * D * B * A
    Ke = zeros(18, 18)
    m_idx = [1,2, 7,8, 13,14]
    Ke[m_idx, m_idx] = K_mem
    G = E / (2*(1+nu))
    k_drill = 1e-3 * G * h * A / 3.0
    for i in [6,12,18]; Ke[i,i] = k_drill; end
    return Ke
end

function stress_strain_tria3(coords, u_elem, E, nu, h)
    x, y = coords[:,1], coords[:,2]
    A = 0.5 * abs(x[1]*(y[2]-y[3]) + x[2]*(y[3]-y[1]) + x[3]*(y[1]-y[2]))
    if A < 1e-12; return zeros(3), zeros(3), zeros(2), zeros(3), zeros(3), zeros(3), zeros(3); end
    b = [y[2]-y[3], y[3]-y[1], y[1]-y[2]] ./ (2*A)
    c = [x[3]-x[2], x[1]-x[3], x[2]-x[1]] ./ (2*A)
    B = zeros(3, 6)
    for i in 1:3; B[1, i*2-1]=b[i]; B[2, i*2]=c[i]; B[3, i*2-1]=c[i]; B[3, i*2]=b[i]; end
    D = (E / (1 - nu^2)) .* [1 nu 0; nu 1 0; 0 0 (1-nu)/2]
    u_mem = [u_elem[1], u_elem[2], u_elem[7], u_elem[8], u_elem[13], u_elem[14]]
    strain = B * u_mem
    stress = D * strain
    return stress*h, zeros(3), zeros(2), stress, stress, strain, strain
end

end
