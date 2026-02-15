"""
    Script with functions for the three body potential
"""

function Upatch(eps,Dp,r)
"""
    Interaction between patches

    Input parameters
    eps ---- Energy of the interaction
    Dp ----- Patch diameter
    r ------ Distance between patches

    Derived parameters:
    r_c ---- Cut-off distance
"""
    r_c=1.5*Dp;
    if r<= r_c
        return 2*eps*( (1/2)*(Dp/r)^4 - 1 )*exp( (Dp/(r-r_c)) + 2 )
    else
        return 0.0
    end
end

function U3(eps_jk,eps,Dp,r)
"""
    Potencial U3, which is an auxiliary function for the swap potential

    Input parameters
    eps_jk -- Interaction energy between particle j and particle k
    eps ----- Interaction energy between particle i and particle j/k
    Dp ------ Patch diameter
    r ------- Distance between patches

    Derived parameters
    r_min --- Distance such that Upatch(eps,Dp,r_min) = -eps
    r_c ----- Cut-off distance of Upatch
"""
    r_min=Dp;
    r_c=1.5*Dp;
    if r>=0 && r<=r_min
        return 1.0
    elseif r>r_min && r<r_c
        return -Upatch(eps,Dp,r)/eps_jk
    else
        return 0.0
    end
end

function Uswap(w,eps_jk,eps_ij,eps_ik,Dp,r_ij,r_ik)
"""
    Swap potential of 3 body interactions

    Input parameters
    w ------- Interpolate the limits of swapping (w=1 swap, w>>1 no swap)
    eps_jk -- Interaction energy between particle j and particle k
    eps_ik -- Interaction energy between particle i and particle k
    eps_ij -- Interaction energy between particle i and particle j
    Dp ------ Patch diameter
    r_ij ---- Distance between patch i and patch j
    r_ik ---- Distance between patch i and patch k
"""
    return w*eps_jk*U3(eps_jk,eps_ij,Dp,r_ij)*U3(eps_jk,eps_ik,Dp,r_ik)
end

function forcePatch(eps,Dp,r)
"""
    Magnitude and components of the force given by the patch-patch interaction potential

    Input parameters
    eps ---- Energy of the interaction
    Dp ----- Patch diameter
    r ------ Distance between patches

    Derived parameters:
    r_c ---- Cut-off distance

    Output
    comp --- Component of the force
    map ---- Magnitude of the force

"""
    r_c=1.5*Dp;
    comp=(1/r - 2*(eps/r)*exp( Dp/(r-r_c) + 2 ) - ( Dp/(r-r_c)^2 ) )*Upatch(eps,Dp,r)
    mag=sqrt(comp^2)
    return (comp,mag) 
end

function forceSwap(w,eps_jk,eps_ij,eps_ik,Dp,r_ij,r_ik)
"""
    Magnitude and components of the force given by the swap 3 body interaction potential

    Input parameters
    w ------- Interpolate the limits of swapping (w=1 swap, w>>1 no swap)
    eps_jk -- Interaction energy between particle j and particle k
    eps_ik -- Interaction energy between particle i and particle k
    eps_ij -- Interaction energy between particle i and particle j
    Dp ------ Patch diameter
    r_ij ---- Distance between patch i and patch j
    r_ik ---- Distance between patch i and patch k

    Output
    comp_ij - ij component of the force
    comp_ik - ik component of the force
    mag ----- The norm of the force
"""
    comp_ij=w*eps_jk*U3(eps_jk,eps_ik,Dp,r_ik)*first(forcePatch(eps_ij,Dp,r_ij))
    comp_ik=w*eps_jk*U3(eps_jk,eps_ij,Dp,r_ij)*first(forcePatch(eps_ik,Dp,r_ik))
    mag=sqrt( comp_ij^2 + comp_ik^2 )
    return [comp_ij comp_ik mag]
end
