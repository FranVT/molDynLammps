"""
    Script with all functions for potentials and forces
"""

function forceSwapTable(w,eps_ij,eps_ik,eps_jk,sig_p,r_ij,r_ik,th)
"""
    Compute the scalars for the proyection of the forces
"""

    th = deg2rad(th);
    r_jk = sqrt(r_ij^2+r_ik^2-2*r_ij*r_ik*cos(th));

#    if r_ij<=r_c && r_ik<=r_c && r_jk<=r_c 
        f_1=forceSwap(w,eps_jk,eps_ij,eps_ik,sig_p,r_ij,r_ik);
        f_2=forceSwap(w,eps_ik,eps_ij,eps_jk,sig_p,r_ij,r_jk);
        f_3=forceSwap(w,eps_ij,eps_ik,eps_jk,sig_p,r_ik,r_jk);

        f_i1=f_1[1]-f_2[1];
        f_i2=f_1[2]-f_3[1];
   
        f_j1=-f_i1; 
        f_j2=f_2[2]-f_3[2]; 

        f_k1=-f_i2; 
        f_k2=-f_j2; 

        eng=0; #Uswap(w,eps_ij,eps_ik,eps_jk,sig_p,r_ij,r_ik) + Uswap(w,eps_ij,eps_ik,eps_jk,sig_p,r_ij,r_jk) + Uswap(w,eps_ij,eps_ik,eps_jk,sig_p,r_ik,r_jk);
        
    return (f_i1,f_i2,f_j1,f_j2,f_k1,f_k2,eng)

end

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
    if r<r_c
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
    if r>r_min && r<=r_c
        return -Upatch(eps,Dp,r)/eps_jk
    elseif r<=r_min
        return 1.0
    elseif r>r_c
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
    if r>=r_c
        comp=0;
    else
        comp=(4*Dp^4*eps)/(r^5)*exp( (Dp/(r-r_c)) + 2 ) + ((2*Dp*eps)/(r-r_c)^2)*exp( (Dp/(r-r_c)) + 2 )*( (1/2)*(Dp/r)^4 - 1 )
    end
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
