function potentialEnergyFig(dat_DF,system_DF)
"""
    Function that stores a comparisson of potential energy
"""
    y_min=round(minimum(map(s->minimum(s.c_ep),system_DF)),sigdigits=1);

    fig=Figure()
    ax_f=Axis(fig[1:1,1:1])
    ax_f2=Axis(fig[1:1,1:1])

    ax=Axis(fig[1:3,1:1],
        title=latexstring("\\mathrm{Potential~Energy}"),
        #subtitle=latexstring(subtitle),
        xlabel=L"\mathrm{Time~units}~[\tau^*]",
        ylabel=L"U~[\mathrm{J}/\epsilon]",
        xminorticksvisible=true,
        xminorgridvisible=true,
        limits=(0,nothing,y_min,nothing),
        #xscale=log10,
        #yscale=log10
    )

    hidespines!(ax_f)
    hidedecorations!(ax_f)
    hidespines!(ax_f2)
    hidedecorations!(ax_f2)

    linkyaxes!(ax,ax_f,ax_f2)
    linkxaxes!(ax,ax_f,ax_f2)

    for it in 1:nrow(dat_DF)
        plot!(ax,dat_DF."time-step"[it].*system_DF[it].TimeStep,system_DF[it].c_ep,
            label=latexstring(100*dat_DF.phi[it]),
            )
        plot!(ax_f,[-1],[0],
            label=latexstring(100*dat_DF."CL-Con"[it])
            )
        plot!(ax_f2,[-1],[0],
            label=latexstring(dat_DF."Temperature"[it])
            )
    end

    Legend(fig[1,2],ax,L"\phi~\%")
    Legend(fig[2,2],ax_f,L"\mathrm{CL}~\%",merge=true)
    Legend(fig[3,2],ax_f2,L"\mathrm{T}",merge=true)

    save(joinpath(pwd(),"ep.png"), fig, px_per_unit = 300/inch)
end


