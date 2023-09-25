

"""
    plotff(F, points_sph; scale="log", normalize=true, type="abs")

Plot the 3D far-field pattern in F at the given spherical points (lying on a sphercial grid).

Optional arguments:
 - scale: "log" or "linear"
 - normalization: true or false
 - type: "abs" or "theta" or "phi"
"""
function SphericalScattering.plotff(F, points_sph; scale="log", normalize=true, type="abs")

    traces = PlotlyJS.GenericTrace[]
    text = ""

    if type == "theta"
        Fsph = SphericalScattering.convertCartesian2Spherical.(F, points_sph)
        FF = abs.([Fsph[i][2] for (i, j) in enumerate(Fsph)])
    elseif type == "phi"
        Fsph = SphericalScattering.convertCartesian2Spherical.(F, points_sph)
        FF = abs.([Fsph[i][3] for (i, j) in enumerate(Fsph)])
    else
        FF = norm.(F)
    end

    if scale == "log"

        text = " in dB"
        FF[FF .< 1e-3] .= 1e-3 # avoid -Inf
        FFdata = 20 * log10.(FF / minimum(FF))      # shift to positive values only

        if normalize
            FFscale = 20 * log10.(FF / maximum(FF))
        else
            FFscale = 20 * log10.(FF)
        end

    else
        text = ""
        FFdata = FF

        if normalize
            FFscale = FF /= maximum(FF)
        else
            FFscale = FF
        end
    end

    FFdata /= maximum(FFdata)

    # --- convert data: surface with (x,y,z) where z is determined by the FF
    data = [[FFdata[ind], points_sph[ind][2], points_sph[ind][3]] for ind in eachindex(FF)]
    data = SphericalScattering.sph2cart.(data)

    # --- bring it in fitting format for Plotly
    x = reshape([data[ind][1] for ind in eachindex(data)], size(points_sph))
    y = reshape([data[ind][2] for ind in eachindex(data)], size(points_sph))
    z = reshape([data[ind][3] for ind in eachindex(data)], size(points_sph))

    t = PlotlyJS.surface(;
        x=x, y=y, z=z, surfacecolor=FFscale, colorscale="Viridis", colorbar=PlotlyJS.PlotlyBase.attr(; title="FF$text")
    )
    push!(traces, t)

    maxmax = maximum([maximum(abs.(x)), maximum(abs.(y)), maximum(abs.(z))])

    # ensure all three axes have the same length
    layout = PlotlyJS.Layout(;
        scene=PlotlyJS.attr(;
            xaxis=PlotlyJS.attr(; visible=true, legend=:none, range=[-maxmax, maxmax]),
            yaxis=PlotlyJS.attr(; visible=true, range=[-maxmax, maxmax]),
            zaxis=PlotlyJS.attr(; visible=true, range=[-maxmax, maxmax]),
            aspectratio=PlotlyJS.attr(; x=1, y=1, z=1),
        ),
    )

    # --- plot it
    PlotlyJS.plot(traces, layout)
end


"""
    plotffcut(F, points; scale="log", normalize=true, format="polar")

Plot a cut of the far-field in a rectangular or a polar plot.

The format can be "polar" or "rectangular".
"""
function SphericalScattering.plotffcut(F, points; scale="log", normalize=true, format="polar")

    traces = PlotlyJS.GenericTrace[]

    if scale == "log"

        if format == "polar"

            FF = F
            FF[FF .< 1e-3] .= 1e-3 # avoid -Inf
            FF = 20 * log10.(FF / minimum(FF))
        else
            FF = 20 * log10.(F)
        end

        if normalize
            FF .-= maximum(FF)
        end
    else
        FF = F
        if normalize
            FF /= maximum(FF)
        end
    end

    if format == "polar"
        t = scatterpolar(; r=FF, theta=rad2deg.(points), mode="lines+markers", marker=attr(; size=5, sizeref=0.03))
    else
        t = PlotlyJS.scatter(; x=rad2deg.(points), y=FF, mode="lines+markers", marker=attr(; size=5, sizeref=0.03))
    end
    push!(traces, t)

    # --- plot it
    PlotlyJS.plot(traces)
end
