

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


# cyclic colour scale for the phase (Plotly.js provides none)
const phaseColorscale = [
    [0.0, "#ff0000"],
    [1 / 6, "#ffff00"],
    [2 / 6, "#00ff00"],
    [3 / 6, "#00ffff"],
    [4 / 6, "#0000ff"],
    [5 / 6, "#ff00ff"],
    [1.0, "#ff0000"],
]


# --- plot the `values` at the `points` of `surfacePoints`, documented with the stub of `plotsurface`
function SphericalScattering.plotsurface(
    sphere::SphericalScattering.Scatterer, values, points; type="abs", scale="linear", normalize=false, clip=1.0, gap=0.02, title=""
)

    size(values) == size(points) || error("The values do not match the points: $(size(values)) and $(size(points)).")
    type in ("abs", "real", "imag", "phase") || error("The type has to be \"abs\", \"real\", \"imag\", or \"phase\".")
    scale in ("linear", "log") || error("The scale has to be \"linear\" or \"log\".")
    scale == "log" && type != "abs" && error("A logarithmic scale is available for the type \"abs\" only.")
    0 < clip <= 1 || error("The quantile `clip` has to lie in (0, 1].")

    data, colorscale, cmin, cmax, unit = surfaceData(values; type, scale, normalize, clip)

    traces = PlotlyJS.GenericTrace[]
    extent = 0.0

    for (i, (rows, offset)) in enumerate(surfaceFaces(sphere, points, gap))

        x = [p[1] + offset[1] for p in points[rows, :]]
        y = [p[2] + offset[2] for p in points[rows, :]]
        z = [p[3] + offset[3] for p in points[rows, :]]

        t = PlotlyJS.surface(;
            x=x,
            y=y,
            z=z,
            surfacecolor=data[rows, :],
            colorscale=colorscale,
            cmin=cmin,
            cmax=cmax,
            showscale=isone(i), # one colour bar for both faces of a disc
            colorbar=PlotlyJS.attr(; title=title * unit),
        )
        push!(traces, t)

        extent = max(extent, maximum(abs, x), maximum(abs, y), maximum(abs, z))
    end

    # equal ranges of the axes in a cube keep the proportions of the scatterer; unlike scaling the axes to the data,
    # they also keep the default view suited to a slender or a flat one
    limits = [-1.05 * extent, 1.05 * extent]

    layout = PlotlyJS.Layout(;
        scene=PlotlyJS.attr(;
            xaxis=PlotlyJS.attr(; title="x", range=limits),
            yaxis=PlotlyJS.attr(; title="y", range=limits),
            zaxis=PlotlyJS.attr(; title="z", range=limits),
            aspectratio=PlotlyJS.attr(; x=1, y=1, z=1),
        ),
    )

    return PlotlyJS.plot(traces, layout)
end


# --- compute the surface density `Q` and plot it
function SphericalScattering.plotsurface(
    sphere::SphericalScattering.Scatterer{<:SphericalScattering.AcousticBoundary},
    excitation::SphericalScattering.AcousticExcitation,
    ::Type{Q};
    total=true,
    resolution=5,
    parameter=SphericalScattering.Parameter(),
    kwargs...,
) where {Q<:Union{PressureTrace,PressureNormalGradient,PressureJump}}

    points = surfacePoints(sphere; resolution=resolution)

    values = surfaceValues(sphere, excitation, Q, points, total, parameter)

    names = Dict(PressureTrace => "γ₀p", PressureNormalGradient => "γ₁p", PressureJump => "[p]")
    title = Q === PressureJump || total ? names[Q] : names[Q] * "ˢ" # the jump of the total and the scattered field agree

    return SphericalScattering.plotsurface(sphere, values, points; title=surfaceTitle(title, get(kwargs, :type, "abs")), kwargs...)
end


# --- the surface density, the normals of a spheroid being its outward normals, and its modal coefficients computed once
function surfaceValues(sphere::SphericalScattering.Scatterer, excitation, Q, points, total, parameter)

    quantity = Q === PressureNormalGradient ? PressureNormalGradient(points) : Q(points)

    return if total
        field(sphere, excitation, quantity; parameter=parameter)
    else
        scatteredfield(sphere, excitation, quantity; parameter=parameter)
    end
end

function surfaceValues(sphere::Spheroid, excitation, Q, points, total, parameter)

    quantity = Q === PressureNormalGradient ? PressureNormalGradient(sphere, points) : Q(points)

    md = SphericalScattering.modes(sphere, excitation; parameter=parameter)

    return if total
        field(sphere, excitation, md, quantity; parameter=parameter)
    else
        scatteredfield(sphere, excitation, md, quantity; parameter=parameter)
    end
end


function surfaceTitle(name, type)

    type == "abs" && return "|$name|"
    type == "real" && return "Re $name"
    type == "imag" && return "Im $name"

    return "arg $name"
end


# --- the data for the colour, its colour scale and range, and the unit appended to the title
function surfaceData(values; type, scale, normalize, clip)

    if type == "phase"
        return rad2deg.(angle.(values)), phaseColorscale, -180.0, 180.0, " in deg"
    end

    magnitude = abs.(values)
    reference = normalize ? maximum(magnitude) : one(eltype(magnitude))

    if type == "abs"

        if scale == "log"
            positive = filter(>(0), magnitude)
            floor = isempty(positive) ? one(eltype(magnitude)) : minimum(positive) # avoid -Inf
            data = 20 * log10.(max.(magnitude, floor) / reference)
            cmax = quantileOf(data, clip)
            return data, "Viridis", minimum(data), cmax, " in dB"
        end

        data = magnitude / reference
        return data, "Viridis", zero(eltype(data)), quantileOf(data, clip), ""
    end

    data = (type == "real" ? real.(values) : imag.(values)) / reference
    cmax = quantileOf(abs.(data), clip)

    return data, "RdBu", -cmax, cmax, "" # symmetric, so that zero is white
end

# the quantile `q` of the finite entries, without depending on Statistics
function quantileOf(data, q)

    finite = sort(filter(isfinite, vec(data)))
    isempty(finite) && return one(eltype(data))

    return finite[clamp(ceil(Int, q * length(finite)), 1, length(finite))]
end


# --- the rows of the points forming a face and its displacement: a disc has two faces, coinciding in space, which
#     are separated by the `gap` along the axis
surfaceFaces(sphere::SphericalScattering.Scatterer, points, gap) = [(axes(points, 1), SVector(0.0, 0.0, 0.0))]

function surfaceFaces(sphere::Spheroid, points, gap)

    isdisc(sphere) || return [(axes(points, 1), SVector(0.0, 0.0, 0.0))]

    shift = gap / 2 * equatorialRadius(sphere) * sphere.axis

    # the face follows from the side of the disc on which the points lie, see `surfacePoints`
    upper = [dot(points[row, 1], sphere.axis) > 0 for row in axes(points, 1)]

    return [(findall(upper), shift), (findall(.!upper), -shift)]
end
