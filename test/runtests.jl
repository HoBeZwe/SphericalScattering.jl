using TestItemRunner

# The test items are discovered in the files below `test/`, each running in a module of its own. The
# setup they share is provided by the snippets below, which the items request individually.

@testsnippet Setup begin

    # the snippets are evaluated before the package under test is brought into scope by the test item,
    # hence they have to import it themselves
    using SphericalScattering

    using StaticArrays
    using LinearAlgebra

    # ----- points on spherical grid
    ϑ = range(0.0; stop=π, length=18)  # 10° steps
    ϕ = range(0.0; stop=2π, length=36) # 10° steps

    function getDefaultPoints(r::Float64)
        points_cart = [SVector{3,Float64}(r * cos(φ) * sin(θ), r * sin(φ) * sin(θ), r * cos(θ)) for θ in ϑ, φ in ϕ]
        points_sph = [SVector{3,Float64}(r, θ, φ) for θ in ϑ, φ in ϕ]

        return points_cart, points_sph
    end

    # ----- variables used in all tests
    spRadius = 1.0 # radius of sphere
    f = 1e8        # frequency

    𝜇 = SphericalScattering.μ0
    𝜀 = SphericalScattering.ε0

    c = 1 / sqrt(𝜇 * 𝜀)

    points_cartFF, points_sphFF = getDefaultPoints(1.0)
    points_cartNF, points_sphNF = getDefaultPoints(5.0)
    points_cartNF_inside, ~ = getDefaultPoints(0.5)
end


@testsnippet BEASTSetup begin

    using SphericalScattering

    using BEAST
    using CompScienceMeshes

    # ----- interface to BEAST
    function (lc::Excitation)(p)

        F_cart = field(lc, ElectricField([p]))

        return F_cart[1]
    end

    BEAST.cross(::BEAST.NormalVector, p::Excitation) = CrossTraceMW(p)
    BEAST.scalartype(p::Excitation) = Complex{typeof(p.embedding.ε)}

    # ----- the mesh and the basis of the method of moments; `spRadius` stems from `Setup`
    Γ = meshsphere(spRadius, 0.45)
    RT = raviartthomas(Γ)
end


@testsnippet Impedance begin

    # the impedance matrix of the method of moments; `f`, `c` and `RT` stem from the snippets above
    κ = 2π * f / c   # Wavenumber

    𝑇 = Maxwell3D.singlelayer(; wavenumber=κ)
    T = assemble(𝑇, RT, RT)
end


@testitem "Formatting of files" begin
    using JuliaFormatter
    pkgpath = pkgdir(SphericalScattering)   # path of this package including name
    @test format(pkgpath, overwrite=false)  # check whether files are formatted according to the .JuliaFormatter.toml
end

@testitem "Method ambiguities" begin
    # the scatterers are dispatched on by geometry, by condition, and by physics; mixing these levels in the
    # signatures of one function is how ambiguities arise, see `Scatterer`
    @test isempty(Test.detect_ambiguities(SphericalScattering; recursive=true))
end


@run_package_tests verbose = true
