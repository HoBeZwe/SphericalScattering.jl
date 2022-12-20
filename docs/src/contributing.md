
# Contributing

In order to contribute to this package directly create a pull request against the `master` branch. Before doing so please: 

- Follow the style of the surrounding code.
- Supplement the documentation.
- Write tests and check that no errors occur.


---
## Style

For a consistent style the [JuliaFormatter.jl](https://github.com/domluna/JuliaFormatter.jl) package is used which enforces the style defined in the *.JuliaFormatter.toml* file. To follow this style simply run
```julia
using JuliaFormatter
format(pkgdir(SphericalScattering))
```

!!! note
    That all files follow the JuliaFormatter style is tested during the unit tests. Hence, do not forget to execute the two lines above. Otherwise, the tests are likely to not pass.


---
## Documentation

Add documentation for any changes or new features following the style of the existing documentation. For more information you can have a look at the [Documenter.jl](https://documenter.juliadocs.org/stable/) documentation.


---
## Tests

Write tests for your code changes and verify that no errors occur, e.g., by running
```julia
using Pkg
Pkg.test("SphericalScattering")
```

For more detailed information on which parts are tested the coverage can be evaluated on your local machine, e.g., by
```julia
using Pkg
Pkg.test("SphericalScattering"; coverage=true, julia_args=`--threads 4`)

# determine coverage
using Coverage
src_folder = pkgdir(SphericalScattering) * "/src"
coverage   = process_folder(src_folder)
LCOV.writefile("path-to-folder-you-like" * "SphericalScattering.lcov.info", coverage)

clean_folder(src_folder) # delete .cov files

# extract information about coverage
covered_lines, total_lines = get_summary(coverage)
@info "Current coverage:\n$covered_lines of $total_lines lines ($(round(Int, covered_lines / total_lines * 100)) %)"
```

In Visual Studio Code the [Coverage Gutters](https://marketplace.visualstudio.com/items?itemName=ryanluker.vscode-coverage-gutters) plugin can be used to visualize the tested lines of the code by inserting the path of the *SphericalScattering.lcov.info* file in the settings.