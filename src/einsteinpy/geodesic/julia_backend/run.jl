"""
This script acts as a bridge between "einsteinpy.geodesic" and "KerrSolver.jl".
It solves the system, using "KerrSolver.jl", and saves the output to .csv files in 
`/tmp/epy_geod_jl_temp/` (on *nix) or `%TEMP%/epy_geod_jl_temp/` (on Windows).
These are accessed by the `_julia_wrapper` utility function in EinsteinPy.

# Stand-alone Usage Example
```
> julia run.jl 10., 1.57, 0., 0., 0., 2.427, 0.99, 0.9, 0.5, 0.05
("Success", "C:\\Users\\UserName\\AppData\\Local\\Temp\\epy_geod_jl")
```
where the argument order is q1, q2, q3, p1, p2, p3, a, E, end_lambda, step_size.
See "einsteinpy.geodesic.Geodesic" or "KerrSolver.jl" for information on the arguments.
The output is a 2-tuple of strings, containing the "return code" from the solver
and "folder-path", where the output is stored.

"""

import Base
import DelimitedFiles: writedlm
import Main: include

Main.include("KerrSolver.jl")

function run(q, p, params)
    # Solving
    retcode, lambdas, vecs = Main.KerrSolver.solveSystem(q, p, params, end_lambda, step_size)

    # Saving to file in 'temp' directory
    temp_dir_name = "epy_geod_jl_temp/"
    cd(tempdir())
    mkpath(temp_dir_name) # Does nothing, if folder exists
    cd(temp_dir_name)

    # Writing to csvs (Overwrites existing files)
    writedlm("lambdas.csv", lambdas, ",")
    writedlm("vecs.csv", vecs, ",")

    # Return retcode and folder-path as a 2-tuple of strings
    return retcode, pwd()
end

# Parsing input
q1, q2, q3, p1, p2, p3, a, E, end_lambda, step_size = map(x -> parse(Float64, x), Base.ARGS)

# Packaging arguments
q = [q1, q2, q3]
p = [p1, p2, p3]
params = [a, E]

# Final output
println(run(q, p, params))
