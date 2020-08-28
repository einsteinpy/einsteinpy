module KerrSolver

using DifferentialEquations
using ODEInterfaceDiffEq

"""
    KerrHamiltonian(q, p, params)::Float64

This function defines the general Hamiltonian for Kerr Spacetime (rotating, uncharged black hole).

# Input
- `q::Array{Float64,1}`: Length-3 Array, containing the initial 3-Position
- `p::Array{Float64,1}`: Length-3 Array, containing the initial Covariant 3-Momentum
- `params::Array{Float64,1}`: Length-3 Array, containing - 
                                Black Hole Spin Parameter, `a`,
                                Test Particle Energy, `E`,
                                Test Particle Rest Mass, `mu`

# Output
- `Float64`: Functional evaluation of the Hamiltonian, for use with `solveSystem()`

"""
function KerrHamiltonian(q, p, params)
    a = params[1]
    E = params[2]

    return -((a^4 * (E^2 - 2 * p[1]^2) - 8 * a * E * p[3] * q[1] - 
        2 * q[1] * (p[2]^2 * (-2 + q[1]) + p[1]^2 * (-2 + q[1])^2 * q[1] - E^2 * q[1]^3) + 
        a^2 * (2 * p[3]^2 - 2 * p[2]^2 - 4 * p[1]^2 * (-2 + q[1]) * q[1] + 
        E^2 * q[1] * (2 + 3 * q[1])) + (a^2 + (-2 + q[1]) * q[1]) * (a^2 * E^2 * cos(2 * q[2]) - 
        2 * p[3]^2 * csc(q[2])^2))/(4 * (a^2 + (-2 + q[1]) * q[1]) * (q[1]^2 + 
        a^2 * cos(q[2])^2)))
end

"""
    solveSystem(state, params, end_lambda, step_size)::(String, Array{Float64,1}, Array{Float64,2})

This function solves for geodesics in Kerr Spacetime, using the Hamiltonian defined by 
`KerrHamiltonian()`. A fixed-step, symplectic integrator, `VerletLeapfrog()`, is used
to integrate the system.

# Input
- `q::Array{Float64,1}`: Length-3 Array, containing the initial 3-Position
- `p::Array{Float64,1}`: Length-3 Array, containing the initial Covariant 3-Momentum
- `params::Array{Float64,1}`: Length-3 Array, containing - 
                                Black Hole Spin Parameter, `a`,
                                Test Particle Energy, `E`,
                                Test Particle Rest Mass, `mu`
- `end_lambda::Float64`: Affine Parameter value, where integration will end
- `step_size::Float64`: Step Size (Fixed)

# Output
- `String`: Return status code of the solver
- `Array{Float64,1}`: Array, containing affine parameter values, where integration was performed
- `Array{Float64,2}`: 2D Array, containing integrated 3-Positions and Covariant 3-Momenta

"""
function solveSystem(q, p, params, end_lambda, step_size)
    a = params[1]
    outer_event_horizon = 1 + sqrt(1 - a ^ 2)

    # Affine Parameter values, where integration will be performed
    lambdas = (0., end_lambda)

    prob = HamiltonianProblem(KerrHamiltonian, q, p, lambdas, params)

    # Condition for terminating integration, when outer event horizon is reached
    reached_event_horizon(u,t,integrator) = u[1][1] - 1.01 * outer_event_horizon
    cb = ContinuousCallback(reached_event_horizon, terminate!)

    sol = solve(prob, VerletLeapfrog(), dt = step_size, callback=cb)
    
    # sol.u is of type, ArrayPartition
    # These operations are necessary for convenient numpy usage
    solu_flat = Array(transpose(reshape(collect(Iterators.flatten(sol.u)), (6, :))))

    return String(sol.retcode), sol.t, solu_flat
end

end # module
