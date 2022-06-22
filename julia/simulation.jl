import Pkg
Pkg.pkg"activate ./julia/StructuralEquationModels.jl"
#Pkg.add.(["SparseArrays", "Symbolics", "CSV", "DataFrames", "Optim", "LineSearches", "Plots", "BenchmarkTools", "FiniteDiff", "LinearAlgebra", "Random", "Distributions", "Statistics", "MKL"])
Pkg.instantiate()
using MKL # on its own line
using StructuralEquationModels, SparseArrays, Symbolics, CSV, DataFrames, Optim, LineSearches, Plots, BenchmarkTools, FiniteDiff, LinearAlgebra, Random, Distributions, Statistics;
import StructuralEquationModels as sem
include("./StructuralEquationModels.jl/src/frontend/parser.jl")# frontend looks very different as of 2022
data_grid = DataFrame(CSV.File("./data_grid.csv"));
#data_grid = data_grid[rand(1:nrow(data_grid),10),:]
#data_grid = filter(row -> (row.n_noise == 25) & (row.coll == .9), data_grid)
function construct_parametertable(csv)
    ParameterTable(
        unique(csv[!, :from][csv.latent]),
        unique(csv[!, :from][.!csv.latent]),
        Vector(csv.from),
        Vector(csv.parameter_type),
        Vector(csv.to),
        csv.free .!= 0,
        csv.estimate,
        [ismissing(l) ? "" : l for l in csv.label],
        csv.start,
        csv.estimate
    )
end

import StructuralEquationModels.RAMMatrices
function RAMMatrices(partable::ParameterTable, parameter_name = :x)
    A, S, F, parameters = get_RAM(partable, parameter_name)
    RAMMatrices(A = Matrix(A), S = Matrix(S), F = Matrix(F), parameters = parameters)
end

function fill_cache!(imply_cache, grid)
    for i in 1:size(grid, 1)
        fill_cache!(imply_cache, grid.par_file[i], grid.model[i])
    end
    nothing
end

function fill_cache!(imply_cache, par_path, model)
    if haskey(imply_cache, model) return nothing
    else
        println("Calc imply")
        par_file = DataFrame(CSV.File(par_path))
        partable = construct_parametertable(par_file)
        ram = RAMMatrices(partable)
        imply = RAMSymbolic(; ram_matrices = ram, start_val = start_simple)
        imply_cache[model] = (imply, imply.start_val)
        return nothing
    end
end

imply_cache = Dict()
which = 1
fill_cache!(imply_cache, data_grid[which, :par_file], data_grid[which, :model])
fill_cache!(imply_cache, data_grid)

sample_cov(true_cov, n = 100) = Statistics.cov(permutedims(rand(MvNormal(Symmetric(true_cov)), n)))

function prepare_sem(cov_path, par_path, model, n, imply_cache)
    par_file = DataFrame(CSV.File(par_path))
    cov_file = Matrix(DataFrame(CSV.File(cov_path)))
    if n != Inf
        cov_sample = sample_cov(true_cov, n)
    else
        cov_sample = cov_file
    end
    partable = construct_parametertable(par_file)
    ram = RAMMatrices(partable)
    imply = imply_cache[model]
    return par_file, cov_sample, partable, ram, imply
end

mse(x, x̂) = mean(skipmissing((x .- x̂).^2))
rmse(x, x̂) = sqrt(mse(x, x̂))

function compare_lavaan(cov_path, par_path, model, imply_cache)
    par_file, cov_sample, partable, ram, imply = prepare_sem(cov_path, par_path, model, Inf, imply_cache)
    model = Sem(
        obs_cov = cov_sample,
        imply = imply[1],
        loss = (SemML,),
        semdiff = SemDiffOptim(
            BFGS(),
            Optim.Options(
                ; f_tol = 1e-10,
                x_tol = 1.5e-8)
        )
    )
    estimate = Optim.minimizer(sem_fit(model))
    rmse(partable.estimate[partable.free], estimate)
end

compare_lavaan(grid, which) = compare_lavaan(grid[which, :cov_file], grid[which, :par_file], grid[which, :model], imply_cache)
compare_lavaan(grid) = [compare_lavaan(grid, i) for i in 1:size(grid, 1)]

compare_lavaan(data_grid, 1)
#deviation_lavaan = mean(compare_lavaan(data_grid))

converged_optim(fit) = Optim.converged(fit)
converged_withing_bounds(fit, bound = 5.0) = all(abs.(Optim.minimizer(fit)) .< bound)
converged_no_inf(fit) = abs(Optim.minimum(fit)) != Inf

function check_converged(fit)
    [converged(fit) for converged in (converged_optim, converged_withing_bounds, converged_no_inf)]
end

function statefull_fit(model::Sem{O, I, L, D}, state = nothing) where {O, I, L, D <: SemDiffOptim}
    method = model.diff.algorithm
    options = model.diff.options
    initial_x = model.imply.start_val
    d = OnceDifferentiable(
        par -> objective!(model, par),
        (grad, par) -> (grad .= StructuralEquationModels.gradient!(model, par)),
        (grad, par) -> StructuralEquationModels.objective_gradient!(grad, model, par),
        initial_x)
    if state === nothing 
        state = Optim.initial_state(method, options, d, initial_x)
    else
        invH = state.invH
        state = Optim.initial_state(method, options, d, initial_x)
        state.invH = invH
    end
    result = Optim.optimize(d, initial_x, method, options, state)
    result, state
end

function warm_start_reg(train, holdout, imply, start_val, which_ridge, alpha; state = nothing)
    imply.start_val .= start_val
    model_train = Sem(
        obs_cov = train,
        imply = imply,
        loss = (SemML, SemRidge),
        which_ridge = which_ridge,
        α_ridge = alpha,
        diff = SemDiffOptim(
            BFGS(),
            Optim.Options(
                ; f_tol = 1e-10,
                x_tol = 1.5e-8)
        )
    )
    fitted, state = statefull_fit(model_train, state)
    converged = check_converged(fitted)
    par = Optim.minimizer(fitted)
    train_loss = Optim.minimum(fitted)
    model_holdout = Sem(
        obs_cov = holdout,
        imply = imply,
        loss = (SemML,)
    )
    model_holdout(par, 1.0, nothing, nothing)
    holdout_loss = model_holdout.loss.F[1]
    DataFrame(converged = Ref(converged), train_loss = train_loss, holdout_loss = holdout_loss, par = Ref(par), to_penalize = Ref(which_ridge), alpha = alpha, state = Ref(state))
end

function warm_start_reg(train, holdout, imply, start_val, which_ridge, alphas::Vector)
    out = warm_start_reg(train, holdout, imply, start_val, which_ridge, alphas[1])
    start_val = out.par[1]
    state = out.state[1]
    #state = nothing
    select!(out, Not(:state))
    for α in alphas[2:end]
        fitted = warm_start_reg(train, holdout, imply, start_val, which_ridge, α, state = state)
        start_val = fitted.par[1]
        state = fitted.state[1]
        #state = nothing
        select!(fitted, Not(:state))
        append!(out, fitted)
    end
    out
end

function get_sem(cov_path, par_path, model, n, alphas, imply_cache)
    par_file, cov_true, partable, ram, imply = prepare_sem(cov_path, par_path, model, Inf, imply_cache)
    cov_train = sample_cov(cov_true, n)
    cov_holdout = sample_cov(cov_true, 1000)
    cov_holdout = cov_true
    # find regressions from M to c1, ..., cn
    to_penalize = findall([((from == "M") & occursin(r"^c\d+$", to)) for (from, to) in zip(par_file[par_file.free .!= 0, :].from, par_file[par_file.free .!= 0, :].to)])
    fits = warm_start_reg(cov_train, cov_holdout, imply_cache[model][1], imply_cache[model][2], to_penalize, alphas)
end

get_sem_grid(grid, which, imply_cache, alphas) = get_sem(grid[which, :cov_file], grid[which, :par_file], grid[which, :model], grid[which, :n_sample], alphas, imply_cache)
get_sem_grid!(grid, imply_cache, alphas) = grid.results = [get_sem_grid(grid, i, imply_cache, alphas) for i in 1:size(grid, 1)]
fitted = get_sem_grid(data_grid, 1, imply_cache, collect(0:0.01:.5))

data_grids = [deepcopy(data_grid) for _ in 1:200]
function fill_data_grids(data_grids)
    for (i, g) in enumerate(data_grids)
        get_sem_grid!(g, imply_cache, collect(0:0.01:.5))
        println(i)
    end
    nothing
end

fill_data_grids(data_grids)

#Arrow.write("simulation_new_interface_results.arrow", data_grids, maxdepth = 10)

#@benchmark fits = [get_sem_grid(data_grid, i, imply_cache, collect(0:0.1:5)) for i in 1:size(data_grid, 1)]

#using ProfileView
#ProfileView.@profview get_sem_grid(data_grid, 1, imply_cache, collect(0:0.01:.5))

#plot(fitted.alpha, hcat(fitted.par...)[fitted.to_penalize[1], :]')
#plot(fitted.alpha, [fitted.train_loss, fitted.holdout_loss])


function flatten_result(result, par_file)
    partable = DataFrame(CSV.File(par_file))
    true_est = partable.estimate[partable.free .!= 0]
    est_rmse = rmse.(Ref(true_est), result.par)
    closest_true = argmin(est_rmse)
    best_holdout = argmin(result.holdout_loss)
    best_train = argmin(result.train_loss)
    smallest_alpha = argmin(result.alpha)
    biggest_alpha = argmax(result.alpha)
    relevant = result[[closest_true, best_holdout, best_train, smallest_alpha, biggest_alpha], :]
    relevant.optimal = ["closest_true", "best_holdout", "best_train", "smallest_alpha", "biggest_alpha"]
    relevant.converged = mean.(relevant.converged)
    relevant.rmse_true = rmse.(Ref(true_est), relevant.par)
    select!(relevant, Not(:par))
    select!(relevant, Not(:to_penalize))
    eachrow(relevant)
end

function flatten_results(results; par_col = :par_file, result_col = :results, new_col = :summary)
    results[!, new_col] = flatten_result.(results[!, result_col], results[!, par_col])
    flattened = flatten(results, new_col)
    select!(flattened, findall(.!in.(names(flattened), Ref(String.([result_col; new_col])))))
    summary = vcat(DataFrame.(flatten(results, new_col)[!, new_col])...)
    hcat(flattened, summary)
end

out = vcat(flatten_results.(data_grids)...)
# how many have converged? (if we count any failed convergence flag as non-converged)
1 - mean(out[:, :converged] .< 1)
CSV.write("simulation_results.csv", out)

