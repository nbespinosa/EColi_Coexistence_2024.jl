## .- .-. - -.-- - - -. .. -.- .-.- .-. - --- - .- -.- - 
@time begin
    using GLPK
    const A = COBREXA.A
    using JSONFBCModels
    using Gurobi
end

## .- .-. - -.-- - - -. .. -.- .-.- .-. - --- - .- -.- - 
# Load model
let
    fn = joinpath(@__DIR__, "iJO1366.json")
    global model = load_model(fn)
end

## .- .-. - -.-- - - -. .. -.- .-.- .-. - --- - .- -.- - 
# get stoichiometry matrix
let
    S = A.stoichiometry(model)
end

## .- .-. - -.-- - - -. .. -.- .-.- .-. - --- - .- -.- - 
# get bounds
let
    lb, ub = A.bounds(model)
end

## .- .-. - -.-- - - -. .. -.- .-.- .-. - --- - .- -.- - 
# fba
let
    solution = flux_balance_analysis(model, optimizer = GLPK.Optimizer)
    solution[:fluxes]
end

## .- .-. - -.-- - - -. .. -.- .-.- .-. - --- - .- -.- - 