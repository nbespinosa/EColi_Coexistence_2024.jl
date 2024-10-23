## .- .-. - -.-- - - -. .. -.- .-.- .-. - --- - .- -.- - 
@time begin
    using GLPK
    const A = COBREXA.A
    using JSONFBCModels
    using GLPK
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
##~~~~~~~~~For ilv~~~~~~~~~~
using COBREXA, GLPK, Tulip, JuMP

fluxes_ilv = flux_balance_analysis_dict(
    model,
    GLPK.Optimizer;
    modifications = [
        change_objective("BIOMASS_Ec_iJO1366_WT_53p95M"), 
        change_constraint("THRD_L"; lb = 0.0, ub = 0.0), 
        knockout(["ilvA", "tdcB"]),
        change_optimizer(Tulip.Optimizer),
        change_optimizer_attribute("IPM_IterationsLimit", 1000), 
        change_sense(JuMP.MAX_SENSE), 
    ],
)

biomass_flux_ilv = fluxes_ilv["BIOMASS_Ec_iJO1366_WT_53p95M"]

#~~To know boundaries~~
sm = load_model(StandardModel, "iJO1366.json")
sm.reactions["THRD_L"].ub
sm.reactions["THRD_L"].lb

##~~~~~~~~~For met~~~~~~~~~~

fluxes_met = flux_balance_analysis_dict(
    model,
    GLPK.Optimizer;
    modifications = [
        change_objective("BIOMASS_Ec_iJO1366_WT_53p95M"), 
        change_constraint("HSST"; lb = 0.0, ub = 0.0), 
        knockout(["metA"]),
        change_optimizer(Tulip.Optimizer),
        change_optimizer_attribute("IPM_IterationsLimit", 1000), 
        change_sense(JuMP.MAX_SENSE), 
    ],
)

biomass_flux_met = fluxes_met["BIOMASS_Ec_iJO1366_WT_53p95M"]
