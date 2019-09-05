using Dierckx
using DifferentialEquations
function wt(params)
   	params
end

include("genotypes.jl")

#getvalue(pars, key) = [p[2] for p in pars if p[1] == key][1]


function makeproblem(modelfile, inputx, inputy, concentration, t,hxt40, genotype = wt)
#making a general spline for every purpose
    global input = Spline1D(inputx, inputy * concentration; k = 1)

#SYSTEM OF ODES WITH A TIME VARYING INPUT
    include(modelfile)

#OUTPUTING A FUNCTION THAT GENERATES A PROBLEM AS A FUNCTION OF ONLY PARAMETERS.
    prob(params) = ODEProblem(model, [hxt40;exp.(params[45:48])], t, [exp(j) for j in genotype(params)])
end