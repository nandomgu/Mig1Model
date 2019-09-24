using Dierckx
using DifferentialEquations
function wt(params)
   	params
end

include("/home/msturrock/Desktop/Mig1/Mig1Model/mechanisticModel/sencillo_model_comparison/genotypes.jl")


function makeproblem(modelfile, inputx,inputy, concentration, t,hxt40, genotype=wt)
    #making a general spline for every purpose
    global input = Spline1D(inputx,inputy*concentration; k=1)
    
    #SYSTEM OF ODES WITH A TIME VARYING INPUT
    include(modelfile)
    model = model_g(input)
    #OUTPUTING A FUNCTION THAT GENERATES A PROBLEM AS A FUNCTION OF ONLY PARAMETERS.
    prob(params)= ODEProblem(model, [hxt40;exp.(genotype(params)[48:51]);0.0], t, [exp(j) for j in genotype(params)])
end