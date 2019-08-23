using Dierckx
Using DifferentialEquations
function wt(params)
	params
end

# parameters
parameters= [
:k3=> 0.2, 
:k2=> 0.6, 
:ksnf1=> 0.4, 
:ksnf1std1=> 0.1, 
:nsnf1=> 2.3, 
:nsnf2=> 2, 
:Snf1tot=> 12, 
:dmth1=> 0.3, 
:nmth1snf3=> 1.2, 
:nmth1rgt2=> 1.4, 
:dmth1snf3=> 0.6, 
:dmth1rgt2=> 0.7, 
:smth1=> 0.3, 
:kmig1mth1=> 0.2, 
:nmig1mth1=> 1.3, 
:kmig2mth1=> 12, 
:nmig2mth1=> 1.3, 
:std1tot=> 2.1, 
:istd1=> 3.4, 
:nstd1=> 1.5, 
:estd1max=> 3.1, 
:mig1tot=> 1, 
:imig1=> 3.4, 
:kmig1snf1=> 0.8, 
:emig1max=> 3.1, 
:dmig2=> 0.3, 
:dmig2snf1=> 2.1, 
:kmig2snf1=> 2.4, 
:smig2=> 0.5, 
:kmig2std1=> std1tot/2, 
:nmig2std1=> 2.1, 
:kmig2mth1std1=> 0.3, 
:nmig2mth1std1=> 1.2, 
:dhxt4=> 0.03, 
:dhxt4max=> 2.1, 
:kdhxt4=> 2.4, 
:ndhxt4=> 1, 
:shxt4=> 0.5, 
:khxt4mth1=> 12, 
:nhxt4mth1=> 1.3, 
:khxt4std1=> std1tot/2, 
:nhxt4std1=> 2.1, 
:khxt4mth1std1=> 0.3, 
:nhxt4mth1std1=> 1.2, 
:Hxt4_0=> 0, 
:Mig1_0=> 0, 
:Mig2_0=> 0, 
:Mth1_0=> smth1/dmth1, 
:Std1_0=> std1tot, 
]



getvalue(pars, key)= [p[2] for p in pars if p[1]==key][1]

function makeproblem(modelfile, inputx,inputy, concentration, t, genotype=wt)
#making a general spline for every purpose
global input = Spline1D(inputx,inputy*concentration; k=1)

#SYSTEM OF ODES WITH A TIME VARYING INPUT
include(modelfile)

#OUTPUTING A FUNCTION THAT GENERATES A PROBLEM AS A FUNCTION OF ONLY PARAMETERS.
prob(params)= ODEProblem(model, exp.([getvalue(params, x) for x in [ :Hxt4_0, :Mig1_0, :Mig2_0, :Mth1_0, :Std1_0]]), t, exp.([j[2] for j in genotype(params)]  ))
end

