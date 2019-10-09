using JSON
using Distributions
using Sundials
using StatsPlots
using DataStructures
using Dierckx
using BlackBoxOptim
using DifferentialEquations

#script assumes you are at the root of the Mig1Model repository
cd("/Users/s1259407/Dropbox/PhD/phd_peter_swain/data/plate_reader_data/PythonScripts/Mig1Model/")
#file thad describes all the modifications of parameters for each genotype.
include("./mechanisticModel/sencillo_model_comparison/genotypes.jl")

subss(x)=  if !isfinite(x) return -100.0 else return x end
subnan(x)=  if typeof(x)==String return missing else return x end
subn(x)=  if x<0 return 0 else return x end
#blackbox optimisation code

getvalue(pars, key)= [p[2] for p in pars if p[1]==key][1]



function solveplotnew(prob, pars, datamean, datastd, outvar=1) #main function to solve a system with two solvers and return a large value if both fail. 
                       solver=CVODE_BDF
                       sol= solve(prob(pars), solver(), saveat = interv,verbose=false,abstol=1e-12,reltol=1e-12)
                       arr=[[j[i] for j in sol.u] for i=1:length(sol.u[1])]
                       #push!(finalarrs, arr[1])
                       #plot!(arr[1], xlims=(0,20))
                       lsqval=sum(((arr[outvar]-datamean)/datastd).^2)
                       p=plot(sol.t, arr[outvar], legend=false)
                       if outvar==1
                               p=plot!(sol.t, datamean, legend=false)
                       end
                       p
       end



function trysolveplot(prob, pars, datamean, datastd, outvar=1) #main function to solve a system with two solvers and return a large value if both fail. 
solver=CVODE_BDF
	try
		sol= solve(prob(pars), solver(), saveat = interv,verbose=false,abstol=1e-12,reltol=1e-12)
		arr=[[j[i] for j in sol.u] for i=1:length(sol.u[1])]
		#push!(finalarrs, arr[1])
		#plot!(arr[1], xlims=(0,20))
		lsqval=sum(((arr[outvar]-datamean)/datastd).^2)
			p=plot(sol.t, arr[outvar])
			p=plot!(sol.t, datamean)
	catch
			try
				sol= solve(prob(pars), solver(), saveat = interv,verbose=false,abstol=1e-12,reltol=1e-12)
				arr=[[j[i] for j in sol.u] for i=1:length(sol.u[1])]
				lsqval=sum(((arr[outvar]-datamean)/datastd).^2)
				p=plot(sol.t, arr[outvar])
				p=plot!(sol.t, datamean)
			catch
				
				lsqval=1000000.0
				p=plot(1:.1:1, 1:.1:1, color="red")
			end
	end
end


function newtrysolve(prob, pars, datamean, datastd) #main function to solve a system with two solvers and return a large value if both fail. 

	try
		 sol=solve(prob(pars), CVODE_BDF(), saveat = interv,verbose=false,abstol=1e-12,reltol=1e-12)
		arr=[[j[i] for j in sol.u] for i=1:length(sol.u[1])]
		#push!(finalarrs, arr[1])
		#plot!(arr[1], xlims=(0,20))
		lsqval=sum(((arr[1]-datamean)/datastd).^2)
	catch
			try
				 sol=solve(prob(pars), CVODE_BDF(), saveat = interv,verbose=false,abstol=1e-12,reltol=1e-12)
				arr=[[j[i] for j in sol.u] for i=1:length(sol.u[1])]
				lsqval=sum(((arr[1]-datamean)/datastd).^2)
			catch
				lsqval=1000000.0
			end
	end
end


###SIMULATION VARIABLES
glucose=[8.32733959421089e-19,2.561481685921e-19,1.04197119608466e-19,5.53738039368354e-20,2.30389033029131e-20,1.01730601065369e-20,9.64607628015416e-21,4.54189858719156e-21,2.62590228463374e-21,1.07030601256657e-21,7.12893122193276e-22,3.77665654101745e-22,1.13613450886179e-22,9.11107194068206e-23,1.68345898188741e-23,4.59571704388656e-23,2.4737145794224e-23,5.06371909794061e-24,5.65818382020064e-24,1.04048279473923e-24,3.94925490074619e-24,3.45909892619726e-25,7.07762515569138e-25,7.59557441484042e-26,2.02035130943482e-27,3.27268072550391e-26,1.56537448533959e-27,1.6320588929363e-27,4.39103892339638e-28,4.80386219709977e-29,3.78228642078199e-28,6.7008774844888e-33,3.38770680654291e-30,3.11390728411053e-36,3.43524261389185e-36,9.0953685838273e-37,0,2.65814991285393e-16,1.40534785861741e-05,0.956449369761021,0.999932720540147,0.999988138202835,0.999994336997356,0.999995688515055,0.99999655799382,0.999997230933311,0.999997349602989,0.999997796320276,0.999997725777247,0.999997825419588,0.999997959294029,0.99999807290326,0.999998148682981,0.999998097810728,0.999998115954266,0.999998085169392,0.999998126522822,0.99999814575262,0.999998172022714,0.999998338656721,0.999998383673738,0.999998504114664,0.999998515535666,0.999998447595621,0.999998528661984,0.99999855188616,0.999998519653104,0.999998565797203,0.999998565177499,0.999998598958849,0.999998669110336,0.999998719774422,0.999998776796943,0.999998794501031,0.99999883534008,0.999998798817409,0.99999880270023,0.999998808299886,0.999998809973147,0.999998846846304,0.999998882337051,0.999998957349787,0.999998958954537,0.999998984982437,0.999998977466028,0.999998955080815,0.999998965848435,0.999998965848435,0.999998955080815,0.999998977466028,0.999998984982437,0.999998958954537,0.999998957349787,0.999998882337051,0.999998846846304,0.999998809973147,0.999998808299886,0.99999880270023,0.999998798817409,0.99999883534008,0.999998794501031,0.999998776796943,0.999998719774422,0.999998669110336,0.999998598958849,0.999998565177499,0.999998565797203,0.999998519653104,0.99999855188616,0.999998528661984,0.999998447595621,0.999998515535666,0.999998504114664,0.999998383673738,0.999998338656721,0.999998172022714,0.99999814575262,0.999998126522822,0.999998085169392,0.999998115954266,0.999998097810728,0.999998148682981,0.99999807290326,0.999997959294029,0.999997825419588,0.999997725777247,0.999997796320276,0.999997349602989,0.999997230933311,0.99999655799382,0.999995688515055,0.999994336997356,0.999988138202835,0.999932720540147,0.956449369761021,1.40534785861741e-05,2.65814991285393e-16,0,9.0953685838273e-37,3.43524261389185e-36,3.11390728411053e-36,3.38770680654291e-30,6.7008774844888e-33,3.78228642078199e-28,4.80386219709977e-29,4.39103892339638e-28,1.6320588929363e-27,1.56537448533959e-27,3.27268072550391e-26,2.02035130943482e-27,7.59557441484042e-26,7.07762515569138e-25,3.45909892619726e-25,3.94925490074619e-24,1.04048279473923e-24,5.65818382020064e-24,5.06371909794061e-24,2.4737145794224e-23,4.59571704388656e-23,1.68345898188741e-23,9.11107194068206e-23,1.13613450886179e-22,3.77665654101745e-22,7.12893122193276e-22,1.07030601256657e-21,2.62590228463374e-21,4.54189858719156e-21,9.64607628015416e-21,1.01730601065369e-20,2.30389033029131e-20,5.53738039368354e-20,1.04197119608466e-19,2.561481685921e-19,8.32733959421089e-19,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
#making sure the shape of the above is right
y=glucose[1:229]
#times at which original glucose is sampled
tim=[0 0.0833333333333333 0.166666666666667 0.25 0.333333333333333 0.416666666666667 0.5 0.583333333333333 0.666666666666667 0.75 0.833333333333333 0.916666666666667 1 1.08333333333333 1.16666666666667 1.25 1.33333333333333 1.41666666666667 1.5 1.58333333333333 1.66666666666667 1.75 1.83333333333333 1.91666666666667 2 2.08333333333333 2.16666666666667 2.25 2.33333333333333 2.41666666666667 2.5 2.58333333333333 2.66666666666667 2.75 2.83333333333333 2.91666666666667 3 3.08333333333333 3.16666666666667 3.25 3.33333333333333 3.41666666666667 3.5 3.58333333333333 3.66666666666667 3.75 3.83333333333333 3.91666666666667 4 4.08333333333333 4.16666666666667 4.25 4.33333333333333 4.41666666666667 4.5 4.58333333333333 4.66666666666667 4.75 4.83333333333333 4.91666666666667 5 5.08333333333333 5.16666666666667 5.25 5.33333333333333 5.41666666666667 5.5 5.58333333333333 5.66666666666667 5.75 5.83333333333333 5.91666666666667 6 6.08333333333333 6.16666666666667 6.25 6.33333333333333 6.41666666666667 6.5 6.58333333333333 6.66666666666667 6.75 6.83333333333333 6.91666666666667 7 7.08333333333333 7.16666666666667 7.25 7.33333333333333 7.41666666666667 7.5 7.58333333333333 7.66666666666667 7.75 7.83333333333333 7.91666666666667 8 8.08333333333333 8.16666666666667 8.25 8.33333333333333 8.41666666666667 8.5 8.58333333333333 8.66666666666667 8.75 8.83333333333333 8.916666666666679 9.08333333333333 9.16666666666667 9.25 9.33333333333333 9.41666666666667 9.5 9.58333333333333 9.66666666666667 9.75 9.83333333333333 9.91666666666667 10 10.0833333333333 10.1666666666667 10.25 10.3333333333333 10.4166666666667 10.5 10.5833333333333 10.6666666666667 10.75 10.8333333333333 10.9166666666667 11 11.0833333333333 11.1666666666667 11.25 11.3333333333333 11.4166666666667 11.5 11.5833333333333 11.6666666666667 11.75 11.8333333333333 11.9166666666667 12 12.0833333333333 12.1666666666667 12.25 12.3333333333333 12.4166666666667 12.5 12.5833333333333 12.6666666666667 12.75 12.8333333333333 12.9166666666667 13 13.0833333333333 13.1666666666667 13.25 13.3333333333333 13.4166666666667 13.5 13.5833333333333 13.6666666666667 13.75 13.8333333333333 13.9166666666667 14 14.0833333333333 14.1666666666667 14.25 14.3333333333333 14.4166666666667 14.5 14.5833333333333 14.6666666666667 14.75 14.8333333333333 14.9166666666667 15 15.0833333333333 15.1666666666667 15.25 15.3333333333333 15.4166666666667 15.5 15.5833333333333 15.6666666666667 15.75 15.8333333333333 15.9166666666667 16 16.0833333333333 16.1666666666667 16.25 16.3333333333333 16.4166666666667 16.5 16.5833333333333 16.6666666666667 16.75 16.8333333333333 16.9166666666667 17 17.0833333333333 17.1666666666667 17.25 17.3333333333333 17.4166666666667 17.5 17.5833333333333 17.6666666666667 17.75 17.8333333333333 17.9166666666667 18 18.0833333333333 18.1666666666667 18.25 18.3333333333333 18.4166666666667 18.5 18.5833333333333 18.6666666666667 18.75 18.8333333333333 18.9166666666667 19 19.0833333333333]
#maing sure the shape of the above is right.
tim=tim[1:229]


solver=ABDF2


#celldata=JSON.parsefile("./json/celldata.json", dicttype=Dict)
#importing data through json
#these are cell traces
celldata=JSON.parsefile("./json/hxtdata.json", dicttype=Dict)
#these are the means for the 18 conditions.
 #each 3 is one genotype, and for each genotype
 # there are three glucose concentrations (0.2%, 0.4%, 1%)
#
datameans=JSON.parsefile("./json/allfitmeans.json", dicttype=Dict)
#ordering the data so that each entry is one phenotype
#dm2=[[datameans[j][k] for j in 1:size(datameans)[1]] for k in 1:size(datameans[1])[1]]
dm2=datameans

subnan(x)=  if typeof(x)==String return missing else return x end # if it is a string(missing value) return a st of 1000  to minimise the relevance of this point.

subsd(x)=  if typeof(x)==String return 1000 else return x end # if it is a string(missing value) return a st of 1000  to minimise the relevance of this point.
subzero(x)=  if x==0 return 1 else return x end # if it is 0, return 1  in order not to divide by zero.

meanerr=JSON.parsefile("./json/meanerr.json", dicttype=Dict)#doesn't need ordering 

me2=[subzero.(subsd.(j)) for j in meanerr]




stdmeans=JSON.parsefile("./json/stdmeans.json", dicttype=Dict)#doesn't need ordering 

sm2=[subzero.(subsd.(j)) for j in stdmeans]

function findrep(val) 
       if val==1000
       return 5
       else 
       return val
       end
       end
       
sm22=[findrep.(j) for j in sm2]


global interv= 0.044543429844097995; # 0.0445... for 450 datapoints   #0.05730659025787966   #0.0573... for 350 datapoints  # 0.0803212851405622  # 0.0803 #spacing to generate the right times comes from linspace(0, 20, 250)

nsteps=convert(Int64, ceil(18.39/interv))
tims= repeat([interv], nsteps)
x=cumsum(tims) #general time vector for matching simulation and data

tt=(0.0, 18.39) #this time is roughly the last timepoint of x
genotypes= repeat([wt,mth1ko, mig1ko,  std1ko, rgt2ko, snf3ko], inner=3) #genotype functions in the order given by data
concs=repeat([0.2, 0.4, 1], 6) # glucose concentrations  in the order given by data.
modelfile="./mechanisticModel/sencillo_model_comparison/onlymodel5.jl"


allprobs=makeproblem.([modelfile], [tim[1:229]], [y], concs, [tt], [j[1] for j in dm2], genotypes)

#parameter template

completepars=
[:k3=>0.1,
:k2=>0.6,
:ksnf1=>0.05,
:ksnf1std1=>160.0,
:nsnf1=>8.0,
:nsnf2=>8.0,
:Snf1tot=>120.0,
:dmth1=>0.3,
:nmth1snf3=>1.2,
:nmth1rgt2=>1.4,
:dmth1snf3=>4.0,
:dmth1rgt2=>2.0,
:smth1=>1.0,
:kmig1mth1=>0.2,
:nmig1mth1=>1.3,
:kmig2mth1=>12.0,
:nmig2mth1=>1.3,
:std1tot=>2.1,
:istd1=>4.0,
:nstd1=>3.0,
:nstd3=>3.0,
:estd1max=>3.1,
:imig1=>9.0,
:kmig1snf1=>3.0,
:emig1max=>7.0,
:dmig2=>0.3,
:dmig2snf1=>2.1,
:kmig2snf1=>2.4,
:smig2=>0.5,
:kmig2std1=>1.05,
:nmig2std1=>2.1,
:kmig2mth1std1=>0.3,
:nmig2mth1std1=>1.2,
:dhxt4=>0.03,
:dhxt4max=>2.1,
:kdhxt4=>2.4,
:ndhxt4=>1.0,
:shxt4=>20.0,
:khxt4mth1=>0.5,
:nhxt4mth1=>3.0,
:khxt4std1=>1.05,
:nhxt4std1=>2.1,
:khxt4mth1std1=>0.3,
:nhxt4mth1std1=>1.2,
:khxt4mig1=>2.0,
:khxt4mig2=>1.0,
:nhxt4mig1=>1.5,
:nhxt4mig2=>2.0,
:Hxt4_0=>0.0,
:Mig1_0=>0.0,
:Mig2_0=>0.0,
:Mth1_0=>1.0,
:Std1_0=>2.1]

##solveall function as a function of param vector (requires completepars)
function solveall(params)  #trying to readapt the parameters to be input as an array

paramscaffold= copy(completepars)

paramscaffold=[paramscaffold[j][1] => params[j] for j in 1:length(paramscaffold)]

sum(newtrysolve.(allprobs,[paramscaffold], [j[1:length(x)] for j in dm2],[k[1:length(x)] for k in sm22])) #used to be me2 but sm22 performs pretty well.

end


#taking the log of parameters
logpars=[p[1]=>log(p[2]) for p in completepars]

bestpars6=[-1.99995, 1.06211, -5.13512, 4.67567, 5.01806, -1.84446, 5.73796, -2.87079, 0.0698773, 0.288116, 0.376694, 3.14651, -4.50517, -0.379322, 5.25098, 0.560657, -1.74271, -0.245076, 1.22922, 0.746876, 4.85039, 5.45447, -2.80427, -0.264936, 3.03049, 5.25449, 5.86184, 2.28341, 5.81895, 1.63477, 2.87004, 2.80064, -1.15641, -2.48345, -3.90747, -4.58572, 1.46517, -4.65211, 3.03817, 0.611965, 1.35178, 6.74873, -0.211267, -1.23877, -2.57697, 1.36639, 2.07541, -100.0, -100.0, -100.0, -4.54752, 5.57239]

bfpairs6=[completepars[j][1]=> bestpars6[j] for j in 1:length(completepars)]

pp=trysolveplot.(allprobs,[bfpairs6], [j[1:length(x)] for j in dm2],[k[1:length(x)] for k in me2]); plot(pp[1:18]..., layout=(6, 3), ylim=[0, 10], legend=false)



#plotting lower and upper bounds of fitting, and where the best parameters fell

bar([p[1] for p in custombounds], legend=:topright, labels=["lower", "upper", "best"], ylim=(-12, 12))

bar!([p[2] for p in custombounds], legend=:topright, labels=["lower", "upper", "best"], ylim=(-12, 12))

bar!([[p[2] for p in  bfpairs6]], legend=:topright, labels=["lower", "upper", "best"], ylim=(-12, 12))

savefig("temp.png")



#opts7=bbsetup(solveall; Method = :adaptive_de_rand_1_bin_radiuslimited, SearchRange = [subss.((p[2]-5,p[2]+5)) for p in  logpars], NumDimensions = 52, MaxSteps = 10000)
#res7=bboptimize(opts6, MaxSteps=1000)


custombounds=
[(-7.30259,2.69741),
(-5.51083,4.48917),
(-7.99573,2.00427),
(0.0751738,10.0752),
(-2.92056,7.07944),
(-2.92056,7.07944),
(-0.212508,9.78749),
(-6.20397,3.79603),
(-4.81768,5.18232),
(-4.66353,5.33647),
(-3.61371,6.38629),
(-4.30685,5.69315),
(-5.0,5.0),
(-6.60944,3.39056),
(-4.73764,5.26236),
(-2.51509,7.48491),
(-4.73764,5.26236),
(-4.25806,5.74194),
(-3.61371,6.38629),
(-3.90139,6.09861),
(-3.8686,6.1314),
(-2.80278,7.19722),
(-3.90139,6.09861),
(-3.05409,6.94591),
(-6.20397,3.79603),
(-4.25806,5.74194),
(-4.12453,5.87547),
(-5.69315,4.30685),
(-4.95121,5.04879),
(-4.25806,5.74194),
(-6.20397,3.79603),
(-4.81768,5.18232),
(-8.50656,1.49344),
(-4.25806,5.74194),
(-4.12453,5.87547),
(-5.0,5.0),
(-2.00427,7.99573),
(-5.69315,4.30685),
(-3.90139,6.09861),
(-4.95121,5.04879),
(-4.25806,5.74194),
(-6.20397,3.79603),
(-4.81768,5.18232),
(-4.30685,5.69315),
(-5.0,5.0),
(-4.59453,5.40547),
(-4.30685,5.69315),
(-100.0,-100.0),
(-100.0,-100.0),
(-100.0,-100.0),
(-5.0,5.0),
(-4.25806,5.74194)]


