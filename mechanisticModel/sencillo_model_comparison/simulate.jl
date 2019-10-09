using JSON
using Distributions
using Sundials
using StatsPlots
#script assumes you are at the root of the Mig1Model repository
cd("/Users/s1259407/Dropbox/PhD/phd_peter_swain/data/plate_reader_data/PythonScripts/Mig1Model/")
#file thad describes all the modifications of parameters for each genotype.
include("./mechanisticModel/sencillo_model_comparison/genotypes.jl")

function trysolve(prob, pars, datamean, datastd) #main function to solve a system with two solvers and return a large value if both fail. 

	try
		sol= solve(prob(pars),ABDF2(), saveat=interv)
		arr=[[j[i] for j in sol.u] for i=1:length(sol.u[1])]
		#push!(finalarrs, arr[1])
		#plot!(arr[1], xlims=(0,20))
		lsqval=sum(((arr[1]-datamean)/datastd).^2)
	catch
			try
				sol= solve(prob(pars),TRBDF2(), saveat=interv)
				arr=[[j[i] for j in sol.u] for i=1:length(sol.u[1])]
				lsqval=sum(((arr[1]-datamean)/datastd).^2)
			catch
				lsqval=1000000.0
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


function trysolveplot(prob, pars, datamean, datastd) #main function to solve a system with two solvers and return a large value if both fail. 

	try
		sol= solve(prob(pars),ABDF2(), saveat=interv)
		arr=[[j[i] for j in sol.u] for i=1:length(sol.u[1])]
		#push!(finalarrs, arr[1])
		#plot!(arr[1], xlims=(0,20))
		lsqval=sum(((arr[1]-datamean)/datastd).^2)
	catch
			try
				sol= solve(prob(pars),TRBDF2(), saveat=interv)
				arr=[[j[i] for j in sol.u] for i=1:length(sol.u[1])]
				lsqval=sum(((arr[1]-datamean)/datastd).^2)
			catch
				lsqval=1000000.0
			end
	end
	plot(sol.t, arr[1])
	plot!(sol.t, datamean)
end

function trysolveplotnew(prob, pars, datamean, datastd, outvar=1) #main function to solve a system with two solvers and return a large value if both fail. 
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


function boldsolve(prob, pars, datamean, datastd) #main function to solve a system with two solvers and return a large value if both fail. 
		solver=CVODE_BDF
		sol= solve(prob(pars), solver(), saveat = interv,verbose=false,abstol=1e-12,reltol=1e-12)
		arr=[[j[i] for j in sol.u] for i=1:length(sol.u[1])]
		#push!(finalarrs, arr[1])
		#plot!(arr[1], xlims=(0,20))
		lsqval=sum(((arr[1]-datamean)/datastd).^2)
		#plot(sol.t, arr[1], legend=false)
		#plot!(sol.t, datamean, legend=false
end


###this function solves the problem but in a straightforward way. used to troubleshoot crashing parameters or situations.
function solvedynamic(prob, pars, datamean, datastd)
global interv
sol= solve(prob(pars),ABDF2(), saveat=interv)
arr=[[j[i] for j in sol.u] for i=1:length(sol.u[1])]
#push!(finalarrs, arr[1])
#plot!(arr[1], xlims=(0,20))
lsqval=sum(((arr[1]-datamean)/datastd).^2)
end



function solvedynamicnew(prob, pars, datamean, datastd, solver=CVODE_BDF)
global interv
sol= solve(prob(pars), solver(), saveat = interv,verbose=false,abstol=1e-12,reltol=1e-12)
arr=[[j[i] for j in sol.u] for i=1:length(sol.u[1])]
#push!(finalarrs, arr[1])
#plot!(arr[1], xlims=(0,20))
lsqval=sum(((arr[1]-datamean)/datastd).^2)
end



function solveall(params)  #trying to readapt the parameters to be input as an array

paramscaffold= copy(pars5d)

paramscaffold=[paramscaffold[j][1] => params[j] for j in 1:length(paramscaffold)]

sum(newtrysolve.(allprobs,[paramscaffold], [j[1:length(x)] for j in dm2],[k[1:length(x)] for k in me2]))

end

function solveall2(params)  #trying to readapt the parameters to be input as an array

paramscaffold= copy(pars5d)

paramscaffold=[paramscaffold[j][1] => params[j] for j in 1:length(paramscaffold)]

newtrysolve.(allprobs,[paramscaffold], [j[1:length(x)] for j in dm2],[k[1:length(x)] for k in me2])


end


function solveall3(params)  #trying to readapt the parameters to be input as an array

paramscaffold= copy(pars5d)

paramscaffold=[paramscaffold[j][1] => params[j] for j in 1:length(paramscaffold)]

boldsolve.(allprobs,[paramscaffold], [j[1:length(x)] for j in dm2],[k[1:length(x)] for k in me2])


end





##  when a parameter is set to 0 and the log is taken, we cap it to -30 infinite with this function instead of leaving it 
subss(x)=  if !isfinite(x) return -30.0 else return x end


####INPUT PROPERTIES
#the mean  can be approximated by interpolation, 
#but not the standard deviation. we need to  stick to the timepoints at which this was originally sampled.


glucose=[8.32733959421089e-19,2.561481685921e-19,1.04197119608466e-19,5.53738039368354e-20,2.30389033029131e-20,1.01730601065369e-20,9.64607628015416e-21,4.54189858719156e-21,2.62590228463374e-21,1.07030601256657e-21,7.12893122193276e-22,3.77665654101745e-22,1.13613450886179e-22,9.11107194068206e-23,1.68345898188741e-23,4.59571704388656e-23,2.4737145794224e-23,5.06371909794061e-24,5.65818382020064e-24,1.04048279473923e-24,3.94925490074619e-24,3.45909892619726e-25,7.07762515569138e-25,7.59557441484042e-26,2.02035130943482e-27,3.27268072550391e-26,1.56537448533959e-27,1.6320588929363e-27,4.39103892339638e-28,4.80386219709977e-29,3.78228642078199e-28,6.7008774844888e-33,3.38770680654291e-30,3.11390728411053e-36,3.43524261389185e-36,9.0953685838273e-37,0,2.65814991285393e-16,1.40534785861741e-05,0.956449369761021,0.999932720540147,0.999988138202835,0.999994336997356,0.999995688515055,0.99999655799382,0.999997230933311,0.999997349602989,0.999997796320276,0.999997725777247,0.999997825419588,0.999997959294029,0.99999807290326,0.999998148682981,0.999998097810728,0.999998115954266,0.999998085169392,0.999998126522822,0.99999814575262,0.999998172022714,0.999998338656721,0.999998383673738,0.999998504114664,0.999998515535666,0.999998447595621,0.999998528661984,0.99999855188616,0.999998519653104,0.999998565797203,0.999998565177499,0.999998598958849,0.999998669110336,0.999998719774422,0.999998776796943,0.999998794501031,0.99999883534008,0.999998798817409,0.99999880270023,0.999998808299886,0.999998809973147,0.999998846846304,0.999998882337051,0.999998957349787,0.999998958954537,0.999998984982437,0.999998977466028,0.999998955080815,0.999998965848435,0.999998965848435,0.999998955080815,0.999998977466028,0.999998984982437,0.999998958954537,0.999998957349787,0.999998882337051,0.999998846846304,0.999998809973147,0.999998808299886,0.99999880270023,0.999998798817409,0.99999883534008,0.999998794501031,0.999998776796943,0.999998719774422,0.999998669110336,0.999998598958849,0.999998565177499,0.999998565797203,0.999998519653104,0.99999855188616,0.999998528661984,0.999998447595621,0.999998515535666,0.999998504114664,0.999998383673738,0.999998338656721,0.999998172022714,0.99999814575262,0.999998126522822,0.999998085169392,0.999998115954266,0.999998097810728,0.999998148682981,0.99999807290326,0.999997959294029,0.999997825419588,0.999997725777247,0.999997796320276,0.999997349602989,0.999997230933311,0.99999655799382,0.999995688515055,0.999994336997356,0.999988138202835,0.999932720540147,0.956449369761021,1.40534785861741e-05,2.65814991285393e-16,0,9.0953685838273e-37,3.43524261389185e-36,3.11390728411053e-36,3.38770680654291e-30,6.7008774844888e-33,3.78228642078199e-28,4.80386219709977e-29,4.39103892339638e-28,1.6320588929363e-27,1.56537448533959e-27,3.27268072550391e-26,2.02035130943482e-27,7.59557441484042e-26,7.07762515569138e-25,3.45909892619726e-25,3.94925490074619e-24,1.04048279473923e-24,5.65818382020064e-24,5.06371909794061e-24,2.4737145794224e-23,4.59571704388656e-23,1.68345898188741e-23,9.11107194068206e-23,1.13613450886179e-22,3.77665654101745e-22,7.12893122193276e-22,1.07030601256657e-21,2.62590228463374e-21,4.54189858719156e-21,9.64607628015416e-21,1.01730601065369e-20,2.30389033029131e-20,5.53738039368354e-20,1.04197119608466e-19,2.561481685921e-19,8.32733959421089e-19,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
#making sure the shape of the above is right
y=glucose[1:229]
#times at which original glucose is sampled
tim=[0 0.0833333333333333 0.166666666666667 0.25 0.333333333333333 0.416666666666667 0.5 0.583333333333333 0.666666666666667 0.75 0.833333333333333 0.916666666666667 1 1.08333333333333 1.16666666666667 1.25 1.33333333333333 1.41666666666667 1.5 1.58333333333333 1.66666666666667 1.75 1.83333333333333 1.91666666666667 2 2.08333333333333 2.16666666666667 2.25 2.33333333333333 2.41666666666667 2.5 2.58333333333333 2.66666666666667 2.75 2.83333333333333 2.91666666666667 3 3.08333333333333 3.16666666666667 3.25 3.33333333333333 3.41666666666667 3.5 3.58333333333333 3.66666666666667 3.75 3.83333333333333 3.91666666666667 4 4.08333333333333 4.16666666666667 4.25 4.33333333333333 4.41666666666667 4.5 4.58333333333333 4.66666666666667 4.75 4.83333333333333 4.91666666666667 5 5.08333333333333 5.16666666666667 5.25 5.33333333333333 5.41666666666667 5.5 5.58333333333333 5.66666666666667 5.75 5.83333333333333 5.91666666666667 6 6.08333333333333 6.16666666666667 6.25 6.33333333333333 6.41666666666667 6.5 6.58333333333333 6.66666666666667 6.75 6.83333333333333 6.91666666666667 7 7.08333333333333 7.16666666666667 7.25 7.33333333333333 7.41666666666667 7.5 7.58333333333333 7.66666666666667 7.75 7.83333333333333 7.91666666666667 8 8.08333333333333 8.16666666666667 8.25 8.33333333333333 8.41666666666667 8.5 8.58333333333333 8.66666666666667 8.75 8.83333333333333 8.916666666666679 9.08333333333333 9.16666666666667 9.25 9.33333333333333 9.41666666666667 9.5 9.58333333333333 9.66666666666667 9.75 9.83333333333333 9.91666666666667 10 10.0833333333333 10.1666666666667 10.25 10.3333333333333 10.4166666666667 10.5 10.5833333333333 10.6666666666667 10.75 10.8333333333333 10.9166666666667 11 11.0833333333333 11.1666666666667 11.25 11.3333333333333 11.4166666666667 11.5 11.5833333333333 11.6666666666667 11.75 11.8333333333333 11.9166666666667 12 12.0833333333333 12.1666666666667 12.25 12.3333333333333 12.4166666666667 12.5 12.5833333333333 12.6666666666667 12.75 12.8333333333333 12.9166666666667 13 13.0833333333333 13.1666666666667 13.25 13.3333333333333 13.4166666666667 13.5 13.5833333333333 13.6666666666667 13.75 13.8333333333333 13.9166666666667 14 14.0833333333333 14.1666666666667 14.25 14.3333333333333 14.4166666666667 14.5 14.5833333333333 14.6666666666667 14.75 14.8333333333333 14.9166666666667 15 15.0833333333333 15.1666666666667 15.25 15.3333333333333 15.4166666666667 15.5 15.5833333333333 15.6666666666667 15.75 15.8333333333333 15.9166666666667 16 16.0833333333333 16.1666666666667 16.25 16.3333333333333 16.4166666666667 16.5 16.5833333333333 16.6666666666667 16.75 16.8333333333333 16.9166666666667 17 17.0833333333333 17.1666666666667 17.25 17.3333333333333 17.4166666666667 17.5 17.5833333333333 17.6666666666667 17.75 17.8333333333333 17.9166666666667 18 18.0833333333333 18.1666666666667 18.25 18.3333333333333 18.4166666666667 18.5 18.5833333333333 18.6666666666667 18.75 18.8333333333333 18.9166666666667 19 19.0833333333333]
#maing sure the shape of the above is right.
tim=tim[1:229]


solver=ABDF2


#importing data through json
#these are cell traces
celldata=JSON.parsefile("./json/hxtmeandata.json", dicttype=Dict)
#=these are the means for the 18 conditions.
 #each 3 is one genotype, and for each genotype
 # there are three glucose concentrations (0.2%, 0.4%, 1%)
=#
datameans=JSON.parsefile("./json/allfitmeans.json", dicttype=Dict)
#ordering the data so that each entry is one phenotype
#dm2=[[datameans[j][k] for j in 1:size(datameans)[1]] for k in 1:size(datameans[1])[1]]
dm2=datameans


subsd(x)=  if typeof(x)==String return 1000 else return x end # if it is a string(missing value) return a st of 1000  to minimise the relevance of this point.
subzero(x)=  if x==0 return 1 else return x end # if it is 0, return 1  in order not to divide by zero.

meanerr=JSON.parsefile("./json/meanerr.json", dicttype=Dict)#doesn't need ordering 

me2=[subzero.(subsd.(j)) for j in meanerr]


stdmeans=JSON.parsefile("./json/stdmeans.json", dicttype=Dict)#doesn't need ordering 

sm2=[subzero.(subsd.(j)) for j in stdmeans]


global interv= 0.044543429844097995; # 0.0445... for 450 datapoints   #0.05730659025787966   #0.0573... for 350 datapoints  # 0.0803212851405622  # 0.0803 #spacing to generate the right times comes from linspace(0, 20, 250)

nsteps=convert(Int64, ceil(18.39/interv))
tims= repeat([interv], nsteps)
x=cumsum(tims) #general time vector for matching simulation and data
#x is the exact timepoints at which the means and standard deviations are originally sampled.
#we use x for all the simulations henceforth.


#Data means for 18 conditions
#=for each data mean, make a
# spline interpolant that can be sampled at any arbitrary time.
#and can be called in a loop.
=#
subs(x)=  if typeof(x)==String return missing else return x end

global splines=Array{Spline1D}(undef, (18)) #making splines for original data
for j in 1:18
if j==7
points=findall(.!ismissing.([subs(x) for x in dm2[j][1:nsteps]]))#not nan points
splines[j]=Spline1D(x[ points], dm2[j][ points])
else
splines[j]=Spline1D(x[1:nsteps], dm2[j][1:nsteps])
end
end


##plotting the data
#interps=[splines[j](x) for j in 1:18]
#plot(x, interps)

##Creating 18 different problems: 6 genotypes x 3 concentrations

tt=(0.0, 18.39) #this time is roughly the last timepoint of x
genotypes= repeat([wt,mth1ko, mig1ko,  std1ko, rgt2ko, snf3ko], inner=3) #genotype functions in the order given by data
concs=repeat([0.2, 0.4, 1], 6) # glucose concentrations  in the order given by data.

#add paths to more files for example.
modelfile="./mechanisticModel/sencillo_model_comparison/hxt4model3d.jl"

#making array of problems ready to receive parameter sets
allprobs=makeproblem.([modelfile], [tim[1:229]], [y], concs, [tt], genotypes)

#taking the log of the parameters and capping the zeros to -30
pars=[p[1]=>subss(log(p[2])) for p in parameters]

##PARAMETER BOUNDS
# values in the data are normalised to the values of mig1  before glucose, which are pretty constant. estimated as 1= 100 molecules

# [1.0e-3, 120] #synthesis of ribosomes
# [0.1, 2.5e3] #degradation
# [0.1, 2.5e4] #nuclear import/export
# [1.0e-3, 1.0e3] #glucose K
# [1.0e-3, 1.0e3] # transcription factor affities in Mig1 units
# [1.0, 10.0] #hill factors

 bounds=[
:k3=>(1.0e-3, 1.0e3), :k2=> (1.0e-3, 1.0e3), :ksnf1=>(1.0e-3, 1.0e3), :ksnf1std1=> (1.0e-3, 1.0e3), :nsnf1=>(1.0, 10.0), :nsnf2=>(1.0, 10.0), :Snf1tot=> (0.1, 1.0e2), #10 times less than mig1 or 100 times more.
:dmth1=>  (0.1, 2.5e3), :nmth1snf3=> (1.0, 10.0), :nmth1rgt2=> (1.0, 10.0), :dmth1snf3=> (0.1, 2.5e3), :dmth1rgt2=>  (0.1, 2.5e3), 
:smth1=> (1.0e-3, 120), :kmig1mth1=>  (1.0e-3, 1.0e3), :nmig1mth1=> (1.0, 10.0), :kmig2mth1=>  (1.0e-3, 1.0e3), :nmig2mth1=> (1.0, 10.0), :std1tot=> (0.1, 1.0e3), :istd1=> (0.1, 2.5e4), :nstd1=> (1.0, 10.0), :estd1max=>(0.1, 2.5e4), :mig1tot=> (.9999999999,1), 
:imig1=> (0.1, 2.5e4) , :kmig1snf1=>  (1.0e-3, 1.0e3), :emig1max=>(0.1, 2.5e4), :dmig2=>  (0.1, 2.5e3), :dmig2snf1=>  (0.1, 2.5e3), :kmig2snf1=>   (1.0e-3, 1.0e3), :smig2=>(1.0e-3, 120), :kmig2std1=>  (1.0e-3, 1.0e3), :nmig2std1=>(1.0, 10.0), :kmig2mth1std1=>  (1.0e-3, 1.0e3), :nmig2mth1std1=> (1.0, 10.0), :dhxt4=>  (0.1, 2.5e3), :dhxt4max=> (0.1, 2.5e3), :kdhxt4=>  (1.0e-3, 1.0e3), :ndhxt4=> (1.0, 10.0), :shxt4=> (1.0e-3, 120), 
:khxt4mth1=> (1.0e-3, 1.0e3), :nhxt4mth1=> (1.0, 10.0), :khxt4std1=>  (1.0e-3, 1.0e3), 
:nhxt4std1=>(1.0, 10.0), :khxt4mth1std1=>  (1.0e-3, 1.0e3), :nhxt4mth1std1=> (1.0, 10.0), #Not exactly sure about what bounds to use for initial conditions.  doing form 0 to 10 in mig1 units.
:Hxt4_0=> (0.0, 1.0e1), 
:Mig1_0=> (0.0, 1.0e1),  
:Mig2_0=> (0.0, 1.0e1), 
:Mth1_0=> (0.0, 1.0e1),  
:Std1_0=> (0.0, 1.0e1), 
]
#from the above universal collection of bounds, make a local collection of bounds to match the parameters in the current model

localbounds=[p[1]=> log.(getvalue(bounds, p[1]))  for p in parameters];

#localbounds=[p[1]=> (p[2]-10,p[2]+10) for p in rg2]

#making (prior) distributions to sample parameters from, for which the parameter value is at the centre.
#must replaced by something more sophisticated here!
pardistributions=[p[1] => Normal(p[2],1) for p in pars];

#making priori Uniform distributions to sample parameters from, based on the above bounds.
#pardistributions=[p[1] => Uniform(p[2][1],p[2][2]) for p in localbounds]

#taking a sample of the parameter distributions. 
#vals=[]
#for j in 1:100
parsample=[d[1]=> rand(d[2])  for d in pardistributions];
#broadcasting the parameters to all predefined systems and conditions, solving and evaluating the cost. 


newtrysolve.(allprobs,[parsample], [j[1:length(x)] for j in dm2],[k[1:length(x)] for k in me2])

trysolve.(allprobs,[parsample], dm2, me2)


#trying a hundred 

vals=[]
for j in 1:1000
parsample=[d[1]=> rand(d[2])  for d in pardistributions]
#trying parsample sometimes sometimes gives reasonable values. trying pars provides a solvable parameter set.

push!(vals, sum(trysolve.(allprobs,[parsample], [j[1:length(x)] for j in dm2],[k[1:length(x)] for k in me2]) ))

end






###bringing parameters



parnames= 
[
:dhxt4, :degmth1, :dmig2, :k3, :k2, :ksnf1, :ksnf1std1, :nsnf1, :nmth1snf3, :nmth1rgt2, :smth1, :kmig1mth1, :nmig1mth1, :kmig2mth1, :nmig2mth1, :std1tot, :istd1, :nstd1, :estd1max, :mig1tot, :imig1, 
 :kmig1snf1, :emig1max, :dmig2snf1, :kmig2snf1, :smig2, :kmig2std1, :nmig2std1, :kmig2mth1std1, :nmig2mth1std1, 
 :dhxt4max, :kdhxt4, :ndhxt4, :shxt4, :khxt4mth1, :nhxt4mth1, :khxt4std1, :nhxt4std1, :khxt4mth1std1, :nhxt4mth1std1, :khxt4mig1, :khxt4mig2, :nhxt4mig1, :nhxt4mig2]

guesspars=[-1.1546942232439559, -1.0887322970066144, 0.7132935756630442, -8.383416706182768, -8.808551521233642, -7.663202682401169, -0.8490911301724443, 0.01942088874544573, 0.8095384545117216, 1.1452570081581448, -2.542702044988643, -9.170542009379027, 0.9808344028633106, 0.28656559774427465, 0.37786454104517064, 0.10068945545460253, -0.9264555562961837, 0.8472390538177788, -0.35235592543987754, 0.22184004676371108, 3.119263346803683, -3.958848702341918, 2.800400741543633, 1.957424585540502, 1.5110268619211527, -0.10241264506686007, 0.5794946667353741, 1.1430127832452404, -0.7549280339651918, 1.7370792587195765, -2.9277235184473414, -10.063799157410548, 0.2902557443376703, 2.137503781773038, -3.045651546769331, 1.0621630786570095, -0.07740907614272813, 1.3161948485904997, -2.369571961625992, 0.5449457722738957, -0.46257286319590596, -1.108626011434643, 2.0978946817087203, 2.5243779161783877]


guessarray=[parnames[j] => guesspars[j] for j in 1:size(parnames)]


function getvalue2(key)

 if key in parnames
 return OrderedDict(guessarray)[key]
  else return 0
end


reorderedguess=[j[1]=> getvalue2(j[1]) for j in pars]
#reorderedguess
reorderedguess=[
:k3=>-8.38342, :k2=>-8.80855, :ksnf1=>-7.6632, :ksnf1std1=>-0.849091, :nsnf1=>0.0194209, :nsnf2=>0.0, :Snf1tot=>1.0, :dmth1=>0.0, :nmth1snf3=>0.809538, :nmth1rgt2=>1.14526, :dmth1snf3=>0.0, :dmth1rgt2=>0.0, :smth1=>-2.5427, :kmig1mth1=>-9.17054, :nmig1mth1=>0.980834, :kmig2mth1=>0.286566, :nmig2mth1=>0.377865, :std1tot=>0.100689, :istd1=>-0.926456, :nstd1=>0.847239, :estd1max=>-0.352356, :mig1tot=>0.22184, :imig1=>3.11926, :kmig1snf1=>-3.95885, :emig1max=>2.8004, :dmig2=>0.713294, :dmig2snf1=>1.95742, :kmig2snf1=>1.51103, :smig2=>-0.102413, :kmig2std1=>0.579495, :nmig2std1=>1.14301, :kmig2mth1std1=>-0.754928, :nmig2mth1std1=>1.73708, :dhxt4=>-1.15469, :dhxt4max=>-2.92772, :kdhxt4=>-10.0638, :ndhxt4=>0.290256, :shxt4=>2.1375, :khxt4mth1=>-3.04565, :nhxt4mth1=>1.06216, :khxt4std1=>-0.0774091, :nhxt4std1=>1.31619, :khxt4mth1std1=>-2.36957, :nhxt4mth1std1=>0.544946, :Hxt4_0=>-10.0, :Mig1_0=>-10.0, :Mig2_0=>-10.0, :Mth1_0=>2.0, :Std1_0=>2.0
]


rg2=parrep(reorderedguess, pars, [3, 14, 25, 27])




function tweakpar(parset)
od=OrderedDict(parset)
od[sample(od)[1]]=0.0
collect(od)
end


function combinepars(par1, par2)
#sample a parameter set from 2 source parameter sets at random.
final=[]
wch=[]
for j in 1:length(par1)

push!(final, par1[j][1]=> sample([par1[j][2], par2[j][2]]) )
push!(wch, final[j][2]==par1[j][2])

end
return final,wch
end

mixes=[]
for j=1:100
push!(mixes, combinepars(reorderedguess, pars))
end

mixarray=[trysolve.(allprobs,[n[1]], [j[1:length(x)] for j in dm2],[k[1:length(x)] for k in me2]) for n in mixes]


unsolved=findall([sum(j)> 1.08e7 for j in mixarray])


function sumup(vec)
res= zeros(49)
for j=1:length(vec)
res[vec[j]]+=1
end
res
end


##single parameter substitution
mixes=[]
for ii in 1:49

parset=OrderedDict(reorderedguess)
p=collect(parset)

p[ii]= par2[ii][1]=> par2[ii][2]


push!(mixes, trysolve.(allprobs,[p], [j[1:length(x)] for j in dm2],[k[1:length(x)] for k in me2]))

end

##two parameter substitution at random
mixes7=[]
sampls7=[]
sps=7
for ii in 1:10000

parset=OrderedDict(reorderedguess)
p=collect(parset)
spl=sample(1:49, sps)
for kk in spl
p[kk]= par2[kk][1]=> par2[kk][2]
end
push!(sampls7, spl)
push!(mixes7, trysolve.(allprobs,[p], [j[1:length(x)] for j in dm2],[k[1:length(x)] for k in me2]))

end
sum([sum(j) < 1.08e7  for j in mixes])



function parrep(par1, par2, spl)

parset=OrderedDict(par1)
p=collect(parset)
for kk in spl
p[kk]= par1[kk][1]=> par2[kk][2]
end
p
end





##bfpairs after fitting with black box

bfpairs[:k3=>0.684139, :k2=>0.10231, :ksnf1=>-2.56811, :ksnf1std1=>4.92665, :nsnf1=>2.20512, :nsnf2=>0.363983, :Snf1tot=>0.822695, :dmth1=>4.70502, :nmth1snf3=>1.21847, :nmth1rgt2=>1.05906, :dmth1snf3=>4.55476, :dmth1rgt2=>-2.10952, :smth1=>0.691909, :kmig1mth1=>5.8695, :nmig1mth1=>0.737479, :kmig2mth1=>0.763351, :nmig2mth1=>1.59547, :std1tot=>-1.36626, :istd1=>-1.21318, :nstd1=>1.46859, :estd1max=>-1.45267, :mig1tot=>-4.12998e-11, :imig1=>0.721844, :kmig1snf1=>-2.09185, :emig1max=>7.90656, :dmig2=>-1.25986, :dmig2snf1=>1.87601, :kmig2snf1=>-4.28422, :smig2=>2.17417, :kmig2std1=>0.160212, :nmig2std1=>2.01052, :kmig2mth1std1=>6.51012, :nmig2mth1std1=>0.198612, :dhxt4=>-2.26107, :dhxt4max=>-2.24995, :kdhxt4=>-6.28087, :ndhxt4=>1.86336, :shxt4=>-1.43697, :khxt4mth1=>5.56162, :nhxt4mth1=>1.03514, :khxt4std1=>4.62843, :nhxt4std1=>1.10152, :khxt4mth1std1=>-6.47512, :nhxt4mth1std1=>2.06227, :Hxt4_0=>-27.0498, :Mig1_0=>-13.8129, :Mig2_0=>-8.29938, :Mth1_0=>-2.098, :Std1_0=>-7.01498]






rg2= #reorderedguess
[:k3=>-8.38342, :k2=>-8.80855, :ksnf1=>-0.916291, :ksnf1std1=>-0.849091, :nsnf1=>0.0194209, :nsnf2=>0.0, :Snf1tot=>1.0, :dmth1=>0.0, :nmth1snf3=>0.809538, :nmth1rgt2=>1.14526, :dmth1snf3=>0.0, :dmth1rgt2=>0.0, :smth1=>-2.5427, :kmig1mth1=>-1.60944, :nmig1mth1=>0.980834, :kmig2mth1=>0.286566, :nmig2mth1=>0.377865, :std1tot=>0.100689, :istd1=>-0.926456, :nstd1=>0.847239, :estd1max=>-0.352356, :mig1tot=>0.22184, :imig1=>3.11926, :kmig1snf1=>-3.95885, :emig1max=>1.1314, :dmig2=>0.713294, :dmig2snf1=>0.741937, :kmig2snf1=>1.51103, :smig2=>-0.102413, :kmig2std1=>0.579495, :nmig2std1=>1.14301, :kmig2mth1std1=>-0.754928, :nmig2mth1std1=>1.73708, :dhxt4=>-1.15469, :dhxt4max=>-2.92772, :kdhxt4=>-10.0638, :ndhxt4=>0.290256, :shxt4=>2.1375, :khxt4mth1=>-3.04565, :nhxt4mth1=>1.06216, :khxt4std1=>-0.0774091, :nhxt4std1=>1.31619, :khxt4mth1std1=>-2.36957, :nhxt4mth1std1=>0.544946, :Hxt4_0=>-10.0, :Mig1_0=>-10.0, :Mig2_0=>-10.0, :Mth1_0=>2.0, :Std1_0=>2.0]




#my manual parameters
pars5d= [
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
:imig1=> 3.4, 
:kmig1snf1=> 0.8, 
:emig1max=> 3.1, 
:dmig2=> 0.3, 
:dmig2snf1=> 2.1, 
:kmig2snf1=> 2.4, 
:smig2=> 0.5, 
:kmig2std1=> 2.1/2, 
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
:khxt4std1=> 2.1/2, 
:nhxt4std1=> 2.1, 
:khxt4mth1std1=> 0.3, 
:nhxt4mth1std1=> 1.2, 
:khxt4mig1=> 2, 
:khxt4mig2=> 1, 
:nhxt4mig1=> 1.5, 
:nhxt4mig2=> 2, 
:Hxt4_0=> 0, 
:Mig1_0=> 0, 
:Mig2_0=> 0, 
:Mth1_0=> 3/.3, 
:Std1_0=> 2.1, 
]

mypars= [
:k3=> 0.1, 
:k2=> 0.6, 
:ksnf1=> 0.05, 
:ksnf1std1=> 160.0, 
:nsnf1=> 8.0, 
:nsnf2=> 8.0, 
:Snf1tot=> 120, 
:dmth1=> 0.3, 
:nmth1snf3=> 1.2, 
:nmth1rgt2=> 1.4, 
:dmth1snf3=> 4, 
:dmth1rgt2=> 2, 
:smth1=> 1, 
:kmig1mth1=> 0.2, 
:nmig1mth1=> 1.3, 
:kmig2mth1=> 12, 
:nmig2mth1=> 1.3, 
:std1tot=> 2.1, 
:istd1=> 4, 
:nstd1=> 3, 
:estd1max=> 3.1, 
:mig1tot=> 1, 
:imig1=> 9, 
:kmig1snf1=> 3, 
:emig1max=> 7, 
:dmig2=> 0.3, 
:dmig2snf1=> 2.1, 
:kmig2snf1=> 2.4, 
:smig2=> 0.5, 
:kmig2std1=> 2.1/2, 
:nmig2std1=> 2.1, 
:kmig2mth1std1=> 0.3, 
:nmig2mth1std1=> 1.2, 
:dhxt4=> 0.03, 
:dhxt4max=> 2.1, 
:kdhxt4=> 2.4, 
:ndhxt4=> 1, 
:shxt4=> 20, 
:khxt4mth1=> .5, 
:nhxt4mth1=> 3, 
:khxt4std1=> 2.1/2, 
:nhxt4std1=> 2.1, 
:khxt4mth1std1=> 0.3, 
:nhxt4mth1std1=> 1.2, 
:Hxt4_0=> 0, 
:Mig1_0=> 0, 
:Mig2_0=> 0, 
:Mth1_0=> 0.3/.3, 
:Std1_0=>2.1, 
]

#function fillparswith()
completepars=[];
for j in 1:length(pars5d)
try
push!(completepars, pars5d[j][1]=> getvalue(mypars,  pars5d[j][1])) 
catch
push!(completepars, pars5d[j][1]=> getvalue(pars5d,  pars5d[j][1])) 
end
end

conc=3
logpars=[p[1]=>log(p[2]) for p in completepars]
pp2=solveplotnew(allprobs[conc], logpars, dm2[1][1:length(x)],me2[1][1:length(x)], 2)
pp3=solveplotnew(allprobs[conc], logpars, dm2[1][1:length(x)],me2[1][1:length(x)], 3)
pp4=solveplotnew(allprobs[conc], logpars, dm2[1][1:length(x)],me2[1][1:length(x)], 4)
pp5=solveplotnew(allprobs[conc], logpars, dm2[1][1:length(x)],me2[1][1:length(x)], 5)
pp6=solveplotnew(allprobs[conc], logpars, dm2[1][1:length(x)],me2[1][1:length(x)], 6)
pp1=solveplotnew(allprobs[conc], logpars, dm2[1][1:length(x)],me2[1][1:length(x)], 1)
plot(pp2,pp3, pp4,pp5, pp6,pp1, layout=(3, 2), ylim=[0,5])

##
##best fit parameters for model hxtmodel5d
bestpars6=[-1.99995, 1.06211, -5.13512, 4.67567, 5.01806, -1.84446, 5.73796, -2.87079, 0.0698773, 0.288116, 0.376694, 3.14651, -4.50517, -0.379322, 5.25098, 0.560657, -1.74271, -0.245076, 1.22922, 0.746876, 4.85039, 5.45447, -2.80427, -0.264936, 3.03049, 5.25449, 5.86184, 2.28341, 5.81895, 1.63477, 2.87004, 2.80064, -1.15641, -2.48345, -3.90747, -4.58572, 1.46517, -4.65211, 3.03817, 0.611965, 1.35178, 6.74873, -0.211267, -1.23877, -2.57697, 1.36639, 2.07541, -100.0, -100.0, -100.0, -4.54752, 5.57239]


bfpairs6=[pars5d[j][1]=> bestpars6[j] for j in 1:length(pars5d)]
