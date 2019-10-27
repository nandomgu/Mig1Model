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


function trysolvedata(prob, pars, datamean, datastd, outvar=1) #main function to solve a system with two solvers and return a large value if both fail. 

	try
		 sol=solve(prob(pars), CVODE_BDF(), saveat = interv,verbose=false,abstol=1e-12,reltol=1e-12)
		arr=[[j[i] for j in sol.u] for i=1:length(sol.u[outvar])]
		#push!(finalarrs, arr[1])
		#plot!(arr[1], xlims=(0,20))
		#lsqval=sum(((arr[1]-datamean)/datastd).^2)
		arr[outvar]
	catch
			try
				 sol=solve(prob(pars), CVODE_BDF(), saveat = interv,verbose=false,abstol=1e-12,reltol=1e-12)
				arr=[[j[i] for j in sol.u] for i=1:length(sol.u[outvar])]
				arr[outvar]
			catch
				arr=[]
				
			end
	arr[outvar]
	end
	
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
				if outvar==1
				p=plot!(sol.t, datamean)
				end
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


function solveone(prob, mn, sdev)  #trying to readapt the parameters to be input as an array

function solvv(params)

paramscaffold= copy(pars5d)

paramscaffold=[paramscaffold[j][1] => params[j] for j in 1:length(paramscaffold)]

sum(newtrysolve(prob,paramscaffold, mn[1:length(x)],sdev[1:length(x)]))

end

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
subss(x)=  if !isfinite(x) return -100.0 else return x end


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

#hxt celldata
celldata=JSON.parsefile("./json/celldata.json", dicttype=Dict)
#importing data through json
#these are cell traces
celldata=JSON.parsefile("./json/hxtdata.json", dicttype=Dict)
#=these are the means for the 18 conditions.
 #each 3 is one genotype, and for each genotype
 # there are three glucose concentrations (0.2%, 0.4%, 1%)
=#
datameans=JSON.parsefile("./json/allfitmeans.json", dicttype=Dict)
#ordering the data so that each entry is one phenotype
#dm2=[[datameans[j][k] for j in 1:size(datameans)[1]] for k in 1:size(datameans[1])[1]]
dm2=datameans

subnan(x)=  if typeof(x)==String return NaN else return x end # if it is a string(missing value) return a st of 1000  to minimise the relevance of this point.

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
:k3=>(1.0e-3, 1.0e3), 
:k2=> (1.0e-3, 1.0e3),
:ksnf1=>(1.0e-3, 1.0e3),
:ksnf1std1=> (1.0e-3, 1.0e3),
:nsnf1=>(1.0, 10.0), 
:nsnf2=>(1.0, 10.0), 
:Snf1tot=> (0.1, 1.0e2), #10 times less than mig1 or 100 times more.
:dmth1=>  (0.1, 2.5e3), 
:nmth1snf3=> (1.0, 10.0), 
:nmth1rgt2=> (1.0, 10.0), 
:dmth1snf3=> (0.1, 2.5e3), 
:dmth1rgt2=>  (0.1, 2.5e3), 
:smth1=> (1.0e-3, 120), 
:kmig1mth1=>  (1.0e-3, 1.0e3), 
:nmig1mth1=> (1.0, 10.0), 
:kmig2mth1=>  (1.0e-3, 1.0e3), 
:nmig2mth1=> (1.0, 10.0), 
:std1tot=> (0.1, 1.0e3), 
:istd1=> (0.1, 2.5e4), 
:nstd1=> (1.0, 10.0), 
:nstd3=> (1.0, 10.0), 
:estd1max=>(0.1, 2.5e4), 
:mig1tot=> (.9999999999,1), 
:imig1=> (0.1, 2.5e4) , 
:kmig1snf1=>  (1.0e-3, 1.0e3), 
:emig1max=>(0.1, 2.5e4), 
:dmig2=>  (0.1, 2.5e3), 
:dmig2snf1=>  (0.1, 2.5e3), 
:kmig2snf1=>   (1.0e-3, 1.0e3), 
:smig2=>(1.0e-3, 120), 
:kmig2std1=>  (1.0e-3, 1.0e3), 
:nmig2std1=>(1.0, 10.0), 
:kmig2mth1std1=>  (1.0e-3, 1.0e3), 
:nmig2mth1std1=> (1.0, 10.0), 
:dhxt4=>  (0.1, 2.5e3), 
:dhxt4max=> (0.1, 2.5e3), 
:kdhxt4=>  (1.0e-3, 1.0e3), 
:ndhxt4=> (1.0, 10.0), 
:shxt4=> (1.0e-3, 120), 
:khxt4mth1=> (1.0e-3, 1.0e3), 
:nhxt4mth1=> (1.0, 10.0), 
:khxt4std1=>  (1.0e-3, 1.0e3), 
:nhxt4std1=>(1.0, 10.0), 
:khxt4mth1std1=>  (1.0e-3, 1.0e3),
:nhxt4mth1std1=> (1.0, 10.0), #Not exactly sure about what bounds to use for initial conditions.  doing form 0 to 10 in mig1 units.
:mutk2=>(1.0e-3, 1.0e3),
:mutk3=>(1.0e-3, 1.0e3),
:khxt4mig1=>(1.0e-3, 1.0e3),
:khxt4mig2=>(1.0e-3, 1.0e3),
:nhxt4mig1=> (1.0, 10.0),
:nhxt4mig2=> (1.0, 10.0),
:Hxt4_0=> (0.0, 1.0e1), 
:Mig1_0=> (0.0, 1.0e1),  
:Mig2_0=> (0.0, 1.0e1), 
:Mth1_0=> (0.0, 1.0e1),  
:Std1_0=> (0.0, 1.0e1), 
]
#from the above universal collection of bounds, make a local collection of bounds to match the parameters in the current model

localbounds=[p[1]=> log.(getvalue(bounds, p[1]))  for p in bfpairs6];




function fitsbounds(p, bnds)

if p>= bnds[1] && p<=bnds[2]
true
else
false
end
end

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



function parimport(par1, par2, spl)

parset=OrderedDict(par1)
p=collect(parset)
for kk in spl
p[kk]= par1[kk][1]=> par2[kk][2]
end
p
end



function parrep(par1,  spl, new)

parset=OrderedDict(par1)
p=collect(parset)
for kk in spl
p[kk]= par1[kk][1]=> new
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

#parameters of the regulators template
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

##best fit parameters for model hxtmodel5d
bestpars6=[-1.99995, 1.06211, -5.13512, 4.67567, 5.01806, -1.84446, 5.73796, -2.87079, 0.0698773, 0.288116, 0.376694, 3.14651, -4.50517, -0.379322, 5.25098, 0.560657, -1.74271, -0.245076, 1.22922, 0.746876, 4.85039, 5.45447, -2.80427, -0.264936, 3.03049, 5.25449, 5.86184, 2.28341, 5.81895, 1.63477, 2.87004, 2.80064, -1.15641, -2.48345, -3.90747, -4.58572, 1.46517, -4.65211, 3.03817, 0.611965, 1.35178, 6.74873, -0.211267, -1.23877, -2.57697, 1.36639, 2.07541, -100.0, -100.0, -100.0, -4.54752, 5.57239]


bfpairs6=[pars5d[j][1]=> bestpars6[j] for j in 1:length(pars5d)]







hillrep=lambda maxv, repressor, k, n: maxv/ (1+ (repressor^n/k^n))



function hillrep(maxv=1, repressor, k, n)
maxv/ (1+ (repressor^n/k^n))
end
plot()

glucose, "k2", ,"k3", "estd1max"
glucose, k2, k3, mth1
glucose +std1 -ksnf1std1, nsnf1, ksnf2, nsnf2| snf1 

snf1 -| mig1, mig2


mig1-| hxt4

mig2-| hxt4
















for f in 1:7

parpair=[pars5d[j][1] =>allhxtpars[f][j] for j in 1:length(pars5d)]
ps=OrderedDict(maskregulators(parpair))
reprange=0:3000
plot!(hillrep.(reprange, [ps[:khxt4mth1]], [ps[:nhxt4mth1]]), label="hxt"*string(f), title="Mth1 repression kinetics", xlabel="Mth1 (relative to Mig1)", show=true)
end







parshxt1=[2.63289,  
3.43367,
-0.178529,
1.66597,
4.9781,
5.16262,
8.1147,
-0.462304,
-1.63929,
-1.18476,
-2.89898,
3.86372,
-1.0823,
0.524441,
-2.148,
1.12131,
4.86626,
-2.36051,
-2.67064,
-3.46295,
4.4916,
7.16251,
-2.08096,
-0.538442,
-5.89613,
3.97264,
-0.514723,
-1.60547,
-4.24597,
2.75386,
3.53546,
1.95889,
-1.70582,
-0.933851,
-2.64086,
1.11174,
2.02209,
-5.60155,
0.516379,
1.9637,
3.25843,
-5.43197,
-4.02793,
1.92474,
4.73574,
1.79486,
-0.329861,
-30.0,
-30.0,
-30.0,
-4.92249,
0.914295]



allhxtpars=[[2.12585, -0.228589, 2.64776, 2.28303, 4.65868, -0.377864, -0.947305, -0.585141, -4.22831, 1.0364, 2.27176, -1.44784, 2.83185, -0.479837, 5.08378, 3.96317, 3.73702, 3.87127, 0.467661, 4.19409, 1.31953, -0.0637325, 0.555823, -1.25645, 0.822131, -3.9334, 2.46758, 2.66472, -1.51326, 4.38586, -4.1655, 4.11043, -1.54677, -1.24495, -2.78571, 4.83177, 2.01964, 4.70408, 1.00814, -0.37189, 2.08752, -4.65876, -2.09565, 4.04565, -1.16734, 5.36862, 4.57054, -100.0, -100.0, -100.0, 6.36287, 1.90093, -2.64908, -1.3065], [-4.05375, -4.38865, -0.387126, -2.09917, -3.45934, -1.85104, 3.00272, -4.63336, -0.367278, -3.96581, -2.58058, -3.29509, -1.10086, -2.25288, 1.33466, -0.997993, -4.14088, 4.91687, 3.57089, 4.97971, 5.26009, 4.27863, -2.05661, 0.600067, -1.37702, -3.97399, -0.135995, -5.64302, -3.38755, 3.12066, 1.94872, 2.29846, -1.69684, -3.8538, 3.50149, 2.6055, 2.26233, 4.1614, 4.5567, 2.0077, 3.58246, -4.36608, 0.993411, -0.421662, -2.49428, 1.58728, 1.91161, -100.0, -100.0, -100.0, -0.933629, 7.16856, 1.49767, -3.59741], [1.51031, -1.60576, -0.227281, 0.720219, -3.19034, -2.88299, -1.25121, 2.88455, -0.0725552, 0.453426, -5.49732, -0.825076, -1.94652, 2.61518, -1.76637, 2.74988, -0.657427, -0.410301, -2.8652, 4.12913, 5.76514, 1.96237, -1.15051, 2.42945, 1.17839, -1.13773, 4.68304, 3.32221, -1.47287, -1.31862, -0.524745, -1.95694, -2.057, -1.87554, 3.04574, -2.81288, 2.178, 4.69048, 0.495032, -2.66179, -2.98485, -6.01788, 1.49537, 2.60623, -1.76547, 0.464955, -3.02913, -100.0, -100.0, -100.0, -0.528297, 4.85383, 1.20051, 0.143849], [-3.4195, -2.35805, 2.05974, 0.368234, -0.0946929, 1.75379, 2.27188, -0.528264, -4.40721, 3.8211, -2.12964, 0.884954, 0.241284, -4.89921, -2.921, 6.43465, -0.447138, 1.16962, 0.516744, -0.0234434, 5.2273, 1.2149, 2.12026, 4.35153, -5.12663, 4.69463, 3.66627, -0.457804, -4.91795, -2.81858, 2.2735, 4.67642, -6.10617, -1.30548, 4.06367, -2.45171, 2.83326, 3.3878, -3.4586, -2.1655, -1.47081, -5.26454, 1.93605, -0.697911, 2.53523, 4.7459, 5.2779, -100.0, -100.0, -100.0, -2.64272, 3.61168, 0.775196, -4.07627], [3.19914, -4.36493, 3.88918, 0.759819, 0.0852745, 1.07591, 2.49312, 1.87091, 4.54656, -2.96536, -0.256022, -5.10586, 2.13312, -0.893856, -1.00505, 7.27344, 3.05992, 3.4837, 2.51924, 5.28089, 4.91284, 4.91682, 1.78848, 0.213955, 2.50124, 4.30451, 1.17171, -5.45006, -3.62831, -0.0148271, 1.432, -0.358004, -1.65598, -2.90821, -1.90595, -4.64083, 3.09701, 7.33114, 0.351706, 1.30168, -0.844006, 3.59283, -0.89614, 2.11009, -3.58634, 0.0538381, 2.02441, -100.0, -100.0, -100.0, 2.71272, 1.08731, -2.51134, 4.45764], [-3.41298, -4.51073, -4.8215, -2.12027, -1.09183, 2.03998, 5.69339, -3.72588, 3.1118, 2.9018, 2.82258, -4.18625, -5.35592, -3.93559, 3.09606, -0.693105, 1.71952, 3.77457, 1.7363, 3.4572, -1.99356, 4.13231, -1.60243, 0.966916, -2.34854, -1.61179, 1.5506, -4.25438, 4.29903, 1.51959, 3.05977, 4.41772, -3.15216, -0.810818, 5.75198, -0.174855, 2.25499, 3.33385, 1.65737, 0.433048, 4.96045, -0.015301, -0.175637, 2.45086, -3.05579, 4.67977, 4.95837, -100.0, -100.0, -100.0, 1.72714, 4.45166, 5.00693, -3.93889], [0.370035, 1.33056, -5.62519, -4.57604, -2.54942, 2.89442, -0.549713, 1.81785, -3.97815, 0.301566, 2.92283, 0.00376165, -1.17438, -0.164133, -3.17186, 4.75173, -1.73025, 1.36004, 3.99965, -2.176, 2.78085, -2.08821, -2.18675, 0.331417, -4.34986, 2.34637, -0.689311, 0.823236, -0.61949, -0.492952, -2.14165, 1.67272, -4.53272, -0.846265, 4.38201, 0.473423, 2.47547, 0.15433, -0.977434, 2.98446, -0.84466, 0.276438, -1.1219, 0.824498, -3.14136, 0.00921242, 4.30019, -100.0, -100.0, -100.0, 4.29784, 1.21602, 5.79175, 1.38541]]


bfpairs6=

[

pars={
"k3":-1.99995,
"k2":1.06211,
"ksnf1":-5.13512,
"ksnf1std1":4.67567,
"nsnf1":5.01806,
"nsnf2":-1.84446,
"Snf1tot":5.73796,
"dmth1":-2.87079,
"nmth1snf3":0.0698773,
"nmth1rgt2":0.288116,
"dmth1snf3":0.376694,
"dmth1rgt2":3.14651,
"smth1":-4.50517,
"kmig1mth1":-0.379322,
"nmig1mth1":5.25098,
"kmig2mth1":0.560657,
"nmig2mth1":-1.74271,
"std1tot":-0.245076,
"istd1":1.22922,
"nstd1":0.746876,
"estd1max":4.85039,
"imig1":5.45447,
"kmig1snf1":-2.80427,
"emig1max":-0.264936,
"dmig2":3.03049,
"dmig2snf1":5.25449,
"kmig2snf1":5.86184,
"smig2":2.28341,
"kmig2std1":5.81895,
"nmig2std1":1.63477,
"kmig2mth1std1":2.87004,
"nmig2mth1std1":2.80064,
"dhxt4":-1.15641,
"dhxt4max":-2.48345,
"kdhxt4":-3.90747,
"ndhxt4":-4.58572,
"shxt4":1.46517,
"khxt4mth1":-4.65211,
"nhxt4mth1":3.03817,
"khxt4std1":0.611965,
"nhxt4std1":1.35178,
"khxt4mth1std1":6.74873,
"nhxt4mth1std1":-0.211267,
"khxt4mig1":-1.23877,
"khxt4mig2":-2.57697,
"nhxt4mig1":1.36639,
"nhxt4mig2":2.07541,
"Hxt4_0":-100.0,
"Mig1_0":-100.0,
"Mig2_0":-100.0,
"Mth1_0":-4.54752,
"Std1_0":5.57239,
"mutk2":log(10.0),
"mutk3":log(10.0)

}



templatenames=
[:k3,
:k2,
:ksnf1,
:ksnf1std1,
:nsnf1,
:nsnf2,
:Snf1tot,
:dmth1,
:nmth1snf3,
:nmth1rgt2,
:dmth1snf3,
:dmth1rgt2,
:smth1,
:kmig1mth1,
:nmig1mth1,
:kmig2mth1,
:nmig2mth1,
:std1tot,
:istd1,
:nstd1,
:estd1max,
:imig1,
:kmig1snf1,
:emig1max,
:dmig2,
:dmig2snf1,
:kmig2snf1,
:smig2,
:kmig2std1,
:nmig2std1,
:kmig2mth1std1,
:nmig2mth1std1,
:Hxt4_0,
:Mig1_0,
:Mig2_0,
:Mth1_0,
:Std1_0]
templatepars= [ j=>getvalue(bfpairs6, j) for j in templatenames] 





function cappars(pars, bounds)

for j in 1:length(pars)

parss=copy(pars)
if parss[j][2]<getvalue( bounds, parss[j][1])[1]

parss[j]=parss[j][1]=>getvalue(bounds,parss[j][1])[1]

elseif parss[j][2]>getvalue( bounds, parss[j][1])[2]

parss[j]=parss[j][1]=> getvalue(bounds,parss[j][1])[2]


end


end
parss
end


function cappar(par, bound)
#check if a parameter is within bounds and cap it to the bounds if not. 
if par<bound[1]
return bound[1]
else return par
if par>bound[2]
return bound[2] else return par
end
end
end

function partype

#names for marc's parameters. initial condition for hxt4 is lost
names=[
"k3",
"k2",
"ksnf1",
"ksnf1std1",
"nsnf1",
"nsnf2",
"Snf1tot",
"dmth1",
"nmth1snf3",
"nmth1rgt2",
"dmth1snf3",
"dmth1rgt2",
"smth1",
"kmig1mth1",
"nmig1mth1",
"kmig2mth1",
"nmig2mth1",
"std1tot",
"istd1",
"nstd1",
"estd1max",
"imig1",
"kmig1snf1",
"emig1max",
"dmig2",
"dmig2snf1",
"kmig2snf1",
"smig2",
"kmig2std1",
"nmig2std1",
"kmig2mth1std1",
"nmig2mth1std1",
"dhxt4",
"dhxt4max",
"kdhxt4",
"ndhxt4",
"shxt4",
"khxt4mth1",
"nhxt4mth1",
"khxt4std1",
"nhxt4std1",
"khxt4mth1std1",
"nhxt4mth1std1",
"khxt4mig1",
"khxt4mig2",
"nhxt4mig1",
"nhxt4mig2",
#"Hxt4_0", hxt4 initial condition is not present in marcs parameters
"Mig1_0",
"Mig2_0",
"Mth1_0",
"Std1_0",
"mutk2",
"mutk3"
]


###MArc's best fit. 320000 simulations with the mean of the standard deviations


d2 = [-1.60724, -1.11798, -7.70941, 9.91471, 7.07934, -0.0734542, 8.70169, -1.59598, 5.18232, 0.551847, 3.9686, 1.6702, -4.23113, 0.680474, 4.64808, 7.48488, 3.65958, 4.35926, -3.53894, 0.412305, -0.707678, 6.1692, 3.85898, 6.94579, 0.993898, -4.25699, -4.12453, 4.08724, 1.81226, 2.52904, -3.75305, -4.81768, -2.88543, -1.25432, 4.27927, 4.98911, 6.1547, -5.69309, 1.51019, -0.620259, 0.491589, -2.32399, 2.23031, -4.30685, -2.41038, 0.63114, 0.368315, 7.84135, 9.99994, 8.05267, -4.72919, 0.517147, -2.11891]

marcpars=ppf.colorDict(names, d2)      
marcpars_sd={'k3': -1.60724, 'k2': -1.11798, 'ksnf1': -7.70941, 'ksnf1std1': 9.91471, 'nsnf1': 7.07934, 'nsnf2': -0.0734542, 'Snf1tot': 8.70169, 'dmth1': -1.59598, 'nmth1snf3': 5.18232, 'nmth1rgt2': 0.551847, 'dmth1snf3': 3.9686, 'dmth1rgt2': 1.6702, 'smth1': -4.23113, 'kmig1mth1': 0.680474, 'nmig1mth1': 4.64808, 'kmig2mth1': 7.48488, 'nmig2mth1': 3.65958, 'std1tot': 4.35926, 'istd1': -3.53894, 'nstd1': 0.412305, 'estd1max': -0.707678, 'imig1': 6.1692, 'kmig1snf1': 3.85898, 'emig1max': 6.94579, 'dmig2': 0.993898, 'dmig2snf1': -4.25699, 'kmig2snf1': -4.12453, 'smig2': 4.08724, 'kmig2std1': 1.81226, 'nmig2std1': 2.52904, 'kmig2mth1std1': -3.75305, 'nmig2mth1std1': -4.81768, 'dhxt4': -2.88543, 'dhxt4max': -1.25432, 'kdhxt4': 4.27927, 'ndhxt4': 4.98911, 'shxt4': 6.1547, 'khxt4mth1': -5.69309, 'nhxt4mth1': 1.51019, 'khxt4std1': -0.620259, 'nhxt4std1': 0.491589, 'khxt4mth1std1': -2.32399, 'nhxt4mth1std1': 2.23031, 'khxt4mig1': -4.30685, 'khxt4mig2': -2.41038, 'nhxt4mig1': 0.63114, 'nhxt4mig2': 0.368315, 'Mig1_0': 7.84135, 'Mig2_0': 9.99994, 'Mth1_0': 8.05267, 'Std1_0': -4.72919, 'mutk2': 0.517147, 'mutk3': -2.11891}


glucose2=[2.0,
 1.9973364358573156,
 1.9962854493829298,
 1.9945856093790901,
 1.9929642587446645,
 1.9918278535304812,
 1.988329747859795,
 1.9850648941782956,
 1.9831827074549198,
 1.9788314443606538,
 1.9728047213638116,
 1.9684869605160598,
 1.9618860478373632,
 1.9546412659421979,
 1.9469908777426785,
 1.9372819807774435,
 1.925614330995987,
 1.9162628520658613,
 1.9015789292045251,
 1.8882855454687553,
 1.8699292151655635,
 1.8513036618491099,
 1.835829189163912,
 1.8123774311441649,
 1.7879240499474358,
 1.7641448443097489,
 1.7384972799513285,
 1.7090858056294143,
 1.6758876522660948,
 1.6394602954411734,
 1.6035251570170175,
 1.5656636369173782,
 1.5282878689788146,
 1.4908755677181449,
 1.4513543332145438,
 1.408299095699638,
 1.3644123191843276,
 1.3142091037264065,
 1.2616067905404136,
 1.1975923893273555,
 1.1287741261382296,
 1.0545482528413546,
 0.97280042216779217,
 0.89122846772703568,
 0.80267200701233254,
 0.71568089762040454,
 0.62102891914447533,
 0.53055771883930802,
 0.44238492604195323,
 0.35394307973641781,
 0.27643896324231076,
 0.20403708542392507,
 0.14856078718325816,
 0.10457613776171515,
 0.07292179298941015,
 0.051584967711607144,
 0.031629122330370629,
 0.020279572892611286,
 0.018560013677649145,
 0.009621834654893302,
 0.0054977583436375266,
 0.0082467818707816498,
 0.0041229723362004034,
 0.0068722256619822009,
 0.0044666601941891759,
 0.0037792902944844275,
 -0.00068721327591525494,
 -0.0027482922017423306,
 0.0]


glucose2time=[  0.        ,   0.209     ,   0.41788889,   0.62675   ,
         0.83558333,   1.04444444,   1.25330556,   1.46213889,
         1.671     ,   1.87986111,   2.08872222,   2.29755556,
         2.50641667,   2.71527778,   2.92411111,   3.13297222,
         3.34183333,   3.55069444,   3.75952778,   3.96838889,
         4.17725   ,   4.38611111,   4.59494444,   4.80380556,
         5.01266667,   5.2215    ,   5.43036111,   5.63922222,
         5.84805556,   6.05691667,   6.26577778,   6.47461111,
         6.68347222,   6.89233333,   7.10119444,   7.31005556,
         7.51888889,   7.72775   ,   7.93661111,   8.14544444,
         8.35430556,   8.56316667,   8.772     ,   8.98086111,
         9.18972222,   9.39858333,   9.60741667,   9.81627778,
        10.02513889,  10.23397222,  10.44283333,  10.65169444,
        10.86055556,  11.06938889,  11.27825   ,  11.48711111,
        11.69597222,  11.90480556,  12.11366667,  12.32252778,
        12.53136111,  12.74022222,  12.94908333,  13.15794444,
        13.36677778,  13.57580556,  13.7845    ,  13.99333333,  14.20219444]


