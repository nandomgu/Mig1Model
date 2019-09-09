using Distributed #- can parallelise code here by adding extra cores...
@everywhere begin
    using StatsBase, KernelDensity, Random
    using Distances, Distributions, DataStructures
    import Distributions.rand
    import Distributions.pdf
    using LsqFit, JSON
    include("/home/msturrock/Desktop/Mig1/Mig1Model/mechanisticModel/sencillo_model_comparison/driver.jl")
    glucose=[8.32733959421089e-19,2.561481685921e-19,1.04197119608466e-19,5.53738039368354e-20,2.30389033029131e-20,1.01730601065369e-20,9.64607628015416e-21,4.54189858719156e-21,2.62590228463374e-21,1.07030601256657e-21,7.12893122193276e-22,3.77665654101745e-22,1.13613450886179e-22,9.11107194068206e-23,1.68345898188741e-23,4.59571704388656e-23,2.4737145794224e-23,5.06371909794061e-24,5.65818382020064e-24,1.04048279473923e-24,3.94925490074619e-24,3.45909892619726e-25,7.07762515569138e-25,7.59557441484042e-26,2.02035130943482e-27,3.27268072550391e-26,1.56537448533959e-27,1.6320588929363e-27,4.39103892339638e-28,4.80386219709977e-29,3.78228642078199e-28,6.7008774844888e-33,3.38770680654291e-30,3.11390728411053e-36,3.43524261389185e-36,9.0953685838273e-37,0,2.65814991285393e-16,1.40534785861741e-05,0.956449369761021,0.999932720540147,0.999988138202835,0.999994336997356,0.999995688515055,0.99999655799382,0.999997230933311,0.999997349602989,0.999997796320276,0.999997725777247,0.999997825419588,0.999997959294029,0.99999807290326,0.999998148682981,0.999998097810728,0.999998115954266,0.999998085169392,0.999998126522822,0.99999814575262,0.999998172022714,0.999998338656721,0.999998383673738,0.999998504114664,0.999998515535666,0.999998447595621,0.999998528661984,0.99999855188616,0.999998519653104,0.999998565797203,0.999998565177499,0.999998598958849,0.999998669110336,0.999998719774422,0.999998776796943,0.999998794501031,0.99999883534008,0.999998798817409,0.99999880270023,0.999998808299886,0.999998809973147,0.999998846846304,0.999998882337051,0.999998957349787,0.999998958954537,0.999998984982437,0.999998977466028,0.999998955080815,0.999998965848435,0.999998965848435,0.999998955080815,0.999998977466028,0.999998984982437,0.999998958954537,0.999998957349787,0.999998882337051,0.999998846846304,0.999998809973147,0.999998808299886,0.99999880270023,0.999998798817409,0.99999883534008,0.999998794501031,0.999998776796943,0.999998719774422,0.999998669110336,0.999998598958849,0.999998565177499,0.999998565797203,0.999998519653104,0.99999855188616,0.999998528661984,0.999998447595621,0.999998515535666,0.999998504114664,0.999998383673738,0.999998338656721,0.999998172022714,0.99999814575262,0.999998126522822,0.999998085169392,0.999998115954266,0.999998097810728,0.999998148682981,0.99999807290326,0.999997959294029,0.999997825419588,0.999997725777247,0.999997796320276,0.999997349602989,0.999997230933311,0.99999655799382,0.999995688515055,0.999994336997356,0.999988138202835,0.999932720540147,0.956449369761021,1.40534785861741e-05,2.65814991285393e-16,0,9.0953685838273e-37,3.43524261389185e-36,3.11390728411053e-36,3.38770680654291e-30,6.7008774844888e-33,3.78228642078199e-28,4.80386219709977e-29,4.39103892339638e-28,1.6320588929363e-27,1.56537448533959e-27,3.27268072550391e-26,2.02035130943482e-27,7.59557441484042e-26,7.07762515569138e-25,3.45909892619726e-25,3.94925490074619e-24,1.04048279473923e-24,5.65818382020064e-24,5.06371909794061e-24,2.4737145794224e-23,4.59571704388656e-23,1.68345898188741e-23,9.11107194068206e-23,1.13613450886179e-22,3.77665654101745e-22,7.12893122193276e-22,1.07030601256657e-21,2.62590228463374e-21,4.54189858719156e-21,9.64607628015416e-21,1.01730601065369e-20,2.30389033029131e-20,5.53738039368354e-20,1.04197119608466e-19,2.561481685921e-19,8.32733959421089e-19,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
    y=glucose[1:229]
    datameans=JSON.parsefile("/home/msturrock/Desktop/Mig1/Mig1Model/json/allfitmeans.json", dicttype=Dict)
    dm2=datameans
    subsd(x)=  if typeof(x)==String return 1000 else return x end # if it is a string(missing value) return a st of 1000  to minimise the relevance of this point.
    subzero(x)=  if x==0 return 1 else return x end 
    meanerr=JSON.parsefile("/home/msturrock/Desktop/Mig1/Mig1Model/json/meanerr.json", dicttype=Dict)
    me2=[subzero.(subsd.(j)) for j in meanerr]
    stdmeans=JSON.parsefile("/home/msturrock/Desktop/Mig1/Mig1Model/json/stdmeans.json", dicttype=Dict)
    sm2=[subzero.(subsd.(j)) for j in stdmeans]
    global interv= 0.044543429844097995; 
    nsteps=convert(Int64, ceil(18.39/interv))
    tims= repeat([interv], nsteps)
    x=cumsum(tims) 
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
    tt=(0.0, 18.39) #this time is roughly the last timepoint of x
    genotypes= repeat([wt,mth1ko, mig1ko, std1ko, rgt2ko, snf3ko], inner=3)
    concs=repeat([0.2, 0.4, 1], 6) # glucose concentrations  in the order given by data.
    modelfile="/home/msturrock/Desktop/Mig1/Mig1Model/mechanisticModel/sencillo_model_comparison/hxt4model3d.jl"

    initial_conds = [splines[j](x)[1] for j in 1:18]
    initial_conds[initial_conds.<0.0] .= 0.0
    allprobs=makeproblem.([modelfile], [tim[1:229]], [y], concs, [tt],initial_conds, genotypes)
    function trysolve(prob, pars, datamean, datastd) #main function to solve a system with two solvers and return a large value if both fail. 

        try
               sol = solve(prob(pars), Rosenbrock32(autodiff=false), d_discontinuities=[3.16,3.25,3.33,11.25,11.33,11.42], saveat = interv,verbose=false)
               arr = [[j[i] for j in sol.u] for i = 1:length(sol.u[1])]
               lsqval = sum(((arr[1] - datamean) / datastd).^2)
        catch
                  try
                         sol = solve(prob(pars), Rodas5(autodiff=false), d_discontinuities=[3.16,3.25,3.33,11.25,11.33,11.42], saveat = interv,verbose=false)
                         arr = [[j[i] for j in sol.u] for i = 1:length(sol.u[1])]
                         lsqval = sum(((arr[1] - datamean) / datastd).^2)
                  catch
                         lsqval = 10000.0
                  end
        end
    end

    include("/home/msturrock/Desktop/Mig1/Mig1Model/mechanisticModel/sencillo_model_comparison/types.jl")
    include("/home/msturrock/Desktop/Mig1/Mig1Model/mechanisticModel/sencillo_model_comparison/functions.jl")
    include("/home/msturrock/Desktop/Mig1/Mig1Model/mechanisticModel/sencillo_model_comparison/abc2.jl")


    function rho_lens(expd,d2)
        error = sum(trysolve.(allprobs,[d2], [j(x) for j in splines],me2))
        if isnan(error) || isinf(error) 
            error = 10000.0
        end
        return error
    end

end

bounds=[
    :k3=>[1.0e-3, 1.0e3],
    :k2=> [1.0e-3, 1.0e3],
    :ksnf1=>[1.0e-3, 1.0e3],
    :ksnf1std1=> [1.0e-3, 1.0e3],
    :nsnf1=>[1.0, 10.0],
    :nsnf2=>[1.0, 10.0],
    :Snf1tot=> [0.1, 1.0e2], #10 times less than mig1 or 100 times more.
    :dmth1=>  [0.1, 2.5e3],
    :nmth1snf3=> [1.0, 10.0],
    :nmth1rgt2=> [1.0, 10.0],
    :dmth1snf3=> [0.1, 2.5e3],
    :dmth1rgt2=>  [0.1, 2.5e3], 
    :smth1=> [1.0e-3, 120],
    :kmig1mth1=>  [1.0e-3, 1.0e3],
    :nmig1mth1=> [1.0, 10.0],
    :kmig2mth1=>  [1.0e-3, 1.0e3],
    :nmig2mth1=> [1.0, 10.0],
    :std1tot=> [0.1, 1.0e3],
    :istd1=> [0.1, 2.5e4],
    :nstd1=> [1.0, 10.0],
    :estd1max=>[0.1, 2.5e4],
    :mig1tot=> [.9999999999,1], 
    :imig1=> [0.1, 2.5e4] ,
    :kmig1snf1=>  [1.0e-3, 1.0e3],
    :emig1max=>[0.1, 2.5e4],
    :dmig2=>  [0.1, 2.5e3],
    :dmig2snf1=>  [0.1, 2.5e3],
    :kmig2snf1=>   [1.0e-3, 1.0e3],
    :smig2=>[1.0e-3, 120],
    :kmig2std1=>  [1.0e-3, 1.0e3],
    :nmig2std1=>[1.0, 10.0],
    :kmig2mth1std1=>  [1.0e-3, 1.0e3],
    :nmig2mth1std1=> [1.0, 10.0],
    :dhxt4=>  [0.1, 2.5e3],
    :dhxt4max=> [0.1, 2.5e3],
    :kdhxt4=>  [1.0e-3, 1.0e3],
    :ndhxt4=> [1.0, 10.0],
    :shxt4=> [1.0e-3, 120], 
    :khxt4mth1=> [1.0e-3, 1.0e3],
    :nhxt4mth1=> [1.0, 10.0],
    :khxt4std1=>  [1.0e-3, 1.0e3], 
    :nhxt4std1=>[1.0, 10.0],
    :khxt4mth1std1=>  [1.0e-3, 1.0e3],
    :nhxt4mth1std1=> [1.0, 10.0]]
    #Not exactly sure about what bounds to use for initial conditions.  doing form 0 to 10 in mig1 units.
    bounds2 = [:Mig1_0=> [0.0, 1.0e1],  
    :Mig2_0=> [0.0, 1.0e1], 
    :Mth1_0=> [0.0, 1.0e1],  
    :Std1_0=> [0.0, 1.0e1]]

pardistributions = [[Uniform(log(p[2][1]), log(p[2][2])) for p in bounds]; [Uniform(p2[2][1],p2[2][2]) for p2 in bounds2]]

np = 1000 # change the number of particles

@time apmc_output =APMC_KDE(np,0.0,[pardistributions],[rho_lens],paccmin=0.01)

#plot best parameter set
d2 = apmc_output.pts[end][1:48,1]
d3 = apmc_output.pts[end][1:48,:]

using Plots
d2 = rand(pardistributions)

inspectdr(grid=false,linewidth=3.0,ylabel="HXT4",labels=["model" "data"])
p=[]
for i = 1:18
    sol = solve(allprobs[i](d2), Rosenbrock32(autodiff=false), d_discontinuities=[3.16,3.25,3.33,11.25,11.33,11.42], saveat = interv,verbose=false)
    push!(p,plot(x,[sol[1,:] splines[i](x)],xlabel="time",title="condition $i"))
end

plot(p[1:18]...,layout=(6,3),size=(1000,1000))
savefig("sampleoutput_hxt4.png")

inspectdr(grid=false,linewidth=3.0,ylabel="mig1")
p=[]
for i = 1:18
    sol = solve(allprobs[i](d2),Rosenbrock32(autodiff=false),d_discontinuities=[3.16,3.25,3.33,11.25,11.33,11.42], saveat = interv,verbose=false)
    push!(p,plot(x,sol[2,:],label=["model" "data"],xlabel="time",title="condition $i"))
end

plot(p[1:18]...,layout=(6,3),size=(1000,1000))
savefig("sampleoutput_mig1.png")

gr(ylabel="frequency",label="",xlabel="parameter value")
anim = @animate for i = 1:48
    if i < 45
        histogram(d3[i,:],title=string(bounds[i][1]),xlims=(log(bounds[i][2][1]),log(bounds[i][2][2])))
    else
        histogram(d3[i,:],title=string(bounds2[i-44][1]),xlims=(bounds2[i-44][2][1],bounds2[i-44][2][2]))
    end
end
gif(anim,"posteriors.gif",fps=1)



out=10000 .* ones(100)
pars = []
for i = 1:100
    @show i
    parsample = rand(pardistributions)
    out[i]=sum(trysolve.(allprobs,[parsample], [j(x) for j in splines],me2))
    if out[i] <= minimum(out)
       global pars = parsample
    end
end


d2 = pars
using Plots
inspectdr(grid=false,linewidth=3.0,ylabel="HXT4")
p=[]
for i = 1:18
    @show i
    sol = solve(allprobs[i](d2), alg_hints=[:stiff,:interpolant], tstops=[3.25,3.33,11.33,11.42], saveat = interv,verbose=true,isoutofdomain=(y,p,t)->any(x->x<0,y),abstol=1e-8,reltol=1e-8)
    @show size(sol[1,:])
    push!(p,plot([sol[1,:] splines[i](x)],label=["model" "data"],xlabel="time",title="condition $i"))
end
plot(p[1:18]...,layout=(6,3),size=(1000,1000))
savefig("sampleoutput.png")
