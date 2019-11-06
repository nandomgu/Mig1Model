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
    tim=[0 0.0833333333333333 0.166666666666667 0.25 0.333333333333333 0.416666666666667 0.5 0.583333333333333 0.666666666666667 0.75 0.833333333333333 0.916666666666667 1 1.08333333333333 1.16666666666667 1.25 1.33333333333333 1.41666666666667 1.5 1.58333333333333 1.66666666666667 1.75 1.83333333333333 1.91666666666667 2 2.08333333333333 2.16666666666667 2.25 2.33333333333333 2.41666666666667 2.5 2.58333333333333 2.66666666666667 2.75 2.83333333333333 2.91666666666667 3 3.08333333333333 3.16666666666667 3.25 3.33333333333333 3.41666666666667 3.5 3.58333333333333 3.66666666666667 3.75 3.83333333333333 3.91666666666667 4 4.08333333333333 4.16666666666667 4.25 4.33333333333333 4.41666666666667 4.5 4.58333333333333 4.66666666666667 4.75 4.83333333333333 4.91666666666667 5 5.08333333333333 5.16666666666667 5.25 5.33333333333333 5.41666666666667 5.5 5.58333333333333 5.66666666666667 5.75 5.83333333333333 5.91666666666667 6 6.08333333333333 6.16666666666667 6.25 6.33333333333333 6.41666666666667 6.5 6.58333333333333 6.66666666666667 6.75 6.83333333333333 6.91666666666667 7 7.08333333333333 7.16666666666667 7.25 7.33333333333333 7.41666666666667 7.5 7.58333333333333 7.66666666666667 7.75 7.83333333333333 7.91666666666667 8 8.08333333333333 8.16666666666667 8.25 8.33333333333333 8.41666666666667 8.5 8.58333333333333 8.66666666666667 8.75 8.83333333333333 8.916666666666679 9.08333333333333 9.16666666666667 9.25 9.33333333333333 9.41666666666667 9.5 9.58333333333333 9.66666666666667 9.75 9.83333333333333 9.91666666666667 10 10.0833333333333 10.1666666666667 10.25 10.3333333333333 10.4166666666667 10.5 10.5833333333333 10.6666666666667 10.75 10.8333333333333 10.9166666666667 11 11.0833333333333 11.1666666666667 11.25 11.3333333333333 11.4166666666667 11.5 11.5833333333333 11.6666666666667 11.75 11.8333333333333 11.9166666666667 12 12.0833333333333 12.1666666666667 12.25 12.3333333333333 12.4166666666667 12.5 12.5833333333333 12.6666666666667 12.75 12.8333333333333 12.9166666666667 13 13.0833333333333 13.1666666666667 13.25 13.3333333333333 13.4166666666667 13.5 13.5833333333333 13.6666666666667 13.75 13.8333333333333 13.9166666666667 14 14.0833333333333 14.1666666666667 14.25 14.3333333333333 14.4166666666667 14.5 14.5833333333333 14.6666666666667 14.75 14.8333333333333 14.9166666666667 15 15.0833333333333 15.1666666666667 15.25 15.3333333333333 15.4166666666667 15.5 15.5833333333333 15.6666666666667 15.75 15.8333333333333 15.9166666666667 16 16.0833333333333 16.1666666666667 16.25 16.3333333333333 16.4166666666667 16.5 16.5833333333333 16.6666666666667 16.75 16.8333333333333 16.9166666666667 17 17.0833333333333 17.1666666666667 17.25 17.3333333333333 17.4166666666667 17.5 17.5833333333333 17.6666666666667 17.75 17.8333333333333 17.9166666666667 18 18.0833333333333 18.1666666666667 18.25 18.3333333333333 18.4166666666667 18.5 18.5833333333333 18.6666666666667 18.75 18.8333333333333 18.9166666666667 19 19.0833333333333]
    #maing sure the shape of the above is right.
    tim=tim[1:229]
    datameans=JSON.parsefile("/home/msturrock/Desktop/Mig1/Mig1Model/json/allfitmeans.json", dicttype=Dict)
    dm2=datameans
    subsd(x)=  if typeof(x)==String return 5 else return x end # if it is a string(missing value) return a st of 1000  to minimise the relevance of this point.
    subzero(x)=  if x==0 return 1 else return x end 
    meanerr=JSON.parsefile("/home/msturrock/Desktop/Mig1/Mig1Model/json/meanerr.json", dicttype=Dict)
    me2=[subzero.(subsd.(j)) for j in meanerr]

    stdmeans=JSON.parsefile("json/stdmeans.json", dicttype=Dict)#doesn't need ordering 

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
    using Sundials
    function trysolve(prob, pars, datamean, datastd) #main function to solve a system with two solvers and return a large value if both fail. 

        try
               sol = solve(prob(pars), CVODE_BDF(), saveat = interv,verbose=false,abstol=1e-10,reltol=1e-10,maxiters=1e3)
               arr = [[j[i] for j in sol.u] for i = 1:length(sol.u[1])]
               lsqval = sum(((arr[1] - datamean) / datastd).^2)
        catch
                lsqval = 1000.0
        end
    end

    include("/home/msturrock/Desktop/Mig1/Mig1Model/mechanisticModel/sencillo_model_comparison/types.jl")
    include("/home/msturrock/Desktop/Mig1/Mig1Model/mechanisticModel/sencillo_model_comparison/functions.jl")
    include("/home/msturrock/Desktop/Mig1/Mig1Model/mechanisticModel/sencillo_model_comparison/abc2.jl")

    function rho_lens(expd,d2)
        error = sum(trysolve.(allprobs,[d2], [j(x) for j in splines],[k[1:length(x)] for k in sm2]))
        if isnan(error) || isinf(error) 
            error = 10000.0
        end
        return error
    end

end

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
(-10.0,3.0),
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
(-2.0,2.0), #mth10
(-2.0,2.0), #std10
(-5.0,5.0)] #mutk3]#mutk2

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

model_lens = [Uniform(p[1],p[2]) for p in custombounds]
priors = [2.64424, 1.23356, -7.98945, 3.10794, 3.04429, -1.00647, 5.73373, 0.771657, -3.8536, 0.777529, 0.250743, 5.10096, -4.1397, -0.422502, 4.40694, 3.21153, 1.88122, -0.100531, 0.729806, 0.850592, 1.47274, 4.87223, 7.13754, 1.825, 5.53788, 0.474569, -1.59972, 0.14398, 1.04363, -1.83138, -0.0639927, -5.4498, 4.2669, -2.40372, -1.30821, 4.40911, -2.51329, 2.23823, -5.06447, 6.09821, -0.122023, 4.12318, -6.20391, -1.83399, -1.73008, -3.10607, 1.23147, 0.424585, 1.89543, 0.122836, 1.45363]

model_lens = [Uniform(p-0.5,p+0.5) for p in priors]


np = 5000 # change the number of particles
@time apmc_output =APMC_KDE(np,0.0,[model_lens],[rho_lens],paccmin=0.01)

#plot best parameter set
d2 = apmc_output.pts[end][1:51,1]
d3 = apmc_output.pts[end][1:51,:]
d2 = [-1.99995, 1.06211, -5.13512, 4.67567, 5.01806, -1.84446, 5.73796, -2.87079, 0.0698773, 0.288116, 0.376694, 3.14651, -4.50517, -0.379322, 5.25098, 0.560657, -1.74271, -0.245076, 1.22922, 0.746876, 4.85039, 5.45447, -2.80427, -0.264936, 3.03049, 5.25449, 5.86184, 2.28341, 5.81895, 1.63477, 2.87004, 2.80064, -1.15641, -2.48345, -3.90747, -4.58572, 1.46517, -4.65211, 3.03817, 0.611965, 1.35178, 6.74873, -0.211267, -1.23877, -2.57697, 1.36639, 2.07541, -100.0, -100.0, -4.54752, 5.57239]
d2 = zeros(51)
for i in 1:51
    d2[i] = mode(apmc_output.pts[end][i,1])
end

sens = inv(StatsBase.cov(transpose(exp.(d3[1:51,1:1000]))))
test=log.(diag(sens))

labels =             
["k3"
"k2"
"ksnf1"
"ksnf1std1"
"nsnf1"
"nsnf2"
"Snf1tot"
"dmth1"
"nmth1snf3"
"nmth1rgt2"
"dmth1snf3"
"dmth1rgt2"
"smth1"
"kmig1mth1"
"nmig1mth1"
"kmig2mth1"
"nmig2mth1"
"std1tot"
"istd1"
"nstd1"
"nstd3"
"estd1max"
"imig1"
"kmig1snf1"
"emig1max"
"dmig2"
"dmig2snf1"
"kmig2snf1"
"smig2"
"kmig2std1"
"nmig2std1"
"kmig2mth1std1"
"nmig2mth1std1"
"dhxt4"
"dhxt4max"
"kdhxt4"
"ndhxt4"
"shxt4"
"khxt4mth1"
"nhxt4mth1"
"khxt4std1"
"nhxt4std1"
"khxt4mth1std1"
"nhxt4mth1std1"
"khxt4mig1"
"khxt4mig2"
"nhxt4mig1"
"nhxt4mig2"
"mth10"
"Std10"
"c"]

using Plots
d2 = rand(pardistributions)


gr(ylabel="frequency",label="",xlabel="parameter value")
anim = @animate for i = 1:51
    histogram(d3[i,:],title="$(labels[i]), sensitivity rank = $(test[i])",xlims=(model_lens[i].a,model_lens[i].b))
end
gif(anim,"posteriors.gif",fps=1)



out=10000 .* ones(100)
pars = []
for i = 1:100
    @show i
    d2 = rand(model_lens)
    out[i]=sum(trysolve.(allprobs,[d2], [j(x) for j in splines],me2))
    if out[i] <= minimum(out)
       global pars = d2
    end
end
findall(out.<10000)

d2 = rand(pardistributions)
using Plots
gr(grid=false,linewidth=3.0,ylabel="HXT4",label=["model" "data"])
p=[]
for i = 1:18
    @show i
    sol = solve(allprobs[i](d2), CVODE_BDF(), saveat = interv,verbose=true,abstol=1e-12,reltol=1e-12)
    push!(p,plot([sol[1,:] splines[i](x)],ribbon=[zeros(length(x)) [k for k in sm2[i][1:length(x)]]],fillalpha=.5,label="",xlabel="time",title="condition $i"))
end
plot(p[1:18]...,layout=(6,3),size=(1000,1000))
savefig("sampleoutput_std.png")


function glucose2(t,a,b,c,n)
    return a.*exp.((-(t.-b).^n)./n*c.^n)
end

gr(grid=false,linewidth=3.0,ylabel="mig1",yscale=:linear)
p=[]
for i = 1:18
     sol = solve(allprobs[i](d2), CVODE_BDF(), saveat = interv,verbose=true,abstol=1e-12,reltol=1e-12)
     push!(p,plot([sol[2,:]],label="",xlabel="time",title="condition $i"))
end
plot(p[1:18]...,layout=(6,3),size=(1000,1000))
savefig("sampleoutputmig1.png")

gr(grid=false,linewidth=3.0,ylabel="mig2",yscale=:linear)
p=[]
for i = 1:18
     sol = solve(allprobs[i](d2), CVODE_BDF(), saveat = interv,verbose=true,abstol=1e-15,reltol=1e-15)
     push!(p,plot([sol[3,:]],label="",xlabel="time",title="condition $i"))
end
plot(p[1:18]...,layout=(6,3),size=(1000,1000))
savefig("sampleoutputmig2.png")

gr(grid=false,linewidth=3.0,ylabel="mth1",yscale=:linear)
p=[]
for i = 1:18
     sol = solve(allprobs[i](d2), CVODE_BDF(), saveat = interv,verbose=true,abstol=1e-15,reltol=1e-15)
     push!(p,plot([sol[4,:]],label="",xlabel="time",title="condition $i"))
end
plot(p[1:18]...,layout=(6,3),size=(1000,1000))
savefig("sampleoutputmth1.png")

gr(grid=false,linewidth=3.0,ylabel="Std1",yscale=:linear)
p=[]
for i = 1:18
     sol = solve(allprobs[i](d2), CVODE_BDF(), saveat = interv,verbose=true,abstol=1e-15,reltol=1e-15)
     push!(p,plot([sol[5,:]],label="",xlabel="time",title="condition $i"))
end
plot(p[1:18]...,layout=(6,3),size=(1000,1000))
savefig("sampleoutputStd1.png")


gr(ylabel="frequency",label="",xlabel="parameter value")
anim = @animate for i = 1:48
    histogram(d3[i,:],title="parameter $i",xlims=(-10,5))
end
gif(anim,"posteriors.gif",fps=1)


@everywhere using BlackBoxOptim
@everywhere function solveall(d2)
    error = sum(trysolve.(allprobs,[d2], [j(x) for j in splines],[k[1:length(x)] for k in sm2]))
    if isnan(error) || isinf(error) 
        error = 10000.0
    end
    return error
end
opts7=bbsetup(solveall; Method=:xnes, SearchRange = custombounds, NumDimensions = 51, MaxSteps = 20000,Workers = workers())
res7=bboptimize(opts7, MaxSteps=40000)

d2 = [-0.955097, -0.402872, -7.99571, 3.66697, 4.19597, 0.300486, 7.3091, -5.58976, -4.1055, 0.646965, 1.0733, 2.53052, -0.00222343, -0.0899301, 5.26236, 7.48491, -4.73763, 0.444363, 0.607859, 0.90996, -8.51219, 1.69052, 2.11599, 1.46336, 3.85118, 1.52119, -0.455843, 3.68919, 0.428872, -0.718575, 0.826405, -5.70954, -1.14255, -1.34946, -2.30892, -1.54427, 2.18582, 2.23409, -1.00402, 2.56725, 0.395004, 3.31004, 2.27537, -4.81768, -4.02057, -3.58627, -0.0918299, 0.93555, 1.54937, 0.753939, 3.30952]

d2 = [2.64424, 1.23356, -7.98945, 3.10794, 3.04429, -1.00647, 5.73373, 0.771657, -3.8536, 0.777529, 0.250743, 5.10096, -4.1397, -0.422502, 4.40694, 3.21153, 1.88122, -0.100531, 0.729806, 0.850592, 1.47274, 4.87223, 7.13754, 1.825, 5.53788, 0.474569, -1.59972, 0.14398, 1.04363, -1.83138, -0.0639927, -5.4498, 4.2669, -2.40372, -1.30821, 4.40911, -2.51329, 2.23823, -5.06447, 6.09821, -0.122023, 4.12318, -6.20391, -1.83399, -1.73008, -3.10607, 1.23147, 0.424585, 1.89543, 0.122836, 1.45363]


using DiffEqSensitivity
m = gsa(allprobs[1](d2), CVODE_BDF(),tim,[[exp(genotypes[1](d2)[i]-0.1),exp(genotypes[1](d2)[i]+0.1)] for i = 1:51],Morris(p_steps=10 .* ones(51),len_trajectory=1500,total_num_trajectory=1000,num_trajectory=150),1e-15,1e-15)

using Plots
i=51
Plots.plot(m.means[i][1,:], ribbon=sqrt.(m.variances[i])[1,:],title="$(labels[i])",label="",xlabel="time",ylabel="HXT4")

gr(ylabel="change in HXT4",label="",xlabel="time")
anim = @animate for i = 1:51
    Plots.plot(m.means[i][1,:], ribbon=sqrt.(m.variances[i])[1,:],title="$(labels[i])")
end
gif(anim,"sens_morris_HXT4.gif",fps=1)

gr(ylabel="change in Mig1",label="",xlabel="time")
anim = @animate for i = 1:51
    Plots.plot(m.means[i][2,:], ribbon=sqrt.(m.variances[i])[2,:],title="$(labels[i])")
end
gif(anim,"sens_morris_mig1.gif",fps=1)

gr(ylabel="change in Mig2",label="",xlabel="time")
anim = @animate for i = 1:51
    Plots.plot(m.means[i][3,:], ribbon=sqrt.(m.variances[i])[3,:],title="$(labels[i])")
end
gif(anim,"sens_morris_mig2.gif",fps=1)

gr(ylabel="change in Mth1",label="",xlabel="time")
anim = @animate for i = 1:51
    Plots.plot(m.means[i][4,:], ribbon=sqrt.(m.variances[i])[4,:],title="$(labels[i])")
end
gif(anim,"sens_morris_mth1.gif",fps=1)

gr(ylabel="change in Std1",label="",xlabel="time")
anim = @animate for i = 1:51
    Plots.plot(m.means[i][5,:], ribbon=sqrt.(m.variances[i])[5,:],title="$(labels[i])")
end
gif(anim,"sens_morris_std1.gif",fps=1)

using Plots
for kk = 4:18
    m = gsa(allprobs[kk](d2), CVODE_BDF(),tim,[[exp(genotypes[kk](d2)[ii]-0.1),max(exp(genotypes[kk](d2)[ii]+0.1),1e-20)] for ii = 1:48],Morris(p_steps=10 .* ones(51),len_trajectory=1500,total_num_trajectory=1000,num_trajectory=150),1e-15,1e-15)
    gr(ylabel="variance in change",label="",xlabel="mean change",title="condition = $kk")
    mean_hxt = zeros(48)
    var_hxt = zeros(48)
    for i = 1:48
        mean_hxt[i] = sum(m.means[i][1,:])
        var_hxt[i] = sum(m.variances[i][1,:])
    end
    scatter(log.(mean_hxt),log.(var_hxt), series_annotations = text.(labels[1:48], :bottom,6))
    savefig("mean_vs_variance_condition_$kk.pdf")
end


N = 10
sobol = gsa(allprobs[1](d2), CVODE_BDF(),tim,[[exp(genotypes[1](d2)[i]-0.1),exp(genotypes[1](d2)[i]+0.1)] for i = 1:48],Sobol(N=N,order=[2]),1e-15,1e-15)

gr(ylabel="sensitivity",label="",xlabel="mode",title="")
scatter(exp.(d2),test, series_annotations = text.(labels[1:51], :bottom,6))
savefig("mode_vs_sensitivity_non_log2.pdf")