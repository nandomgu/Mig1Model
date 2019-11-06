using Distributed #- can parallelise code here by adding extra cores...
addprocs(22)
@everywhere begin 
    using Polynomials 
    using StatsBase, KernelDensity, Random
    using Distances, Distributions
    import Distributions.rand
    import Distributions.pdf

    include("/home/msturrock/Desktop/Mig1/Mig1Model/mechanisticModel/sencillo_model_comparison/types.jl")
    include("/home/msturrock/Desktop/Mig1/Mig1Model/mechanisticModel/sencillo_model_comparison/functions.jl")
    include("/home/msturrock/Desktop/Mig1/Mig1Model/mechanisticModel/sencillo_model_comparison/abc2.jl")
    params = 5
    data = randn(params)
    times = range(0,stop=1,length=1000)
    poly_data = Poly(data,:x)
    poly_datas = poly_data(times)

    function rho_lens(expd,d2)
        model_data = Poly(d2,:y)
        model_datas = model_data(times)
        error = Distances.sqeuclidean(poly_datas,model_datas)
        return error
    end

end

A = Uniform(-5,5)
model_lens = repeat([A];outer=params+1)

np = 10000

@time apmc_output =APMC_KDE(np,0.0,[model_lens],[rho_lens],paccmin=0.01)

d2 = apmc_output.pts[end][1:(params+1),1]
d3 = apmc_output.pts[end][1:(params+1),1:1000]

using StatsPlots
pyplot()
StatsPlots.corrplot(d3')