#genotype functions alter parameters to simulate mutants.

function wt2(params)
    params
end

function mth1ko2(params)
    pars = copy(params)
    pars[13] = log(0.0)
    pars[43] = log(0.0)
    return collect(pars)
end

function std1ko2(params)
    pars = copy(params)
    pars[14] = log(0.0)
    pars[44] = log(0.0)
    return collect(pars)
end	

function mig2ko2(params)
    pars = copy(params)
    pars[25] = log(0.0)
    return collect(pars)
end	

function mig1ko2(params)
    pars = copy(params)
    pars[19] = log(0.0)
    return collect(pars)
end	

function rgt2ko2(params)
    pars = copy(params)
     pars[12] = log(0.0)
     mutk2 = pars[2] + pars[45]
     pars[2] = mutk2
    return collect(pars)
 end    

 function snf3ko2(params)
    pars = copy(params)
    pars[11] = log(0.0)
    mutk3 = pars[1] + pars[45]
    pars[1] = mutk3
    return collect(pars)
 end