#genotype functions alter parameters to simulate mutants.

function wt5(params)
    params
end

function mth1ko5(params)
    pars = copy(params)
    pars[11] = log(0.0)
    pars[41] = log(0.0)
    return collect(pars)
end

function std1ko5(params)
    pars = copy(params)
    pars[12] = log(0.0)
    pars[42] = log(0.0)
    return collect(pars)
end	

function mig2ko5(params)
    pars = copy(params)
    pars[23] = log(0.0)
    return collect(pars)
end	

function mig1ko5(params)
    pars = copy(params)
    pars[17] = log(0.0)
    return collect(pars)
end	

function rgt2ko5(params)
    pars = copy(params)
     pars[10] = log(0.0)
     mutk2 = pars[2] + pars[43]
     pars[2] = mutk2
    return collect(pars)
 end    

 function snf3ko5(params)
    pars = copy(params)
    pars[9] = log(0.0)
    mutk3 = pars[1] + pars[43]
    pars[1] = mutk3
    return collect(pars)
 end