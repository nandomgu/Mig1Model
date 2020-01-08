#genotype functions alter parameters to simulate mutants.

function wt6(params)
    params
end

function mth1ko6(params)
    pars = copy(params)
    pars[11] = log(0.0)
    pars[39] = log(0.0)
    return collect(pars)
end

function std1ko6(params)
    pars = copy(params)
    pars[12] = log(0.0)
    pars[40] = log(0.0)
    return collect(pars)
end	

function mig2ko6(params)
    pars = copy(params)
    pars[23] = log(0.0)
    return collect(pars)
end	

function mig1ko6(params)
    pars = copy(params)
    pars[17] = log(0.0)
    return collect(pars)
end	

function rgt2ko6(params)
    pars = copy(params)
     pars[10] = log(0.0)
     mutk2 = pars[2] + pars[41]
     pars[2] = mutk2
    return collect(pars)
 end    

 function snf3ko6(params)
    pars = copy(params)
    pars[9] = log(0.0)
    mutk3 = pars[1] + pars[41]
    pars[1] = mutk3
    return collect(pars)
 end