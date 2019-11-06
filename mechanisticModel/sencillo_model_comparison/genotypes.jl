#genotype functions alter parameters to simulate mutants.

function wt(params)
    params
end

function mth1ko(params)
    pars = copy(params)
    pars[13] = log(0.0)
    pars[49] = log(0.0)
    return collect(pars)
end

function std1ko(params)
    pars = copy(params)
    pars[18] = log(0.0)
    pars[50] = log(0.0)
    return collect(pars)
end	

function mig2ko(params)
    pars = copy(params)
    pars[29] = log(0.0)
    return collect(pars)
end	

function mig1ko(params)
    pars = copy(params)
    pars[23] = log(0.0)
    return collect(pars)
end	

function rgt2ko(params)
    pars = copy(params)
     pars[12] = log(0.0)
     mutk2 = pars[2] + pars[51]
     pars[2] = mutk2
    return collect(pars)
 end    

 function snf3ko(params)
    pars = copy(params)
    pars[11] = log(0.0)
    mutk3 = pars[1] + pars[51]
    pars[1] = mutk3
    return collect(pars)
 end