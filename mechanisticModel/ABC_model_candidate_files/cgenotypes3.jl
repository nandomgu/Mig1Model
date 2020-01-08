#genotype functions alter parameters to simulate mutants.

function wt3(params)
    params
end

function mth1ko3(params)
    pars = copy(params)
    pars[11] = log(0.0)
    pars[45] = log(0.0)
    return collect(pars)
end

function std1ko3(params)
    pars = copy(params)
    pars[16] = log(0.0)
    pars[46] = log(0.0)
    return collect(pars)
end	

function mig2ko3(params)
    pars = copy(params)
    pars[27] = log(0.0)
    return collect(pars)
end	

function mig1ko3(params)
    pars = copy(params)
    pars[21] = log(0.0)
    return collect(pars)
end	

function rgt2ko3(params)
    pars = copy(params)
     pars[10] = log(0.0)
     mutk2 = pars[2] + pars[47]
     pars[2] = mutk2
    return collect(pars)
 end    

 function snf3ko3(params)
    pars = copy(params)
    pars[9] = log(0.0)
    mutk3 = pars[1] + pars[47]
    pars[1] = mutk3
    return collect(pars)
 end