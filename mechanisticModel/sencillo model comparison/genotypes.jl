#genotype functions alter parameters to simulate mutants.

function wt(params)
	params
end
 
function mth1ko(params)
 pars = copy(params)
	pars[13] = 0.0
	pars[47] = 0.0
	return collect(pars)
end
 
function std1ko(params)
 pars = copy(params)
	pars[18] = 0.0
	pars[48] = 0.0
	return collect(pars)
end	

function mig2ko(params)
 pars = copy(params)
	pars[29] = 0.0
 pars[46] = 0.0
 return collect(pars)
end	

function mig1ko(params)
 pars = copy(params)
	pars[22] = 0.0
 pars[45] = 0.0
 return collect(pars)
end	

function rgt2ko(params)
 pars = copy(params)
 pars[12] = 0.0
 pars[2] = 100.0 #should this be another parameter?
 return collect(pars)
end	

function snf3ko(params)
 pars = copy(params)
 pars[11] = 0.0
 return collect(pars)
end