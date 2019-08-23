#genotype functions alter parameters to simulate mutants.

function wt(params)
	params
end
	
function mth1ko(params)
	pars=OrderedDict(params)
	pars[:smth1]=0.0
	pars[:Mth1_0]=0.0
	collect(pars)
end
	
function std1ko(params)
	pars=OrderedDict(params)
	pars[:std1tot]=0.0
	pars[:Std1_0]=0.0
	collect(pars)
end	
function mig2ko(params)
	pars=OrderedDict(params)
	pars[:smig2]=0.0
	pars[:Mig2_0]=0.0
	collect(pars)
end	
function mig1ko(params)
	pars=OrderedDict(params)
	pars[:mig1tot]=0.0
	pars[:Mig1_0]=0.0
	collect(pars)
end	
function rgt2ko(params)
	pars=OrderedDict(params)
	pars[:dmth1rgt2] = 0.0
	pars[:dstd1rgt2] = 0.0
	collect(pars)
end	
function snf3ko(params)
	pars=OrderedDict(params)
	pars[:dmth1snf3] = 0.0
	pars[:dstdtsnf3] = 0.0
	collect(pars)
end	

