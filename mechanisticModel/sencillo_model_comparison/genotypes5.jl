#genotype functions alter parameters to simulate mutants.

function wt(params)
	params
end
	
function mth1ko(params)
	pars=OrderedDict(params)
	pars[:smth1]=log(0.0)
	pars[:Mth1_0]=log(0.0)
	collect(pars)
end
	
function std1ko(params)
	pars=OrderedDict(params)
	pars[:std1tot]=log(0.0)
	pars[:Std1_0]=log(0.0)
	collect(pars)
end	
function mig2ko(params)
	pars=OrderedDict(params)
	pars[:smig2]=log(0.0)
	pars[:Mig2_0]=log(0.0)
	collect(pars)
end	
function mig1ko(params)
	pars=OrderedDict(params)
	pars[:mig1tot]=log(0.0)
	pars[:imig1]=log(0.0)
	pars[:Mig1_0]=log(0.0)
	collect(pars)
end	
function rgt2ko(params)
	pars=OrderedDict(params)
	pars[:dmth1rgt2] = log(0.0)
	if :dstd1rgt2 in OrderedDict(params).keys
	pars[:dstd1rgt2] =log(0.0)
	end
	pars[:k2]=log(10)
	collect(pars)
end	
function snf3ko(params)
	pars=OrderedDict(params)
	pars[:dmth1snf3] = log(0.0)
	pars[:k3] = log(10) #introduced in hxt4model6
	collect(pars)
end	


function makehxt4(params, val)
	pars=OrderedDict(params)
	pars[:Hxt4_0] = dm2[val][1]
	
	collect(pars)
end	



function maskregulators(params)
	template=bfpairs6[1:32]
	pars=OrderedDict(params)
	
	[pars[x[1]] = x[2] for x in template]
	
	collect(pars)
end	






