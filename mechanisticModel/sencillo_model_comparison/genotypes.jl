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
	pars[:k2]=pars[:mutk2]
	collect(pars)
end	
function snf3ko(params)
	pars=OrderedDict(params)
	pars[:dmth1snf3] = pars[:mutk3]
	pars[:k3] = log(10) #introduced in hxt4model6
	collect(pars)
end	

function simplesnf1(params)
	pars=OrderedDict(params)
	pars[ksnf2] = log(0.0)
	pars[] = log(10) #introduced in hxt4model6
	collect(pars)
end	



function composemutant(mut1, mut2)

	function composed(pars)
	
	mut1(mut2(pars))
	
	end
	
	end
	

function makehxt4(params, val)
	pars=OrderedDict(params)
	pars[:Hxt4_0] = dm2[val][1]
	
	collect(pars)
end	


templatenames=
[:k3,
:k2,
:ksnf1,
:ksnf1std1,
:nsnf1,
:nsnf2,
:Snf1tot,
:dmth1,
:nmth1snf3,
:nmth1rgt2,
:dmth1snf3,
:dmth1rgt2,
:smth1,
:kmig1mth1,
:nmig1mth1,
:kmig2mth1,
:nmig2mth1,
:std1tot,
:istd1,
:nstd1,
:estd1max,
:imig1,
:kmig1snf1,
:emig1max,
:dmig2,
:dmig2snf1,
:kmig2snf1,
:smig2,
:kmig2std1,
:nmig2std1,
:kmig2mth1std1,
:nmig2mth1std1,
:Mig1_0,
:Mig2_0,
:Mth1_0,
:Std1_0]
templatepars= [ j=>getvalue(bfpairs6, j) for j in templatenames] 





function maskregulators(params)
	template=templatepars
	pars=OrderedDict(params)
	
	[pars[x[1]] = x[2] for x in template]
	
	collect(pars)
end	






