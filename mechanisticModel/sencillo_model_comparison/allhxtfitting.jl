
g=1
pph=[]
allhxtpars=[]
hxtmeans=[]
hxtstds=[]
for h in ["hxt1", "hxt2", "hxt3", "hxt4", "hxt5", "hxt6", "hxt7"]

r="rep2"

hxtprocessed=[subnan.(j) for j in celldata[h]["g1percent"][r]]

hxtmean=[ mean(skipmissing([j[x] for j in hxtprocessed])) for x in 1:length(hxtprocessed[1])]
hxtstd=[ std(skipmissing([j[x] for j in hxtprocessed])) for x in 1:length(hxtprocessed[1])]
push!(hxtmeans, hxtmean)
push!(hxtstds, hxtstd)

#calling the genotype maskregulators, which only leaves hxt parameters (and the regulator affinities) as free parameters.
hxtprob=makeproblem(modelfile, tim[1:229], y, 1.0, tt, hxtmean[1], maskregulators) #concentration 1 percent

solvehxt=solveone(hxtprob, hxtmean, hxtstd) #function to make a function that depends on param array

#defining the blackbox optimisation
optshxt=bbsetup(solvehxt; Method = :adaptive_de_rand_1_bin_radiuslimited, SearchRange = [subss.((p[2]-5,p[2]+5)) for p in  logpars], NumDimensions = 52, MaxSteps = 100000)
#running the blackbox optimisation
reshxt=bboptimize(optshxt, MaxSteps=2000)

parshxt=best_candidate(reshxt); #parameter vector

push!(allhxtpars, parshxt)

#making paired array for the plotting
parpair=[bfpairs6[j][1] =>parshxt[j] for j in 1:length(bfpairs6)]

push!(pph, solveplotnew(hxtprob,maskregulators(parpair), hxtmean[1:length(x)],hxtstd[1:length(x)], 1))

end
plot(pph..., layout=7)



parpair=[bfpairs6[j][1] =>allhxtpars[4][j] for j in 1:length(bfpairs6)]
solveplotnew(hxtprob,maskregulators(parpair), hxtmean[1:length(x)]