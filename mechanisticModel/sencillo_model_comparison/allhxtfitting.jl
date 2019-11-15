pars5d= [ :k3=> 0.2,  :k2=> 0.6,  :ksnf1=> 0.4,  :ksnf1std1=> 0.1, 
 :nsnf1=> 2.3,  :nsnf2=> 2,  :Snf1tot=> 12,  :dmth1=> 0.3, 
  :nmth1snf3=> 1.2,  :nmth1rgt2=> 1.4,  :dmth1snf3=> 0.6,  
  :dmth1rgt2=> 0.7,  :smth1=> 0.3,  :kmig1mth1=> 0.2, 
   :nmig1mth1=> 1.3,  :kmig2mth1=> 12,  :nmig2mth1=> 1.3, 
    :std1tot=> 2.1,  :istd1=> 3.4,  :nstd1=> 1.5,  :estd1max=> 3.1, 
     :imig1=> 3.4,  :kmig1snf1=> 0.8,  :emig1max=> 3.1,  :dmig2=> 0.3, 
      :dmig2snf1=> 2.1,  :kmig2snf1=> 2.4,  :smig2=> 0.5,  :kmig2std1=> 2.1/2,  
      :nmig2std1=> 2.1,  :kmig2mth1std1=> 0.3,  :nmig2mth1std1=> 1.2,  :dhxt4=> 0.03,  
      :dhxt4max=> 2.1,  :kdhxt4=> 2.4,  :ndhxt4=> 1,  :shxt4=> 0.5,  :khxt4mth1=> 12,  
      :nhxt4mth1=> 1.3,  :khxt4std1=> 2.1/2,  :nhxt4std1=> 2.1,  :khxt4mth1std1=> 0.3,  
      :nhxt4mth1std1=> 1.2,  :khxt4mig1=> 2,  :khxt4mig2=> 1,  :nhxt4mig1=> 1.5,  
      :nhxt4mig2=> 2,  :Hxt4_0=> 0,  :Mig1_0=> 0,  :Mig2_0=> 0,  :Mth1_0=> 3/.3,  :Std1_0=> 2.1, 
]

subnan(x)=  if typeof(x)==String return missing else return x end 
optim=true
g=1
pph=[]
allhxtpars=[]
hxtmeans=[]
hxtstds=[]
allopts=[]
allres=[]
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
if optim==true
#defining the blackbox optimisation

optshxt=bbsetup(solvehxt; Method = :adaptive_de_rand_1_bin_radiuslimited, SearchRange = [subss.((log(p[2])-5,log(p[2])+5)) for p in  pars5d], NumDimensions = 52, MaxSteps = 100000)

#running the blackbox optimisation
reshxt=bboptimize(optshxt, MaxSteps=15000)


parshxt=best_candidate(reshxt); #parameter vector


#making paired array for the plotting
parpair=[pars5d[j][1] =>parshxt[j] for j in 1:length(pars5d)]


push!(allopts, optshxt)
push!(allres, reshxt)
push!(allhxtpars, parshxt)
push!(pph, solveplotnew(hxtprob,maskregulators(parpair), hxtmean[1:length(x)],hxtstd[1:length(x)], 1))
else

hnum=findall(hnames.==h)[1] #get the number of the hxt in turn

parpair=[bfpairs9[j][1] =>allhxtpars[hnum][j] for j in 1:length(bfpairs9)]
push!(pph, solveplotnew(hxtprob,maskregulators(parpair), hxtmeans[hnum][1:length(x)],hxtstds[hnum][1:length(x)], 1))

end 


end
plot(pph..., layout=7)




 hnames=["hxt1", "hxt2", "hxt3", "hxt4", "hxt5", "hxt6", "hxt7"]
 j=1
 pb=[]

bfpairs9=[:k3=>2.63508, :k2=>2.17225, :ksnf1=>-5.98817, :ksnf1std1=>9.05465, :nsnf1=>4.50887, :nsnf2=>-0.600449, :Snf1tot=>8.95374, :dmth1=>-0.815013, :nmth1snf3=>-0.221673, :nmth1rgt2=>-2.11001, :dmth1snf3=>3.05302, :dmth1rgt2=>-0.624956, :smth1=>-1.79843, :kmig1mth1=>0.743028, :nmig1mth1=>2.14294, :kmig2mth1=>0.535124, :nmig2mth1=>2.17213, :std1tot=>-3.11169, :istd1=>0.794552, :nstd1=>0.611865, :nstd3=>4.94466, :estd1max=>5.16118, :imig1=>6.1056, :kmig1snf1=>-1.66105, :emig1max=>0.105971, :dmig2=>1.26108, :dmig2snf1=>-3.94351, :kmig2snf1=>3.0439, :smig2=>-2.12053, :kmig2std1=>-3.38116, :nmig2std1=>0.45985, :kmig2mth1std1=>-5.18045, :nmig2mth1std1=>3.68821, :dhxt4=>-2.26598, :dhxt4max=>-1.42215, :kdhxt4=>2.2143, :ndhxt4=>-2.29222, :shxt4=>2.47129, :khxt4mth1=>-1.70927, :nhxt4mth1=>6.09682, :khxt4std1=>-3.14254, :nhxt4std1=>3.98171, :khxt4mth1std1=>-5.75592, :nhxt4mth1std1=>-0.601246, :khxt4mig1=>-3.43196, :khxt4mig2=>-4.51407, :nhxt4mig1=>0.412203, :nhxt4mig2=>1.52485, :Hxt4_0=>-104.157, :Mig1_0=>-97.565, :Mig2_0=>-102.962, :Mth1_0=>-1.47144, :Std1_0=>-2.7336, :mutk3=>1.83416, :mutk2=>4.14787] 
 

pb=[]

parpairs=[]

for ii in 34:48

push!(pb, bar([allhxtpars[1][ii],allhxtpars[2][ii],allhxtpars[3][ii],allhxtpars[4][ii],allhxtpars[6][ii],allhxtpars[7][ii]], title=bfpairs9[ii][1]))

end

plot(pb..., layout=(5,5), legend=false)
plot(pb...,legend=false)




parpair=[bfpairs6[j][1] =>allhxtpars[4][j] for j in 1:length(bfpairs6)]
solveplotnew(hxtprob,maskregulators(parpair), hxtmean[1:length(x)]











