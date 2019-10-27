pars5d= [
:k3=> 0.2, 
:k2=> 0.6, 
:ksnf1=> 0.4, 
:ksnf1std1=> 0.1, 
:nsnf1=> 2.3, 
:nsnf2=> 2, 
:Snf1tot=> 12, 
:dmth1=> 0.3, 
:nmth1snf3=> 1.2, 
:nmth1rgt2=> 1.4, 
:dmth1snf3=> 0.6, 
:dmth1rgt2=> 0.7, 
:smth1=> 0.3, 
:kmig1mth1=> 0.2, 
:nmig1mth1=> 1.3, 
:kmig2mth1=> 12, 
:nmig2mth1=> 1.3, 
:std1tot=> 2.1, 
:istd1=> 3.4, 
:nstd1=> 1.5, 
:estd1max=> 3.1, 
:imig1=> 3.4, 
:kmig1snf1=> 0.8, 
:emig1max=> 3.1, 
:dmig2=> 0.3, 
:dmig2snf1=> 2.1, 
:kmig2snf1=> 2.4, 
:smig2=> 0.5, 
:kmig2std1=> 2.1/2, 
:nmig2std1=> 2.1, 
:kmig2mth1std1=> 0.3, 
:nmig2mth1std1=> 1.2, 
:dhxt4=> 0.03, 
:dhxt4max=> 2.1, 
:kdhxt4=> 2.4, 
:ndhxt4=> 1, 
:shxt4=> 0.5, 
:khxt4mth1=> 12, 
:nhxt4mth1=> 1.3, 
:khxt4std1=> 2.1/2, 
:nhxt4std1=> 2.1, 
:khxt4mth1std1=> 0.3, 
:nhxt4mth1std1=> 1.2, 
:khxt4mig1=> 2, 
:khxt4mig2=> 1, 
:nhxt4mig1=> 1.5, 
:nhxt4mig2=> 2, 
:Hxt4_0=> 0, 
:Mig1_0=> 0, 
:Mig2_0=> 0, 
:Mth1_0=> 3/.3, 
:Std1_0=> 2.1, 
]

mypars= [
:k3=> 0.1, 
:k2=> 0.6, 
:ksnf1=> 0.2, 
:ksnf1std1=> 160.0, 
:nsnf1=> 8.0, 
:nsnf2=> 8.0, 
:Snf1tot=> 1000, 
:dmth1=> 0.3, 
:nmth1snf3=> 1.2, 
:nmth1rgt2=> 1.4, 
:dmth1snf3=> 4, 
:dmth1rgt2=> 2, 
:smth1=> 1, 
:kmig1mth1=> 0.2, 
:nmig1mth1=> 1.3, 
:kmig2mth1=> 12, 
:nmig2mth1=> 1.3, 
:std1tot=> 2.1, 
:istd1=> 4, 
:nstd1=> 3, 
:estd1max=> 100, 
:imig1=>100, 
:kmig1snf1=> 3, 
:emig1max=>15, 
:dmig2=> 0.3, 
:dmig2snf1=> 2.1, 
:kmig2snf1=> 2.4, 
:smig2=> 10, 
:kmig2std1=> 10, 
:nmig2std1=> 2.1, 
:kmig2mth1std1=> 0.3, 
:nmig2mth1std1=> 1.2, 
:dhxt4=> 0.03, 
:dhxt4max=> 2.1, 
:kdhxt4=> 2.4, 
:ndhxt4=> 1, 
:shxt4=> 10, 
:khxt4mth1=> .5, 
:nhxt4mth1=> 5, 
:khxt4std1=>5 , 
:nhxt4std1=> 2.1, 
:khxt4mth1std1=> 20, 
:nhxt4mth1std1=> 1.2, 
:khxt4mig1=> 10, 
:khxt4mig2=> 10, 
:nhxt4mig1=> 1.5, 
:nhxt4mig2=> 2, 
:Hxt4_0=> 0, 
:Mig1_0=> 0, 
:Mig2_0=> 0, 
:Mth1_0=> 0.3/.3, 
:Std1_0=>2.1, 
]

#function fillparswith()
completepars=[];
for j in 1:length(pars5d)
try
push!(completepars, pars5d[j][1]=> getvalue(mypars,  pars5d[j][1])) 
catch
push!(completepars, pars5d[j][1]=> getvalue(pars5d,  pars5d[j][1])) 
end
end

conc=3
logpars=[p[1]=>log(p[2]) for p in completepars]
pp2=solveplotnew(allprobs[conc], logpars, dm2[conc][1:length(x)],me2[conc][1:length(x)], 2)
pp3=solveplotnew(allprobs[conc], logpars, dm2[conc][1:length(x)],me2[conc][1:length(x)], 3)
pp4=solveplotnew(allprobs[conc], logpars, dm2[conc][1:length(x)],me2[conc][1:length(x)], 4)
pp5=solveplotnew(allprobs[conc], logpars, dm2[conc][1:length(x)],me2[conc][1:length(x)], 5)
pp6=solveplotnew(allprobs[conc], logpars, dm2[conc][1:length(x)],me2[conc][1:length(x)], 6)
pp1=solveplotnew(allprobs[conc], logpars, dm2[conc][1:length(x)],me2[conc][1:length(x)], 1)
plot(pp2,pp3, pp4,pp5, pp6,pp1, layout=(3, 2), ylim=[0,5])











