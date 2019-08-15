using DifferentialEquations, Sundials, Dierckx, JSON, JLD2, BlackBoxOptim, StatsPlots, MATLAB


#importing data through json
#these are cell traces
celldata=JSON.parsefile("/Users/s1259407/Dropbox/PhD/phd_peter_swain/R/hxtmeandata.json", dicttype=Dict)
#=these are the means for the 18 conditions.
 each 3 is one genotype, and for each genotype
  there are three glucose concentrations (0.2%, 0.4%, 1%)
=#
datameans=JSON.parsefile("/Users/s1259407/Documents/MATLABGIT/fitdatameans.json", dicttype=Dict)
#ordering the data so that each entry is one phenotype
dm2=[[datameans[j][k] for j in 1:size(datameans)[1]] for k in 1:size(datameans[1])[1]]

#glucose data and creating its interpolant
glucose=[8.32733959421089e-19,2.561481685921e-19,1.04197119608466e-19,5.53738039368354e-20,2.30389033029131e-20,1.01730601065369e-20,9.64607628015416e-21,4.54189858719156e-21,2.62590228463374e-21,1.07030601256657e-21,7.12893122193276e-22,3.77665654101745e-22,1.13613450886179e-22,9.11107194068206e-23,1.68345898188741e-23,4.59571704388656e-23,2.4737145794224e-23,5.06371909794061e-24,5.65818382020064e-24,1.04048279473923e-24,3.94925490074619e-24,3.45909892619726e-25,7.07762515569138e-25,7.59557441484042e-26,2.02035130943482e-27,3.27268072550391e-26,1.56537448533959e-27,1.6320588929363e-27,4.39103892339638e-28,4.80386219709977e-29,3.78228642078199e-28,6.7008774844888e-33,3.38770680654291e-30,3.11390728411053e-36,3.43524261389185e-36,9.0953685838273e-37,0,2.65814991285393e-16,1.40534785861741e-05,0.956449369761021,0.999932720540147,0.999988138202835,0.999994336997356,0.999995688515055,0.99999655799382,0.999997230933311,0.999997349602989,0.999997796320276,0.999997725777247,0.999997825419588,0.999997959294029,0.99999807290326,0.999998148682981,0.999998097810728,0.999998115954266,0.999998085169392,0.999998126522822,0.99999814575262,0.999998172022714,0.999998338656721,0.999998383673738,0.999998504114664,0.999998515535666,0.999998447595621,0.999998528661984,0.99999855188616,0.999998519653104,0.999998565797203,0.999998565177499,0.999998598958849,0.999998669110336,0.999998719774422,0.999998776796943,0.999998794501031,0.99999883534008,0.999998798817409,0.99999880270023,0.999998808299886,0.999998809973147,0.999998846846304,0.999998882337051,0.999998957349787,0.999998958954537,0.999998984982437,0.999998977466028,0.999998955080815,0.999998965848435,0.999998965848435,0.999998955080815,0.999998977466028,0.999998984982437,0.999998958954537,0.999998957349787,0.999998882337051,0.999998846846304,0.999998809973147,0.999998808299886,0.99999880270023,0.999998798817409,0.99999883534008,0.999998794501031,0.999998776796943,0.999998719774422,0.999998669110336,0.999998598958849,0.999998565177499,0.999998565797203,0.999998519653104,0.99999855188616,0.999998528661984,0.999998447595621,0.999998515535666,0.999998504114664,0.999998383673738,0.999998338656721,0.999998172022714,0.99999814575262,0.999998126522822,0.999998085169392,0.999998115954266,0.999998097810728,0.999998148682981,0.99999807290326,0.999997959294029,0.999997825419588,0.999997725777247,0.999997796320276,0.999997349602989,0.999997230933311,0.99999655799382,0.999995688515055,0.999994336997356,0.999988138202835,0.999932720540147,0.956449369761021,1.40534785861741e-05,2.65814991285393e-16,0,9.0953685838273e-37,3.43524261389185e-36,3.11390728411053e-36,3.38770680654291e-30,6.7008774844888e-33,3.78228642078199e-28,4.80386219709977e-29,4.39103892339638e-28,1.6320588929363e-27,1.56537448533959e-27,3.27268072550391e-26,2.02035130943482e-27,7.59557441484042e-26,7.07762515569138e-25,3.45909892619726e-25,3.94925490074619e-24,1.04048279473923e-24,5.65818382020064e-24,5.06371909794061e-24,2.4737145794224e-23,4.59571704388656e-23,1.68345898188741e-23,9.11107194068206e-23,1.13613450886179e-22,3.77665654101745e-22,7.12893122193276e-22,1.07030601256657e-21,2.62590228463374e-21,4.54189858719156e-21,9.64607628015416e-21,1.01730601065369e-20,2.30389033029131e-20,5.53738039368354e-20,1.04197119608466e-19,2.561481685921e-19,8.32733959421089e-19,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
#times at which glucose is sampled
tim=[0 0.0833333333333333 0.166666666666667 0.25 0.333333333333333 0.416666666666667 0.5 0.583333333333333 0.666666666666667 0.75 0.833333333333333 0.916666666666667 1 1.08333333333333 1.16666666666667 1.25 1.33333333333333 1.41666666666667 1.5 1.58333333333333 1.66666666666667 1.75 1.83333333333333 1.91666666666667 2 2.08333333333333 2.16666666666667 2.25 2.33333333333333 2.41666666666667 2.5 2.58333333333333 2.66666666666667 2.75 2.83333333333333 2.91666666666667 3 3.08333333333333 3.16666666666667 3.25 3.33333333333333 3.41666666666667 3.5 3.58333333333333 3.66666666666667 3.75 3.83333333333333 3.91666666666667 4 4.08333333333333 4.16666666666667 4.25 4.33333333333333 4.41666666666667 4.5 4.58333333333333 4.66666666666667 4.75 4.83333333333333 4.91666666666667 5 5.08333333333333 5.16666666666667 5.25 5.33333333333333 5.41666666666667 5.5 5.58333333333333 5.66666666666667 5.75 5.83333333333333 5.91666666666667 6 6.08333333333333 6.16666666666667 6.25 6.33333333333333 6.41666666666667 6.5 6.58333333333333 6.66666666666667 6.75 6.83333333333333 6.91666666666667 7 7.08333333333333 7.16666666666667 7.25 7.33333333333333 7.41666666666667 7.5 7.58333333333333 7.66666666666667 7.75 7.83333333333333 7.91666666666667 8 8.08333333333333 8.16666666666667 8.25 8.33333333333333 8.41666666666667 8.5 8.58333333333333 8.66666666666667 8.75 8.83333333333333 8.916666666666679 9.08333333333333 9.16666666666667 9.25 9.33333333333333 9.41666666666667 9.5 9.58333333333333 9.66666666666667 9.75 9.83333333333333 9.91666666666667 10 10.0833333333333 10.1666666666667 10.25 10.3333333333333 10.4166666666667 10.5 10.5833333333333 10.6666666666667 10.75 10.8333333333333 10.9166666666667 11 11.0833333333333 11.1666666666667 11.25 11.3333333333333 11.4166666666667 11.5 11.5833333333333 11.6666666666667 11.75 11.8333333333333 11.9166666666667 12 12.0833333333333 12.1666666666667 12.25 12.3333333333333 12.4166666666667 12.5 12.5833333333333 12.6666666666667 12.75 12.8333333333333 12.9166666666667 13 13.0833333333333 13.1666666666667 13.25 13.3333333333333 13.4166666666667 13.5 13.5833333333333 13.6666666666667 13.75 13.8333333333333 13.9166666666667 14 14.0833333333333 14.1666666666667 14.25 14.3333333333333 14.4166666666667 14.5 14.5833333333333 14.6666666666667 14.75 14.8333333333333 14.9166666666667 15 15.0833333333333 15.1666666666667 15.25 15.3333333333333 15.4166666666667 15.5 15.5833333333333 15.6666666666667 15.75 15.8333333333333 15.9166666666667 16 16.0833333333333 16.1666666666667 16.25 16.3333333333333 16.4166666666667 16.5 16.5833333333333 16.6666666666667 16.75 16.8333333333333 16.9166666666667 17 17.0833333333333 17.1666666666667 17.25 17.3333333333333 17.4166666666667 17.5 17.5833333333333 17.6666666666667 17.75 17.8333333333333 17.9166666666667 18 18.0833333333333 18.1666666666667 18.25 18.3333333333333 18.4166666666667 18.5 18.5833333333333 18.6666666666667 18.75 18.8333333333333 18.9166666666667 19 19.0833333333333]

#making a general glucose spline for every purpose
gspline = Spline1D(tim[1,:], glucose[1:229]; k=1)

#function to substitute nans
subs(x)=  if typeof(x)==String return missing else return x end

#UNUSED: hxt4 mean data splines
hspline=Spline1D(tim[1:229], dm2[3][1:229]) #hxt4 mean in 1% glucose
hspline2=Spline1D(tim[1:229], dm2[1][1:229]) #hxt4 mean in 0.2% glucose
hspline4=Spline1D(tim[1:229], dm2[2][1:229]) #hxt4 mean in 0.4% glucose

#=for each data mean, make a
 spline interpolant that can be sampled at any arbitrary time.
and can be called in a loop.
=#
splines=Array{Spline1D}(undef, (18)) #making splines for original data
for j in 1:18
if j==7
splines[j]=Spline1D(tim[ findall(.!ismissing.([subs(x) for x in dm2[j][1:229]]))], dm2[j][ findall(.!ismissing.([subs(x) for x in dm2[j][1:229]]))])
else
splines[j]=Spline1D(tim[1:229], dm2[j][1:229])
end
end

##plotting the interpolations to see if they worked.
interps=[[splines[j](x) for x in 0.:.1:20.0] for j in 1:18]
plot(0.:.1:20.0, interps)


function ablate(theta, ind)
	#=for theta being a binary array where each entry implies whether a specific
	allele/ species/connection is present,
	ablate transforms entry ind () into a zero to simulate a mutation or different genotype.
=#
abarray=ones(size(theta))
abarray[ind]=0
return theta.*abarray
end


#=DECLARING THE GENOTYPE ARRAY THETA
we fix the entries in an arbitrary order
values are the default values but this way we can wuickly redefine the array.
=#
thetamig1=1
thetamth1=1
thetastd1=1
thetasnf3=1
thetargt2=1
thetamig2=1
thetamig1mth1=0
thetamig2mth1=1
thetasnf1mig1=1
thetasnf1mig2=1
thetastd1hxt4=1
thetastd1mig2=1
thetamth1mig2=1

theta=[thetamig1, thetamth1, thetastd1, thetasnf3, thetargt2, thetamig2, thetamig1mth1, thetamig2mth1,thetasnf1mig1,thetasnf1mig2, thetastd1hxt4, thetastd1mig2, thetamth1mig2]


#= generating thetas for each of the conditions we have data for,
and according to theorder of the condition means set we imported.
=#
function makethetamatrix(theta)
allthetas=[
theta,#wt
theta,
theta,

ablate(theta,1), #mig1ko
ablate(theta,1),
ablate(theta,1),

ablate(theta,2),#mth1ko
ablate(theta,2),
ablate(theta,2),

ablate(theta,3),#std1ko
ablate(theta,3),
ablate(theta,3),

ablate(theta,5),#rgt2ko
ablate(theta,5),
ablate(theta,5),

ablate(theta,4),#snf3ko
ablate(theta,4),
ablate(theta,4),
]

end


#bundling parameters into groups.
mig1params=[20 21 22 23]
hxt4repression=[35:44]
hxt4synthesis =[31,32,33,34]


#=FUNCTION TO MAKE AN ODE PROBLEM
AS A F(GLUCOSE_TRACE, G_LEVEL, GENOTYPE, INITIAL_CONDITIONS, TSPAN)
=#
function makeproblem(input, concentration, theta, inits, t)

	#IMPORTING THE GENOTYPES
	thetamig1=theta[1]
	thetamth1=theta[2]
	thetastd1=theta[3]
	thetasnf3=theta[4]
	thetargt2=theta[5]
	thetamig2=theta[6]
	thetamig1mth1=theta[7]
	thetamig2mth1=theta[8]
	thetasnf1mig1=theta[9]
	thetasnf1mig2=theta[10]
	thetastd1hxt4=theta[11]
	thetastd1mig2=theta[12]
	thetamth1mig2=theta[13]

#MODIFYING THE INITIAL CONDITIONS IN CASE THE GENOTYPE INCLUDES A MUTANT
inits2=copy(inits)
 if thetamig1==0
	inits2[2]=inits[2]*thetamig1
end
if thetamig2==0
	inits2[3]=inits[3]*thetamig2
end
if thetamth1==0
	inits2[4]=inits[4]*thetamth1
end
if thetastd1==0
	inits2[5]=inits[5]*thetastd1
end


#SYSTEM OF ODES WITH A TIME VARYING INPUT
	function hxtsim(dydt, y, parameters, t)
		g=input(t)*concentration
		dhxt4= parameters[1]
		degmth1= parameters[2]
		dmig2= parameters[3]
		ksnf3= parameters[4]
		krgt2= parameters[5]
		ksnf1= parameters[6]
		ksnf1std1= parameters[7]
		nsnf1= parameters[8]
		nmth1snf3= parameters[9]
		nmth1rgt2= parameters[10]
		smth1= parameters[11]
		kmig1mth1= parameters[12]
		nmig1mth1= parameters[13]
		kmig2mth1= parameters[14]
		nmig2mth1= parameters[15]
		std1tot= parameters[16]
		istd1= parameters[17]
		nstd1= parameters[18]
		estd1max= parameters[19]
		mig1tot= parameters[20]
		imig1= parameters[21]
		kmig1snf1= parameters[22]
		emig1max= parameters[23]
		dmig2snf1= parameters[24]
		kmig2snf1= parameters[25]
		smig2= parameters[26]
		kmig2std1= parameters[27]
		nmig2std1= parameters[28]
		kmig2mth1std1= parameters[29]
		nmig2mth1std1= parameters[30]
		dhxt4max= parameters[31]
		kdhxt4= parameters[32]
		ndhxt4= parameters[33]
		shxt4= parameters[34]
		khxt4mth1= parameters[35]
		nhxt4mth1= parameters[36]
		khxt4std1= parameters[37]
		nhxt4std1= parameters[38]
		khxt4mth1std1= parameters[39]
		nhxt4mth1std1= parameters[40]
		khxt4mig1= parameters[41]
		khxt4mig2= parameters[42]
		nhxt4mig1= parameters[43]
		nhxt4mig2= parameters[44]

		Hxt4= y[1]
		Mig1= y[2]
		Mig2= y[3]
		Mth1= y[4]
		Std1= y[5]

		Snf1= (ksnf1 + ksnf1std1*Std1^nsnf1)/((ksnf1 + ksnf1std1*Std1)^nsnf1 + g^nsnf1)
		dydt[1]= +(shxt4/(1 + thetamig1* (Mig1/khxt4mig1)^nhxt4mig1 +  thetamig2* (Mig2/khxt4mig2)^nhxt4mig2 + thetamth1*(Mth1/khxt4mth1)^nhxt4mth1 + thetastd1hxt4*thetastd1*(Std1/khxt4std1)^nhxt4std1 + thetamth1*thetastd1*(Mth1*Std1/ khxt4mth1std1)^nhxt4mth1std1))-dhxt4*Hxt4-(dhxt4max/(1 + (g/kdhxt4)^ndhxt4)*Hxt4)
		dydt[2]= thetamig1*(imig1*(mig1tot - Mig1))-(emig1max*Snf1*Mig1/(kmig1snf1 + Mig1))
		dydt[3]= +(smig2/(1 + thetamth1mig2*thetamth1*(Mth1/kmig2mth1)^nmig2mth1 + thetastd1mig2*thetastd1*(Std1/kmig2std1)^nmig2std1 + thetamth1mig2*thetastd1mig2*thetastd1*thetamth1*(Mth1*Std1/ kmig2mth1std1)^nmig2mth1std1))-dmig2*Mig2-(dmig2snf1*Snf1*Mig2/(kmig2snf1 + Mig2))
		#dydt[4]= +(smth1/(1 + (Mig1/kmig1mth1)^nmig1mth1 + (Mig2/kmig2mth1)^nmig2mth1))-degmth1*Mth1-(g^nmth1snf3/(ksnf3^nmth1snf3 + g^nmth1snf3)*Mth1)-(g^nmth1rgt2/(krgt2^nmth1rgt2 + g^nmth1rgt2)*Mth1)
		dydt[4]= thetamth1* (smth1/(1 + thetamig1mth1*(Mig1/kmig1mth1)^nmig1mth1 + thetamig2mth1*(Mig2/kmig2mth1)^nmig2mth1)-degmth1*Mth1- (thetasnf3*(g^nmth1snf3/(ksnf3^nmth1snf3 + (g^nmth1snf3)))*Mth1)- thetargt2*(g^nmth1rgt2/(krgt2^nmth1rgt2 + g^nmth1rgt2))*Mth1)
		dydt[5]= thetastd1*(istd1*(std1tot - Std1)- thetargt2*(estd1max*g^nstd1/(krgt2^nstd1 + g^nstd1)*Std1)) #+ thetasnf3*(estd1max*g^nstd1/(ksnf3^nstd1 + g^nstd1)*Std1))

	end
#OUTPUTING A FUNCTION THAT GENERATES A PROBLEM AS A FUNCTION OF ONLY PARAMETERS.
prob(params)= ODEProblem(hxtsim, inits2, t, exp.(params))

end



####### TUNING PARAMETERS MANUALLY


# define rate constants
dhxt4= 0.03
degmth1= 0.3
dmig2= 1

# parameters
ksnf3= 0.2
krgt2= 0.6

ksnf1= 0.4
ksnf1std1= 0.1
nsnf1= 2.3

nmth1snf3= 1.2
nmth1rgt2= 1.4
smth1= 0.3
kmig1mth1= 0.2
nmig1mth1= 1.3
kmig2mth1= 2
nmig2mth1= 1.3

std1tot= 1
istd1= .5
nstd1= 1.5
estd1max= 2

mig1tot= 1
imig1= 3
kmig1snf1= 2
emig1max= 12

dmig2snf1= 10
kmig2snf1= 5
smig2= 0.5
kmig2std1= std1tot/2
nmig2std1= 2.1
kmig2mth1std1= 0.3
nmig2mth1std1= 1.2
dhxt4max= 2.1
kdhxt4= 2.4
ndhxt4= 1
shxt4= 0.5
khxt4mth1= .6
nhxt4mth1= 2
khxt4std1= std1tot/2
nhxt4std1= 2.1
khxt4mth1std1= 0.3
nhxt4mth1std1= 1.2
khxt4mig1= .3
khxt4mig2= .1
nhxt4mig1= 2
nhxt4mig2= 1.5
#              1      2        3      4     5       6      7         8      9          10         11       12         13        14          15        16       17    18      19        20      21      22         23
parameters= [dhxt4, degmth1, dmig2, ksnf3, krgt2, ksnf1, ksnf1std1, nsnf1, nmth1snf3, nmth1rgt2, smth1, kmig1mth1, nmig1mth1, kmig2mth1, nmig2mth1, std1tot, istd1, nstd1, estd1max, mig1tot, imig1, kmig1snf1, emig1max, dmig2snf1, kmig2snf1, smig2, kmig2std1, nmig2std1, kmig2mth1std1, nmig2mth1std1, dhxt4max, kdhxt4, ndhxt4, shxt4, khxt4mth1, nhxt4mth1, khxt4std1, nhxt4std1, khxt4mth1std1, nhxt4mth1std1, khxt4mig1, khxt4mig2, nhxt4mig1, nhxt4mig2]
parnames= ["dhxt4", "degmth1", "dmig2", "ksnf3", "krgt2", "ksnf1", "ksnf1std1", "nsnf1", "nmth1snf3", "nmth1rgt2", "smth1", "kmig1mth1", "nmig1mth1", "kmig2mth1", "nmig2mth1", "std1tot", "istd1", "nstd1", "estd1max", "mig1tot", "imig1" , "kmig1snf1", "emig1max", "dmig2snf1", "kmig2snf1", "smig2", "kmig2std1", "nmig2std1", "kmig2mth1std1", "nmig2mth1std1" , "dhxt4max", "kdhxt4", "ndhxt4", "shxt4", "khxt4mth1", "nhxt4mth1", "khxt4std1", "nhxt4std1", "khxt4mth1std1", "nhxt4mth1std1","khxt4mig1","khxt4mig2", "nhxt4mig1", "nhxt4mig2"]

#CREATING A PARAMETER SEARCH RANGE FOR PARAMETERS BASED ON PRIOR KNOWLEDGE
lb=0.00001 #instead of zero
sb=50  # max synthesis rates
sk=20
sn=6
	srange= [
		(lb,.3), #dhxt4,
	 	(lb,.3),# degmth1,
	 	(lb, .3),#dmig2,
		(lb, 1.0),#ksnf3,
		(lb,1.0),#krgt2,
		(lb,1.0),#ksnf1,
		(lb,1.0),#ksnf1std1,
		(1.0,6.0),#nsnf1,
		(1.0,6.0),#nmth1snf3,
		(1.0,6.0),#nmth1rgt2,
		(lb,1),#smth1,
		(lb, sk),#kmig1mth1,
		(1.0,6.0),#nmig1mth1,
		(lb,1.0),#kmig2mth1,
		(1.0,6.0),#nmig2mth1,
		(0.99,1.0),#std1tot,
		(lb, sb),#istd1,
		(1.0,6.0),# nstd1,
		(lb,sb),# estd1max,
		(0.99,1.0),# mig1tot,
		(lb,sb),# imig1,
		(lb,1.0),# kmig1snf1,
		(lb, sb),# emig1max,
		(lb,sb),# dmig2snf1,
		(lb, sk),# kmig2snf1,
		(lb,  sb),# smig2,
		(lb, sk),# kmig2std1,
		(1.0,6.0),# nmig2std1,
		(lb,sk),# kmig2mth1std1,
		(1.0,6.0),# nmig2mth1std1,
		(lb,1.),# dhxt4max,
		(lb,1.),# kdhxt4,
		(1.0,6.0),# ndhxt4,
		(lb,sb),# shxt4,
		(lb, sk),# khxt4mth1,
		(1.0,6.0),#  nhxt4mth1,
		(lb,sk),#  khxt4std1,
		(1.0,6.0),#  nhxt4std1,
		(lb,sk),#  khxt4mth1std1,
		(1.0,6.0),#   nhxt4mth1std1,
		(lb,sk),#    khxt4mig1,
		(lb,sk),#	khxt4mig2,
		(1.0,6.0),#	nhxt4mig1,
		(1.0,6.0)#	nhxt4mig2,
	]
	#SUBSTITUTING INFINITE VALUES FOR A VERY SMALL NUMBER INSTEAD.
subss(x)=  if !isfinite(x) return -30.0 else return x end
lrange=[log.(j) for j in  srange] # getting some infs because of zeros
trange=[subss.(j) for j in lrange]


# define initial conditions
Hxt4_0= 0
Mig1_0= 0
Mig2_0= 0
Mth1_0= smth1/degmth1
Std1_0= std1tot
inits= [ Hxt4_0, Mig1_0, Mig2_0, Mth1_0, Std1_0]

#TIME SPAN INPUT
t= (0., 20.)
species= [:Hxt4, :Mig1, :Mig2, :Mth1, :Std1]



function squares1(params)
               prob= ODEProblem(hxtsim_1, inits, t, exp.(params))
               try
               sol= solve(prob,TRBDF2(), saveat=.5)
               arr=[[j[i] for j in sol.u] for i=1:length(sol.u[1])]
               return sum((arr[1]-hspline(sol.t)).^2)
       catch
               return 1000000.0
       end
       end


function squares_serial(params)
	#solving the systems for 3 different glucose levels one after the other
	prob= ODEProblem(hxtsim_1, inits, t, exp.(params))
	try

	sol= solve(prob,TRBDF2(), saveat=.5)
	arr=[[j[i] for j in sol.u] for i=1:length(sol.u[1])]
	lsq1=sum((arr[1]-hspline(sol.t)).^2)
	catch
		return 1000000.0
	end
	try
	prob= ODEProblem(hxtsim_2, inits, t, exp.(params))
	sol= solve(prob,TRBDF2(), saveat=.5)
	arr=[[j[i] for j in sol.u] for i=1:length(sol.u[1])]
	lsq2=sum((arr[1]-hspline(sol.t)).^2)
catch
	lsq2=1000000.0
end
try
	prob= ODEProblem(hxtsim_4, inits, t, exp.(params))
	sol= solve(prob,TRBDF2(), saveat=.5)
	arr=[[j[i] for j in sol.u] for i=1:length(sol.u[1])]
	lsq4=sum((arr[1]-hspline(sol.t)).^2)
catch
lsq4=1000000.0
end

	return lsq1+lsq2+lsq4

end



#declaring theta to activate or inactivate components in model

###
#For default tolerances, AutoTsit5(Rosenbrock23())
#is a good choice. For lower tolerances,
#using AutoVern7 or AutoVern9 with Rodas4, KenCarp4,
#or Rodas5 can all be good choices depending on the problem.
#For very large systems (>1000 ODEs?), consider using lsoda.
###

#=
For stiff problems at high tolerances
(>1e-2?) it is recommended that you use Rosenbrock23
or TRBDF2. These are robust to oscillations and massive stiffness
is needed, though are only efficient when low accuracy is needed.
Rosenbrock23 is more efficient for small systems where re-evaluating and
 re-factorizing the Jacobian is not too costly,
 and for sufficiently large systems TRBDF2 will be more efficient.
  ABDF2 can be the most efficient the largest systems
  or most expensive f.

At medium tolerances (>1e-8?) it is recommended you use Rodas5,
 Rodas4P (the former is more efficient but
  the later is more reliable), Kvaerno5, or KenCarp4.
   As native DifferentialEquations.jl solvers, many
   Julia numeric types (such as BigFloats, ArbFloats, or DecFP)
   will work. When the equation is defined via the @ode_def macro,
   these will be the most efficient.

For faster solving at low tolerances (<1e-9)
 but when Vector{Float64} is used, use radau.

For asymptotically large systems of ODEs (N>1000?) where f is very costly and the complex eigenvalues are minimal (low oscillations), in that case CVODE_BDF will be the most efficient but requires Vector{Float64}. CVODE_BDF will also do surprisingly well if the solution is smooth. However, this method can be less stiff than other methods and stuff may fail at low accuracy situations. Another good choice for this regime is lsoda.
=#

#parameters must be in log
#serial simulation ready for optimisation
solver=ABDF2
function serialall(pars)
#plot(1)
finalarrs=[]
lsq=zeros(18)
concs=repeat([0.2, 0.4, 1], 6)

for k in 1:size(concs,1)

prob=makeproblem(gspline,concs[k], allthetas[k], inits, t)
try
sol= solve(prob(pars),ABDF2(), saveat=0.1)
arr=[[j[i] for j in sol.u] for i=1:length(sol.u[1])]
#push!(finalarrs, arr[1])
#plot!(arr[1], xlims=(0,20))
lsq[k]=sum((arr[1]-splines[k](sol.t)).^2)
catch
	try
		sol= solve(prob(pars),TRBDF2(), saveat=0.1)
		arr=[[j[i] for j in sol.u] for i=1:length(sol.u[1])]
		push!(finalarrs, arr[1])
		lsq[k]=sum((arr[1]-splines[k](sol.t)).^2)
	catch
	lsq[k]=1000000.0
end
end
end


return sum(lsq)
end

#simulate and calculate the square difference for each condition
#using 2 simulators, and pick the one with the best score
function makesimulator( theta, gspline, inits, t, outvar=1)
allthetas=makethetamatrix(theta)
function serialall3(pars, allthetas)
#plot(1)
#=
type 1 outputs the lsq sum
type 2 outputs the simulated data
type 3 outputs the lsq for each condition
=#
type=2
finalarr= a=fill([], 18)
finalt=fill([], 18)
lsq=zeros(18)
concs=repeat([0.2, 0.4, 1], 6)
xx=[]
yy=[]
for k in 1:size(concs,1)

prob=makeproblem(gspline,concs[k], allthetas[k], inits, t)
try
sol= solve(prob(pars),ABDF2(), saveat=0.1)
arr=[[j[i] for j in sol.u] for i=1:length(sol.u[1])]
finalarr[k]=arr[outvar]
finalt[k]=sol.t
#push!(finalarrs, arr[1])
#plot!(arr[1], xlims=(0,20))
xx=sum((arr[outvar]-splines[k](sol.t)).^2)
catch
	xx=1000000.0
end
	try
		sol= solve(prob(pars),TRBDF2(),  saveat=0.1)
		arr=[[j[i] for j in sol.u] for i=1:length(sol.u[1])]
		#push!(finalarrs, arr[1])
		finalarr[k]=arr[outvar]
		finalt[k]=sol.t
		yy=sum((arr[outvar]-splines[k](sol.t)).^2)
	catch
	yy=1000000.0
end

lsq[k]= min(xx,yy)
end

if type==1
	return sum(lsq)
elseif type==2
	return finalt, finalarr
elseif type==3
	return lsq
end


end

func(pars)= serialall3(pars, allthetas)
end

#troubleshooting simulation.
function serialallts(pars)
#plot(1)
finalarrs=[]
lsq=zeros(18)
concs=repeat([0.2, 0.4, 1], 6)
xx=[]
yy=[]
zz=[]
for k in 1:size(concs,1)
println("condition ", k)
println("solver ABDF2")
prob=makeproblem(gspline,concs[k], allthetas[k], inits, t)
try
sol= solve(prob(pars),ABDF2(), saveat=0.1)
arr=[[j[i] for j in sol.u] for i=1:length(sol.u[1])]
#push!(finalarrs, arr[1])
#plot!(arr[1], xlims=(0,20))
xx=sum((arr[1]-splines[k](sol.t)).^2)
println("result: ", xx)
catch err
	println("simulation failed")
	println(err)
	xx=1000000.0
end
try
println("solver TRBDF2")
		sol= solve(prob(pars),TRBDF2(), saveat=0.1)
		arr=[[j[i] for j in sol.u] for i=1:length(sol.u[1])]
		#push!(finalarrs, arr[1])
		yy=sum((arr[1]-splines[k](sol.t)).^2)
println("result: ", yy)
catch err
	println("simulation failed")
	println(err)
	yy=1000000.0
end

try
println("solver Rosenbrock23")
		sol= solve(prob(pars),Rosenbrock23(), saveat=0.1)
		arr=[[j[i] for j in sol.u] for i=1:length(sol.u[1])]
		#push!(finalarrs, arr[1])
		zz=sum((arr[1]-splines[k](sol.t)).^2)
println("result: ", yy)
catch err
	println("simulation failed")
	println(err)
	zz=1000000.0
end



lsq[k]= min(xx,yy,zz)
end

return sum(lsq), lsq
end




###simulate and plot  all simulation results
function serialallplot(pars)
#plot(1)
finalarrs=[]
lsq=zeros(18)
concs=repeat([0.2, 0.4, 1], 6)
xx=[]
yy=[]
for k in 1:size(concs,1)
arr=[]
prob=makeproblem(gspline,concs[k], allthetas[k], inits, t)
println(k, concs[k])
try
println("trying ABDF2\n")
sol= solve(prob(pars),ABDF2(), saveat=0.1)
global tt=sol.t #making sure we define a time that can be used in the next for loop
arr=[[j[i] for j in sol.u] for i=1:length(sol.u[1])]

#plot!(arr[1], xlims=(0,20))

xx=sum((arr[1]-splines[k](sol.t)).^2)
catch err
	print("simulation failed")
	xx=1000000.0
end
#println(xx)
	try
		println("trying TRBDF2\n")
		sol= solve(prob(pars),TRBDF2(), saveat=0.1)
		arr=[[j[i] for j in sol.u] for i=1:length(sol.u[1])]
		yy=sum((arr[1]-splines[k](sol.t)).^2)
	catch err
		print("simulation failed")
	yy=1000000.0
end
println("printing")
if isempty(arr)
	println(finalarrs[k-1])
push!(finalarrs, zeros(size(finalarrs[1])))
else
push!(finalarrs, arr[1])
end
#println(yy)
lsq[k]= min(xx,yy)

end
strainnames=["WT","mig1ko","mth1ko","std1ko","rgt2ko","snf3ko"]
plots=[];
c=1
      for j=1:3:16
       push!(plots, plot( tt, finalarrs[j:j+2], color=[:cyan :blue :purple],title=strainnames[c],legend=false, titlefontsize=5))
c+=1
       end
plot(plots[1],plots[2], plots[3], plots[4], plots[5], plots[6], layout=(3,2))
end




#serial simulation to test solvers
solvers=[Rosenbrock23, TRBDF2,ABDF2, Rodas4P ,Kvaerno5 ,KenCarp4]
function serialall2(pars, solver)
#plot(1)
finalarrs=[]
lsq=zeros(18)
concs=repeat([0.2, 0.4, 1], 6)
cols=repeat(["cyan", "blue", "purple"], outer=6)
for k in 1:size(concs,1)
try
prob=makeproblem(gspline,concs[k], allthetas[k], inits, t)
sol= solve(prob(pars),solver(), saveat=0.1)
arr=[[j[i] for j in sol.u] for i=1:length(sol.u[1])]
push!(finalarrs, arr[1])
#plot!(arr[1], xlims=(0,20))
lsq[k]=sum((arr[1]-splines[k](sol.t)).^2)
catch
	lsq[k]=1000000.0
end

end

return sum(lsq), lsq
end
#=
lss=[]
for j in 1:size(solvers)[1]
a,b= serialall(theta2, solvers[j])
push!(lss,b)
end
=#




solver=ABDF2
function simplot(params, prob)
	readyprob=prob(params)
	try
	sol= solve(readyprob,TRBDF2(), saveat=0.1)
	arr=[[j[i] for j in sol.u] for i=1:length(sol.u[1])]
	plot(sol.t, arr, xlims=(0,20), label=["hxt4", "mig1", "mig2", "mth1", "std1"])
catch
	sol= solve(readyprob,ABDF2(), saveat=0.1)
	arr=[[j[i] for j in sol.u] for i=1:length(sol.u[1])]
	plot(sol.t, arr, xlims=(0,20), label=["hxt4", "mig1", "mig2", "mth1", "std1"])


end

end


#thetax= optimizing only hxt4 in 1%
#thetay=optimizing for 3 concentrations of Glucose
#thetaz=optimizing to means of 18 conditions. mig1 and mth1 not active in theta

#nice parameter sets.

#thetax contains the result tor a model fit towards 1% glucose, yielding a great mig1 signal
thetax=[-1.64592, -11.0823, -6.7671, -11.3897, -11.5007, -9.22642, -0.0383132, 0.0130669, 0.134756, 0.356679, -5.90065, -8.20684, 1.22262, -7.71933, 0.72516, -0.00869444, -6.52181, 0.000410178, 0.945155, -0.00343608, 0.600393, -9.05274, 0.564541, 0.362751, -0.131216, 2.16333, -2.99267, 1.34913, -1.85758, 1.08103, -2.5563, -2.30758, 0.0265852, 2.18615, 0.819032, 1.3156, 0.112461, 1.18497, -11.4782, 1.38626, -1.26902, 2.25975, 1.38305, 1.21805]


#improving parameters step by step
#adding the mig1params
niceparams=[]
push!(niceparams, [0.03, 0.3, 1.0, 0.2, 0.6, 0.4, 0.1, 2.3, 1.2, 1.4, 0.3, 0.2, 1.3, 2.0, 1.3, 1.0, 0.5, 1.5, 2.0, 3.0, 1.82284, 1.0e-5, 5.0, 10.0, 5.0, 0.5, 0.5, 2.1, 0.3, 4.0, 2.1, 2.4, 1.0, 20.0, 0.25, 3.0, 0.5, 2.1, 0.3, 1.2, 0.3, 0.1, 2.0, 1.5])
#adjusted mig2 and its regulation by std1
push!(niceparams, [0.03, 0.3, 1.0, 0.2, 0.6, 0.4, 0.1, 2.3, 1.2, 1.4, 0.3, 0.2, 1.3, 2.0, 1.3, 1.0, 0.5, 1.5, 2.0, 3.0, 1.82284, 1.0e-5, 5.0, 10.0, 5.0, 0.5, 0.7, 2.1, 0.3, 4.0, 2.1, 2.4, 1.0, 20.0, 0.5, 3.0, 0.5, 2.1, 0.3, 1.2, 0.3, 0.08, 2.0, 1.5])
#after adjusting mig2 levels we could bring up the synthesis rate of hxt4 by quite a lot
push!(niceparams, [0.03, 0.3, 1.0, 0.2, 0.6, 0.4, 0.1, 2.3, 1.2, 1.4, 0.3, 0.2, 1.3, 2.0, 1.3, 1.0, 0.5, 1.5, 2.0, 3.0, 1.82284, 1.0e-5, 5.0, 10.0, 5.0, 0.5, 0.7, 2.1, 0.3, 4.0, 2.1, 2.4, 1.0, 80.0, 0.5, 3.0, 0.5, 2.1, 0.3, 1.2, 0.3, 0.08, 2.0, 1.5])
#probably best fit of complete model (without mig1 and mig2)#theta3
push!(niceparams,[-1.20397, -1.22114, -10.1282, -8.80975, -9.42375, -7.75983, -1.40067, 0.00419071, 1.72936, 1.68641, -2.25353, -5.14219, 0.111147, -0.121496, 1.70058, -0.00637065, -1.01043, 0.694172, -0.366817, -0.00421826, 3.42508, -8.44178, 3.25086, -9.15642, 0.936593, -9.82395, -10.1233, 1.73947, -6.62262, 0.52591, -3.02366, -10.8725, 0.30445, 1.6401, -3.14036, 0.463717, -0.00867208, 1.44995, 2.48249, 0.306599, -0.661812, 1.6845, 1.76194, 1.62666])
# theta4. initial values reach zero, apaprently through abuse of std1
push!(niceparams, [-1.22299, -1.22403, -11.4048, -8.74496, -9.00203, -8.35082, -0.872288, 0.398006, 1.70334, 1.26115, -2.36072, -9.21626, 1.0542, -6.33612, 1.01746, -0.00328004, -0.772175, 1.15681, 0.46746, -0.00709036, 3.41669, -4.5556, 3.25259, -4.8109, 2.14409, -11.4927, -8.06542, 0.723956, -7.44838, 1.50541, -3.06343, -9.94579, 1.18679, 1.59893, -3.18582, 0.382377, 0.0152555, 1.66978, -2.23489, 1.41341, -0.61033, -4.78659, 1.76044, 1.68482])

### thetamm, theta 4 with a corrected mig2 dynamics which is actually responsive to glucose

push!(niceparams,[-1.22299, -1.22403, 0.0, -8.74496, -9.00203, -8.35082, -0.872288, 0.398006, 1.70334, 1.26115, -2.36072, -9.21626, 1.0542, 0.693147, 0.262364, -0.00328004, -0.772175, 1.15681, 0.46746, -0.00709036, 3.41669, -4.5556, 3.25259, 2.30259, 1.60944, -0.693147, -0.356675, 0.741937, -1.20397, 1.38629, -3.06343, -9.94579, 1.18679, 2.5, -3.18582, 0.382377, 0.0152555, 1.66978, -2.23489, 1.41341, -0.61033, -1.38, 1.76044, 1.68482])


#= theta7. outstanding score. of 2424 or so. skills mig2 dynamics as this undermines expression levels
therefore different levels of the mig1ko are not well fit.
importantly, the std1 ko
 would bring mig2 levels to
 the roof which would punish hxt4 dynamics.
naturally this means that  the snf3 and rgt2
 were not fit very well either since std1 levels
 are probably not responding much to glucose
=#
push!(niceparams,[-1.22299, -1.22403, -11.4048, -8.74496, -9.00203, -8.35082, -0.872288, 0.398006, 1.70334, 1.26115, -2.36072, -9.21626, 1.0542, -6.33612, 1.01746, -0.00328004, -0.772175, 1.15681, 0.46746, -0.00709036, 3.41669, -4.5556, 3.25259, -4.8109, 2.14409, -11.4927, -8.06542, 0.723956, -7.44838, 1.50541, -3.06343, -9.94579, 1.18679, 1.59893, -3.18582, 0.382377, 0.0152555, 1.66978, -2.23489, 1.41341, -0.61033, -4.78659, 1.76044, 1.68482])
#= theta8. very promising set.
regulation of mth1 over mig2 was reconstituted,
which was the reason why mig2 was going crazy in a std1 deletion.
after reconstitution the fit makes almost complete sense.
sensor mutants have overall lower levels than the wild type, but
the qualitative phenomenon is still not greatly fit
=#
push!(niceparams,[-0.999852, -1.26372, 0.403835, -9.24132, -9.42964, -8.19311, -0.486726, 0.176302, 1.44706, 1.47452, -2.55226, -9.45124, 1.42551, 0.411843, 0.288631, -0.220279, -0.551536, 1.16724, -0.0210462, 0.481953, 3.33342, -4.49532, 3.38703, 2.12734, 1.27722, -0.59074, 0.139484, 1.15084, -0.783974, 0.880301, -2.7309, -10.146, 0.47693, 2.02273, -2.97313, 0.850846, -0.213677, 1.3482, -1.76923, 0.829166, -0.222927, -1.41487, 1.95319, 2.08884])
#theta 9. maxing out the score from previous round. with thetamig2mth2 on. 2395
push!(niceparams,[-0.999852, -1.26372, 0.403832, -9.24132, -9.42964, -8.19326, -0.486726, 0.176302, 1.44706, 1.47452, -2.55226, -9.06174, 0.259888, 0.411843, 0.288631, -0.220279, -0.551536, 1.16724, -0.0210462, 0.489727, 3.32591, -4.50128, 3.38703, 2.12734, 1.27722, -0.59074, 0.139484, 1.15084, -0.783974, 0.880301, -2.7309, -10.146, 0.476931, 2.02273, -2.97313, 0.850846, -0.213677, 1.34808, -1.76923, 0.829166, -0.215866, -1.41487, 1.95319, 2.08884])

#=theta 10. this contains a functional version of the rgt2 mutant (thetargt2) over std1,
which was not included before. this yields a qualitatively accurate model.
this version still does not capture the effect of the snf3 mutant,
showing low glucose expression when it should be gravely diminished.
=#
push!(niceparams,[-1.09883, -1.2407, 0.597548, -8.85889, -9.09216, -8.00027, -0.797185, 0.0149666, 1.10511, 1.84656, -2.50323, -9.14348, 0.537932, 0.346894, 0.321236, 0.0224742, -0.818667, 1.22954, -0.390288, 0.234259, 3.18123, -4.12202, 2.99951, 2.1878, 1.59932, -0.418731, 0.611756, 1.10609, -0.563993, 1.36311, -2.74879, -10.0942, 0.350835, 2.10463, -3.04331, 0.893875, -0.0754482, 1.55275, -2.00151, 1.015, -0.590292, -1.36041, 1.8196, 2.49419])

#theta11. reoptimisation and retrieval of full float numbers since some of these are quite crucial for the very best fit.
push!(niceparams, [-1.1546942232439559, -1.0887322970066144, 0.7132935756630442, -8.383416706182768, -8.808551521233642, -7.663202682401169, -0.8490911301724443, 0.01942088874544573, 0.8095384545117216, 1.1452570081581448, -2.542702044988643, -9.170542009379027, 0.9808344028633106, 0.28656559774427465, 0.37786454104517064, 0.10068945545460253, -0.9264555562961837, 0.8472390538177788, -0.35235592543987754, 0.22184004676371108, 3.119263346803683, -3.958848702341918, 2.800400741543633, 1.957424585540502, 1.5110268619211527, -0.10241264506686007, 0.5794946667353741, 1.1430127832452404, -0.7549280339651918, 1.7370792587195765, -2.9277235184473414, -10.063799157410548, 0.2902557443376703, 2.137503781773038, -3.045651546769331, 1.0621630786570095, -0.07740907614272813, 1.3161948485904997, -2.369571961625992, 0.5449457722738957, -0.46257286319590596, -1.108626011434643, 2.0978946817087203, 2.5243779161783877])



#same as above but we give it an accessible name. parameters are in log.
pars=[-1.1546942232439559, -1.0887322970066144, 0.7132935756630442, -8.383416706182768, -8.808551521233642, -7.663202682401169, -0.8490911301724443, 0.01942088874544573, 0.8095384545117216, 1.1452570081581448, -2.542702044988643, -9.170542009379027, 0.9808344028633106, 0.28656559774427465, 0.37786454104517064, 0.10068945545460253, -0.9264555562961837, 0.8472390538177788, -0.35235592543987754, 0.22184004676371108, 3.119263346803683, -3.958848702341918, 2.800400741543633, 1.957424585540502, 1.5110268619211527, -0.10241264506686007, 0.5794946667353741, 1.1430127832452404, -0.7549280339651918, 1.7370792587195765, -2.9277235184473414, -10.063799157410548, 0.2902557443376703, 2.137503781773038, -3.045651546769331, 1.0621630786570095, -0.07740907614272813, 1.3161948485904997, -2.369571961625992, 0.5449457722738957, -0.46257286319590596, -1.108626011434643, 2.0978946817087203, 2.5243779161783877]


##creating ranges centered on a parameter set
#find which parameter is not a hill factor, otherwise lower bound is 1 (zero in log)
isnothill=convert.(Float64,  [j[1]!='n' for j in  parnames])
#pars=copy(theta3)
#we assume parameters are in log
nt=[((pars[j]-.5)*isnothill[j],pars[j]+.5) for j in 1:size(pars)[1]]


###making prior distributions array based on the best parameters found
#pars=copy(theta11)
#degradation mapped to normal distribution
dss=[]
for j in 1:size(pars)[1]
if parnames[j][1]=='d' || parnames[j][1]=='k' || parnames[j][1]=='s' || parnames[j][1]=='e' || parnames[j][1]=='i'  && !occursin("tot",  parnames[1])
	push!(dss, Normal(pars[j], 1))
elseif parnames[j][1]=='n'
	push!(dss,   Uniform(0,4))
else
	push!(dss, Uniform(0.99, 1.001))

end
end


#########TEST CODE
#########test run of an individual simulation
#best fit parameters
pars=[-1.1546942232439559, -1.0887322970066144, 0.7132935756630442, -8.383416706182768, -8.808551521233642, -7.663202682401169, -0.8490911301724443, 0.01942088874544573, 0.8095384545117216, 1.1452570081581448, -2.542702044988643, -9.170542009379027, 0.9808344028633106, 0.28656559774427465, 0.37786454104517064, 0.10068945545460253, -0.9264555562961837, 0.8472390538177788, -0.35235592543987754, 0.22184004676371108, 3.119263346803683, -3.958848702341918, 2.800400741543633, 1.957424585540502, 1.5110268619211527, -0.10241264506686007, 0.5794946667353741, 1.1430127832452404, -0.7549280339651918, 1.7370792587195765, -2.9277235184473414, -10.063799157410548, 0.2902557443376703, 2.137503781773038, -3.045651546769331, 1.0621630786570095, -0.07740907614272813, 1.3161948485904997, -2.369571961625992, 0.5449457722738957, -0.46257286319590596, -1.108626011434643, 2.0978946817087203, 2.5243779161783877]
#genotype. we start with the wild type.

#first six entries define the genes that are present or absent.
#following entries define presence of regulatory branches.
#e.g. thetamig1mth1=0 means that mig1 is set to not not repress mth1.
#full theta is the full model, with all branches present (except for the mig1mth1 branch which just doesn't fit the data)
fulltheta=[thetamig1, thetamth1, thetastd1, thetasnf3, thetargt2, thetamig2, thetamig1mth1, 1,thetasnf1mig1,thetasnf1mig2, thetastd1hxt4, thetastd1mig2, thetamth1mig2]
#this line sets up an array of 18 genotypes so that all mutants can be simulated for one parameter set.
#allthetas is called within serialall3

simulate=makesimulator(fulltheta, gspline, inits, t)
lsqsumfull= simulate(pars)
#plotting the result for all simulations
#serialallplot(pars)
#generating all model variants in the wt version.
#the above theta from 1 to 7, followed by an array between 000000 and 111111
allvariants=[vcat(theta[1:7], k[:])  for k in [convert.(Int, mat"de2bi($j,6)") for j in 0:63]]

lsqsummut=ones(size(allvariants)[1])*NaN
global c=1
#for this parameter set, simulate all model variants
for j in 1:size(allvariants)[1]
global c
#simulating some variant of the model
#making a simulator function for this model variant
simulate=makesimulator(allvariants[j], gspline, inits, t)
lsqsummut[c]= simulate(pars)
c=c+1
end


###make sure type is set to 2 in makesimulator
#make sure to set the right outvar in  makesimulator
for outvar in 1:5 #we loop through all species
sim=makesimulator(theta, gspline, inits, t, outvar)
a,b=sim(pars) #a is the times and b is the simulated datapoints
pl=[]
speciesnames= ["Hxt4", "Mig1", "Mig2", "Mth1", "Std1"]
strainnames=["WT", "mig1ko", "mth1ko", "std1ko", "rgt2ko", "snf3ko"]
strainnamesall=repeat(strainnames, inner=3)

[push!(pl, plot(a[k:k+2], b[k:k+2], color= [:cyan :blue :purple], w=2,legend=false, title=speciesnames[outvar]*" in "*strainnamesall[k]      ) ) for k in 1:3:16]
#adding scatter plots of real data
plot(pl[1], pl[2], pl[3], pl[4], pl[5], pl[6])
if outvar==1
[plot!(pl[k], tim[1:229], [splines[j](tim[1:229]) for j in ((k*3)-2):(k*3)], color= [:cyan :blue :purple], alpha=.5, m=([:cyan :blue :purple], 2, Plots.stroke(0, :cyan)), legend=false, title=strainnames[k]) for k in 1:6]
end
plot(pl[1], pl[2], pl[3], pl[4], pl[5], pl[6])
savefig(speciesnames[outvar]*"_dynamics.pdf")
end
