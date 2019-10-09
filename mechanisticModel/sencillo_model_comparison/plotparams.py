###PYTHON
bfpairs6=   {'k3': -1.99995, 'k2': 1.06211, 'ksnf1': -5.13512, 'ksnf1std1': 4.67567, 'nsnf1': 5.01806, 'nsnf2': -1.84446, 'Snf1tot': 5.73796, 'dmth1': -2.87079, 'nmth1snf3': 0.0698773, 'nmth1rgt2': 0.288116, 'dmth1snf3': 0.376694, 'dmth1rgt2': 3.14651, 'smth1': -4.50517, 'kmig1mth1': -0.379322, 'nmig1mth1': 5.25098, 'kmig2mth1': 0.560657, 'nmig2mth1': -1.74271, 'std1tot': -0.245076, 'istd1': 1.22922, 'nstd1': 0.746876, 'estd1max': 4.85039, 'imig1': 5.45447, 'kmig1snf1': -2.80427, 'emig1max': -0.264936, 'dmig2': 3.03049, 'dmig2snf1': 5.25449, 'kmig2snf1': 5.86184, 'smig2': 2.28341, 'kmig2std1': 5.81895, 'nmig2std1': 1.63477, 'kmig2mth1std1': 2.87004, 'nmig2mth1std1': 2.80064, 'dhxt4': -1.15641, 'dhxt4max': -2.48345, 'kdhxt4': -3.90747, 'ndhxt4': -4.58572, 'shxt4': 1.46517, 'khxt4mth1': -4.65211, 'nhxt4mth1': 3.03817, 'khxt4std1': 0.611965, 'nhxt4std1': 1.35178, 'khxt4mth1std1': 6.74873, 'nhxt4mth1std1': -0.211267, 'khxt4mig1': -1.23877, 'khxt4mig2': -2.57697, 'nhxt4mig1': 1.36639, 'nhxt4mig2': 2.07541, 'Hxt4_0': -100.0, 'Mig1_0': -100.0, 'Mig2_0': -100.0, 'Mth1_0': -4.54752, 'Std1_0': 5.57239, 'mutk2': 2.3025850929940459, 'mutk3': 2.3025850929940459}
marcpars_sd={'k3': -1.60724, 'k2': -1.11798, 'ksnf1': -7.70941, 'ksnf1std1': 9.91471, 'nsnf1': 7.07934, 'nsnf2': -0.0734542, 'Snf1tot': 8.70169, 'dmth1': -1.59598, 'nmth1snf3': 5.18232, 'nmth1rgt2': 0.551847, 'dmth1snf3': 3.9686, 'dmth1rgt2': 1.6702, 'smth1': -4.23113, 'kmig1mth1': 0.680474, 'nmig1mth1': 4.64808, 'kmig2mth1': 7.48488, 'nmig2mth1': 3.65958, 'std1tot': 4.35926, 'istd1': -3.53894, 'nstd1': 0.412305, 'estd1max': -0.707678, 'imig1': 6.1692, 'kmig1snf1': 3.85898, 'emig1max': 6.94579, 'dmig2': 0.993898, 'dmig2snf1': -4.25699, 'kmig2snf1': -4.12453, 'smig2': 4.08724, 'kmig2std1': 1.81226, 'nmig2std1': 2.52904, 'kmig2mth1std1': -3.75305, 'nmig2mth1std1': -4.81768, 'dhxt4': -2.88543, 'dhxt4max': -1.25432, 'kdhxt4': 4.27927, 'ndhxt4': 4.98911, 'shxt4': 6.1547, 'khxt4mth1': -5.69309, 'nhxt4mth1': 1.51019, 'khxt4std1': -0.620259, 'nhxt4std1': 0.491589, 'khxt4mth1std1': -2.32399, 'nhxt4mth1std1': 2.23031, 'khxt4mig1': -4.30685, 'khxt4mig2': -2.41038, 'nhxt4mig1': 0.63114, 'nhxt4mig2': 0.368315,  'Hxt4_0': -100.0, 'Mig1_0': 7.84135, 'Mig2_0': 9.99994, 'Mth1_0': 8.05267, 'Std1_0': -4.72919, 'mutk2': 0.517147, 'mutk3': -2.11891}

pars=bfpairs6






hillrep=lambda maxv, repressor, k, n: maxv/ (1+ (repressor**n/k**n))

hilldivisor=lambda repressor, k, n: (1+ (repressor**n/k**n))

hillrep2=lambda maxv, parameters1, parameters2: maxv/ (hilldivisor(*parameters1)+ hilldivisor(*parameters2))

fun= lambda g, Snf1tot, Std1, ksnf1std1, nsnf1, ksnf1, nsnf2: Snf1tot*(1 + Std1/ksnf1std1)**nsnf1/((1 + Std1/ksnf1std1)**nsnf1 + (g/ksnf1)**nsnf2)   



plt.figure()
c=1
nrow=4
numplots=3


### mth1 and std1 repression over hxt4
plt.subplot(nrow, numplots,c)
mth1min=0
mth1max=.015 
std1min=0
std1max=6
resp=np.zeros([tot, tot])*np.nan
mth1=np.linspace(mth1max, mth1min, tot)
for j in range(0, len(mth1)):
	resp[j, :]=[hillrep2(exp(pars["shxt4"]), [std1, exp(pars["khxt4std1"]),exp(pars["nhxt4std1"])], [mth1[j],exp(pars["khxt4mth1"]),exp(pars["nhxt4mth1"])]) for std1 in  np.linspace(std1min, std1max, tot)]
plt.imshow(resp, cmap='gnuplot2', extent=[std1min, std1max, mth1min, mth1max], aspect='auto'); plt.show()   
title('hxt4 exp. rate') 
plt.xlabel('Std1')
plt.ylabel('mth1') 


### mth1 and std1 dimer repression over hxt4

c+=1;
plt.subplot(nrow, numplots,c)
mth1min=0
mth1max=.015 
std1min=0
std1max=6
resp=np.zeros([tot, tot])*np.nan
mth1=np.linspace(mth1max, mth1min, tot)
for j in range(0, len(mth1)):
	resp[j, :]=[hillrep(exp(pars["shxt4"]), std1*mth1[j], exp(pars["khxt4mth1std1"]),exp(pars["nhxt4mth1std1"])) for std1 in  np.linspace(std1min, std1max, tot)]
plt.imshow(resp, cmap='gnuplot2', extent=[std1min, std1max, mth1min, mth1max], aspect='auto'); plt.show()   
title('hxt4 exp. rate \n(dimer)') 
plt.xlabel('Std1')
plt.ylabel('mth1') 


#mig1 mig2 repression over hxt4
c+=1;
plt.subplot(nrow, numplots,c)
mig1min=0
mig1max=1
mig2min=0
mig2max=.2
resp=np.zeros([tot, tot])*np.nan
mig2=np.linspace(mig2max, mig2min, tot)
for j in range(0, len(mth1)):
	resp[j, :]=[hillrep2(exp(pars["shxt4"]), [mig1, exp(pars["khxt4mig1"]),exp(pars["nhxt4mig1"])], [mig2[j],exp(pars["khxt4mig2"]),exp(pars["nhxt4mig2"])]) for mig1 in  np.linspace(mig1min, mig1max, tot)]
plt.imshow(resp, cmap='gnuplot2', extent=[mig1min, mig1max, mig2min, mig2max], aspect='auto'); plt.show()   
title('Hxt4 exp.rate') 
plt.xlabel('Mig1')
plt.ylabel('Mig2') 


#### mig1 and mig2 over mth1
c+=1;
plt.subplot(nrow, numplots,c)
mig1min=0
mig1max=1
mig2min=0
mig2max=1000
resp=np.zeros([tot, tot])*np.nan
mig2=np.linspace(mig2max, mig2min, tot)
for j in range(0, len(mth1)):
	resp[j, :]=[hillrep2(exp(pars["smth1"]), [mig1, exp(pars["kmig1mth1"]),exp(pars["nmig1mth1"])], [mig2[j],exp(pars["kmig2mth1"]),exp(pars["nmig2mth1"])]) for mig1 in  np.linspace(mig1min, mig1max, tot)]
plt.imshow(resp, cmap='gnuplot2', extent=[mig1min, mig1max, mig2min, mig2max], aspect='auto'); plt.show()   
title('Mth1 exp.rate') 
plt.xlabel('Mig1')
plt.ylabel('Mig2') 


hillinput=lambda vmax, input, k, n: vmax* (input**n)/(k**n+ input**n)  
#glucose over  mth1 and std1

#mth1 and std1 over mig2
c+=1;
plt.subplot(nrow, numplots,c)
mth1min=0
mth1max=10
std1min=0
std1max=400
resp=np.zeros([tot, tot])*np.nan
mth1=np.linspace(mth1max, mth1min, tot)
for j in range(0, len(mth1)):
	resp[j, :]=[hillrep2(exp(pars["smig2"]), [std1, exp(pars["kmig2std1"]),exp(pars["nmig2std1"])], [mth1[j],exp(pars["kmig2mth1"]),exp(pars["nmig2mth1"])]) for std1 in  np.linspace(std1min, std1max, tot)]
plt.imshow(resp, cmap='gnuplot2', extent=[std1min, std1max, mth1min, mth1max], aspect='auto'); plt.show()   
title('Mig2 exp. rate') 
plt.xlabel('Std1')
plt.ylabel('mth1') 


### mth1 and std1 dimer repression over mig2

c+=1;
plt.subplot(nrow, numplots,c)
mth1min=0
mth1max=10
std1min=0
std1max=10
resp=np.zeros([tot, tot])*np.nan
mth1=np.linspace(mth1max, mth1min, tot)
for j in range(0, len(mth1)):
	resp[j, :]=[hillrep(exp(pars["smig2"]), std1*mth1[j], exp(pars["kmig2mth1std1"]),exp(pars["nmig2mth1std1"])) for std1 in  np.linspace(std1min, std1max, tot)]
plt.imshow(resp, cmap='gnuplot2', extent=[std1min, std1max, mth1min, mth1max], aspect='auto'); plt.show()   
title('mig2 exp. rate\n(dimer)') 
plt.xlabel('Std1')
plt.ylabel('mth1') 


##snf1 activity
c+=1;
plt.subplot(nrow, numplots,c)
cmap='gnuplot2'
tot=500
glucmin=.00000001
glucmax=.000001
std1min=0
std1max=1
resp=np.zeros([tot, tot])*np.nan
glus=np.linspace(glucmax, glucmin, tot)
for j in range(0, len(glus)):
	resp[j, :]=[fun(glus[j], exp(pars['Snf1tot']), std1, exp(pars['ksnf1std1']), exp(pars['nsnf1']) , exp(pars['ksnf1']), exp(pars['nsnf2'])) for std1 in np.linspace(std1min, std1max, tot)]
plt.imshow(resp, cmap=cmap, extent=[std1min, std1max, glucmin, glucmax], aspect='auto' ); plt.show()   
title('Snf1 activity') 
plt.xlabel('Std1')
plt.ylabel('Glucose') 

######mig1 nuclear exclusion as dependent on  std1 and glucose
c+=1;
plt.subplot(nrow, numplots,c)
mth1min=0
mth1max=50
std1min=0
std1max=1.5
glucmin=0
glucmax=1
resp=np.zeros([tot, tot])*np.nan
mth1=np.linspace(mth1max, mth1min, tot)
glus=np.linspace(glucmax, glucmin, tot)
for j in range(0, len(glus)):
	#getting snf1
	s=[fun(glus[j], exp(pars['Snf1tot']), std1, exp(pars['ksnf1std1']), exp(pars['nsnf1']) , exp(pars['ksnf1']), exp(pars['nsnf2'])) for std1 in np.linspace(std1min, std1max, tot)]	
	resp[j, :]=[(exp(pars["emig1max"])*snf1*1)/(exp(pars["kmig1snf1"]) + 1) for snf1 in s]
plt.imshow(resp, cmap='gnuplot2', extent=[std1min, std1max, glucmin, glucmax], aspect='auto'); plt.show()   
title('Mig1 deloc. rate') 
plt.xlabel('Std1')
plt.ylabel('glucose') 


######mig2 degradation rate as dependent on  std1 and glucose
c+=1;
plt.subplot(nrow, numplots,c)
mth1min=0
mth1max=50
std1min=0
std1max=1.5
glucmin=0
glucmax=1
resp=np.zeros([tot, tot])*np.nan
mth1=np.linspace(mth1max, mth1min, tot)
glus=np.linspace(glucmax, glucmin, tot)
for j in range(0, len(glus)):
	#1%glucose
	s=[fun(glus[j], exp(pars['Snf1tot']), std1, exp(pars['ksnf1std1']), exp(pars['nsnf1']) , exp(pars['ksnf1']), exp(pars['nsnf2'])) for std1 in np.linspace(std1min, std1max, tot)]	
	resp[j, :]=[hillinput(exp(pars["dmig2snf1"]),snf1, exp(pars["kmig2snf1"]), exp(0)) for snf1 in s]
plt.imshow(resp, cmap='gnuplot2', extent=[std1min, std1max, glucmin, glucmax], aspect='auto'); plt.show()   
title('Mig2 deg. rate') 
plt.xlabel('Std1')
plt.ylabel('glucose') 



###repressor inactivation as a function of glucose

gsteps=200000
c+=1;
plt.subplot(nrow, numplots,c)
glucmin=0
glucmax=1

glus=np.linspace(glucmin, glucmax, gsteps)
degmth1snf3=[hillinput(exp(pars["dmth1snf3"]),glu, exp(pars["k3"]), exp(pars["nmth1snf3"]))  for glu in glus]
degmth1rgt2=[hillinput(exp(pars["dmth1rgt2"]),glu, exp(pars["k2"]), exp(pars["nmth1rgt2"]))  for glu in glus]
exitstd1rgt2=[hillinput(exp(pars["estd1max"]),glu, exp(pars["k2"]), exp(pars["nstd1"]))  for glu in glus]
deghxt4=[hillrep(exp(pars["dhxt4max"]), glu, exp(pars["kdhxt4"]),exp(pars["ndhxt4"])) for glu in glus]
plt.figure; plt.plot(glus, degmth1snf3, label='Mth1 deg via Snf3') 
plt.plot(glus, degmth1rgt2, label='Mth1 deg via Rgt2') 
plt.plot(glus, exitstd1rgt2, label='Std1 nuclear exit via Rgt2')  
plt.plot(glus, deghxt4, label='Hxt4 deg via glucose')
plt.xlabel('glucose')                                                                                                                                                                                                                                                
plt.ylabel('inactivation rate')                                                                                                                            
plt.title('glucose inactivation') 
plt.legend( loc='upper right', bbox_to_anchor=(2.5, 0.3))       
#`figure.bbox`

    