import json
from scipy import interpolate
locdata=json.load(open('stepsdata.json'))  
t=np.around(np.linspace(0,20,250),3);
allvars=['cy5', 'times', 'centeredtimes', 'type']
strains=list(locdata.keys())
for j in allvars:
	try:
		strains.remove(j) ;
	except:
		print('something went wrong while getting rid of '+j)
		continue
formatnan= lambda x: np.nan if x=='_NaN_' else float(x)
for strain in strains:
	for conc in list(locdata[strain].keys()):
		for rep in list(locdata[strain][conc].keys()):
			locdata[strain][conc][rep]=np.array([[formatnan(j) for j in locdata[strain][conc][rep][y]] for y in range(0, np.size(locdata[strain][conc][rep],0))]);
plt.figure; plt.plot(locdata[strain][conc][rep].T)    

#function to fill nans by interpolation
def fill_nan(A):
    '''
    interpolate to fill nan values
    '''
    inds = np.arange(A.shape[0])
    good = np.where(np.isfinite(A))
    f = interpolate.interp1d(inds[good], A[good],bounds_error=False)
    B = np.where(np.isfinite(A),A,f(inds))
    return B
    
#cells that have less than 10 nans, we will try to inteprolate them
#plt.plot(locdata[strain][conc][rep][np.sum(~np.isnan(locdata[strain][conc][rep]),1)>240,:].T)



mat250=locdata[strain][conc][rep][np.sum(~np.isnan(locdata[strain][conc][rep]),1)==250,:]
mat240=locdata[strain][conc][rep][np.sum(~np.isnan(locdata[strain][conc][rep]),1)>240,:]
filled=[fill_nan(mat240[j, :]) for j in range(0, np.size(mat240, 0))]
#getting the median 
colmedian=np.nanmedian(mat240, 0)


##importing hxt mean data
cols=['grey' ,'orange','cyan', 'magenta', 'green', 'red', 'blue'];
hxtmeandata=json.load(open('hxtmeandata.json'))
strains=list(hxtmeandata.keys())
c=0;
plt.figure()
for strain in strains:
	for conc in list(hxtmeandata[strain].keys()):
		#hxtmeandata[strain][conc]=np.array(hxtmeandata[strain][conc])
	plt.plot(hxtmeandata[strain][conc].T, color=cols[c])     
	c+=1;	
	
	
	