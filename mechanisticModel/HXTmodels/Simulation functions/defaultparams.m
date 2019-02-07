function pars= defaultparams(model)
paramvals=struct;
paramvals.('VHXT'  		)	=2.6528;
paramvals.('Vloc'    		)	=8.3923;
paramvals.('basaldeg'  	)	 =-1.1092;
paramvals.('maxdegMT'   	)	= 0.8044;
paramvals.('KdegMT'    	)	=0.0007;
paramvals.('hilldegMT'   	)	=		4.1181;
paramvals.('basaldegMT'  	)	=	-1.0028;
paramvals.('maxdegST'   	)	= 0.8044;
paramvals.('KdegST'    	)	=0.0007;
paramvals.('hilldegST'   	)	=		4.1181;
paramvals.('basaldegST'  	)	=	-1.0028;
paramvals.('deloc'    	)	=	1.2254;
paramvals.('KrepMT'   	)	=-2.8048;
paramvals.('hillrepMT'   	)	=	2.8987;
paramvals.('KrepHMG'   	)	=-4.0021;
paramvals.('hillHMG'    	)	=0.8293;
paramvals.('KrepHMT'   	)	=-0.5924;
paramvals.('hillHMT'    	)	=1.3689;
paramvals.('KMG'    		)	=0.4901;
paramvals.('hillMG'    	)	=2.7869;
paramvals.('VdegHXT'   	)	=-0.3587;
paramvals.('threshdegHXT'	)	=  -5.0099;
paramvals.('hilldegHXT'  	)	= 	1.3232;
paramvals.('KinhSMG')	= 	-.5;
paramvals.('KrepHST')	= 	-2;
paramvals.('hillHST')	= 	2;

pars=[];
mf=extractModelFeatures(model);
for j =1:numel(mf.paramNames)
    pars(j)= paramvals.(mf.paramNames{j});
end
    

