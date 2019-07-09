
function pars= defaultparams(model, paramvals)


paramvalsdef=struct;
paramvalsdef.('VHXT'  		)	=2.6528;
paramvalsdef.('Vloc'    		)	=8.3923;
paramvalsdef.('basaldeg'  	)	 =-1.1092;
paramvalsdef.('maxdegMT'   	)	= 0.8044;
paramvalsdef.('KdegMT'    	)	=0.0007;
paramvalsdef.('hilldegMT'   	)	=		4.1181;
paramvalsdef.('basaldegMT'  	)	=	-1.0028;
paramvalsdef.('maxdegST'   	)	= 0.8044;
paramvalsdef.('KdegST'    	)	=0.0007;
paramvalsdef.('hilldegST'   	)	=		4.1181;
paramvalsdef.('basaldegST'  	)	=	-1.0028;
paramvalsdef.('deloc'    	)	=	1.2254;
paramvalsdef.('KrepMT'   	)	=-2.8048;
paramvalsdef.('hillrepMT'   	)	=	2.8987;
paramvalsdef.('KrepHMG'   	)	=-4.0021;
paramvalsdef.('hillHMG'    	)	=0.8293;
paramvalsdef.('KrepHMT'   	)	=-0.5924;
paramvalsdef.('hillHMT'    	)	=1.3689;
paramvalsdef.('KMG'    		)	=-1.2;
paramvalsdef.('hillMG'    	)	=2.7869;
paramvalsdef.('VdegHXT'   	)	=-0.3587;
paramvalsdef.('threshdegHXT'	)	=  -5.0099;
paramvalsdef.('hilldegHXT'  	)	= 	1.3232;
paramvalsdef.('KinhSMG')	= 	-.5;
paramvalsdef.('KrepHST')	= 	-2;
paramvalsdef.('hillrepHST')	= 	2;
paramvalsdef.('KrepMTST')	= 	-4;
paramvalsdef.('hillrepMTST')	= 	2;
paramvalsdef.('basaldegM2')	= 	.05;
paramvalsdef.('KrepSTM2')	= 	log(.5);
paramvalsdef.('hillrepSTM2')	= 	log(2);
paramvalsdef.('KrepHM2')	= 	2;
paramvalsdef.('hillrepHM2')	= 	2;
paramvalsdef.('KinhSM2')	= 	100;
paramvalsdef.('KrepM2MT')	= 	.04;
paramvalsdef.('hillrepM2MT')	= 	4;
paramvalsdef.('VinactGM2')= log(2);
paramvalsdef.('KinactGM2')= log(.3);
paramvalsdef.('hillinactGM2')=log(3);


defval=1;
pars=[];
mf=extractModelFeatures(model);


for j =1:numel(mf.paramNames)
   %disp(mf.paramNames{j})
    %if the parameter is in the new params, use it
    if nargin>1 && isfield(paramvals, mf.paramNames{j})
        %disp([ mf.paramNames{j} 'is in the list of custom parameters'])
        %pause(2)
    pars(j)= paramvals.(mf.paramNames{j});
    %%if the parameter is not in the new set but there is  a dfault, use
    %%default
    else if (nargin <2 && isfield(paramvalsdef, mf.paramNames{j})) || (nargin >1 && ~isfield(paramvals, mf.paramNames{j}))
         disp([ 'there is a predefined value  for ' mf.paramNames{j} '.will use ' num2str(paramvalsdef.(mf.paramNames{j}))])

            pars(j)= paramvalsdef.(mf.paramNames{j}); 
        else % if it is nowhere, just default the parameter to 1 
             disp([ 'there is NO predefined value  for ' mf.paramNames{j} '.will use ' num2str(defval)])
            
            pars(j)= defval;
        end
    end
    
        
end
end
    

