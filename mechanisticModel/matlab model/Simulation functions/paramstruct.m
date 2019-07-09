function s=paramstruct(modelname, paramvalues)
%% generate a structure of parameters from  an array of parameters
%%modelname is the original model to which the param order corresponds.
mf=extractModelFeatures(modelname);
%s=defaultparams(modelname);
s=struct;
for j=1:numel(mf.paramNames)
   s.(mf.paramNames{j})=paramvalues(j);
    
    
end