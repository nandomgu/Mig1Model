function s=paramstruct(modelname, paramvalues)
%% generate a structure of parameters from a list
s=struct;
mf=extractModelFeatures(modelname);
for j=1:numel(mf.paramNames)
   s.(mf.paramNames{j})=paramvalues(j);
    
    
end