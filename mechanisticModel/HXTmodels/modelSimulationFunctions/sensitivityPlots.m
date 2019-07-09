function sensitivityPlots(modelname, params, paramNum,numIter,integer, initialConditions, outvar)
%%explores simulations of model modelname with parameters params, when each
%%of the parameters paramNum is given random values using exploreParamComb.
%% we start with 16 fixed subplots for simplicity. 

modelFeatures= extractModelFeatures(modelname);
%figure; 


if nargin<7
    initialConditions=[0 0 0 0];
end
    
    
    
for j=1:numel(paramNum)
    disp(j)
    subplot(4, 4, j);
    
title(modelFeatures.paramNames{j})
exploreParamComb(modelname, params, paramNum(j), numIter, integer, initialConditions, outvar);

end





