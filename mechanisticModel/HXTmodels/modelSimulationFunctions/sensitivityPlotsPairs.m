function sensitivityPlotsPairs(modelname, params, paramNum,numIter,integer, initialConditions, outvar)
%%explores simulations of model modelname by systematically modifying pairs
%%of parameters in the set paramNum, starting from the starting set params. if no
%%paramNum is given, then it does pairs of all parameters.
modelFeatures= extractModelFeatures(modelname)

if isempty(paramNum)
    paramNum= 1:numel(params)
end

figure; 


if nargin<7
    initialConditions=[0 0 0 0];
end
    
    
    c=1
for j=1:numel(paramNum)
    disp(j)
    
   for ii=1:(numel(paramNum))
       
    
    c=c+1;
    if ii>=j
        continue
    else
        subplot(numel(paramNum), numel(paramNum), c);
        ylabel(modelFeatures.paramNames{paramNum(j)})
        title(modelFeatures.paramNames{paramNum(ii)})
        exploreParamComb(modelname, params, [paramNum(j) paramNum(ii)], numIter, integer, initialConditions, outvar);
    end
    
    
   end
end





