function pars=generateRandomParams(model,num, hillTop)
%%this function generates a long list of random parameters specific for
%%model. to do so it uses extractModelFeatures.
%%%hillTop generates a 
modelFeatures= extractModelFeatures(model);

pars=normrnd(0,2,num,numel(modelFeatures.paramNum));

if ~isempty(hillTop)
    
for j = 1:numel(modelFeatures.hillIndices)
    
    
    pars(:, modelFeatures.hillIndices(j))= log(randsample(hillTop,num, true));

end

end

end


