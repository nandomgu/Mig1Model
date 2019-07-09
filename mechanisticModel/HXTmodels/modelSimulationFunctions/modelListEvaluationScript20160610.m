allModels={'HXTmodel110L'...
'HXTmodel111L' 'HXTmodel112L' 'HXTmodel113L' 'HXTmodel114L'...
'HXTmodel115L' 'HXTmodel116L' 'HXTmodel117L' 'HXTmodel118L'...
'HXTmodel119L' 'HXTmodel120L' 'HXTmodel121L' 'HXTmodel300corL'...
'HXTmodel300L' 'HXTmodel301L' 'HXTmodel302L' 'HXTmodel400L'...
'HXTmodel401L' 'HXTmodel601L' 'HXTmodel602L' 'HXTmodel603L'...
'HXTmodel603bL''HXTmodel604L' 'HXTmodel605L' 'HXTmodel606L'...
'HXTmodel607L' 'HXTmodel608L' 'HXTmodel703L' 'HXTmodel1161L'...
'HXTmodel1162L' 'HXTmodel1163L' 'HXTmodelFFL1L' 'HXTmodelFFL1TF'...
'HXTmodelFFL1TF2' 'HXTmodelFFL2L' 'HXTmodelFFL3L''HXTmodelFFL4L'...
'HXTmodelFFL6L' 'HXTmodelFFL7L' 'HXTmodelFFL8L' 'HXTmodelFFL12L'...
'HXTmodelFFL13L' 'HXTmodelFFL32L' 'HXTmodelFFL121L'};



for j=8:numel(allModels)
    allModels{j}
mf=extractModelFeatures(allModels{j});
expSim = multiExpSimulationCor(allModels{j}, 3, zeros(1,mf.numInit) ,0, [.2 .4 1 1 1],1, multichamber20160223.hxt4, multichamber20160301.hxt4,multichamber20160206.hxt4, multichamber20160117.hxt4, multichamber20160123.hxt4)
[fparams fval]= fminsearch(expSim,log(exp(normrnd(0,2,1,numel(mf.paramNum)))), minsearchoptions);
end 

system(sprintf(['echo ''' '#' 'MODEL_NAME' '\t' 'PARAMETERS' '\t' 'INITIAL_CONDITIONS'  '\t' 'CORRELATION_COST' '\t' 'LSQDIFF_COST' '''>> ' reverseDate '_modelExplorations.txt']));
for j=1:numel(allModels)
model=allModels{j};
expSim = multiExpExploreCor(model, 3, [0 0 0],0, [.2 .4 1 1 1],1, multichamber20160223.hxt4, multichamber20160301.hxt4,multichamber20160206.hxt4, multichamber20160117.hxt4, multichamber20160123.hxt4)
mf=extractModelFeatures(model); 
pars=generateRandomParams(model,10,[]);
pars2=generateRandomParams(model,10, 8);
pars= [pars; pars2];
expSim(pars);
end

%%%random running of models
for j=1:numel(allModels)
model=allModels{randsample(numel(allModels), 1)};
mf=extractModelFeatures(model); 
expSim = multiExpExploreCor(model, 3, repmat(0, 1, mf.numInit),0, [.2 .4 1 1 1],1, multichamber20160223.hxt4, multichamber20160301.hxt4,multichamber20160206.hxt4, multichamber20160117.hxt4, multichamber20160123.hxt4);
pars=generateRandomParams(model,10,[]);
pars2=generateRandomParams(model,10, 8);
pars= [pars; pars2];
expSim(pars);
end