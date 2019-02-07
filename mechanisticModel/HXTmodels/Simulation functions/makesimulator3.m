function handle=makesimulator3(modelName,params )
%make makesimulator 3 is for mechanistic models, 
%just as make simulator 2, but it allows to
%do the fitting only on a subset of parameters while keeping
%others fixed
input=params.input;
data=params.data;
lsq=params.lsq;
mig1=params.mig1;
%data is a m replicate x n timepoint matrix
%% run specific params
modelFeatures=extractModelFeatures(modelName);
times= linspace(0,20, 250);
times=times(1:numel(input));
params.times= times;


func=eval(['@' modelName ';'])
handle= @rampSim;
if lsq
    outfun=@(y,d) sum((y-d).^2) ;
else
    outfun=@(y,d) y-d;
end

function [lsqdiff, tf,yf,d]=rampSim(pars) 
%
temppars=pars; %save temp copy of parameters
pars=params.defaultparams; %make all parameter values default
%now plug in the values of the parameters we care about,
%whose indices are in params.onlyparams
%this way we provide genearal paramsfor,  eg. mth1 and mig1
%but leave the other parameters free.
%to leave all parameters free, just make onlyparams=1:16
pars(params.onlyparams)=temppars(params.onlyparams);
pars=exp(pars);
d=nanmean(data)';
sim=func(pars,params);
[tf, yf]=ode23(sim, times, [d(1), 1,0], modelFeatures.options);

%s=1:size(yf(:, params.outvar),1);
d=d(1:numel(times));
yf=yf(:, params.outvar);
lsqdiff=outfun(yf,d);

end

end





% plot(tf, yf, 'Color', cols(c,:), 'LineWidth', 1.5, 'DisplayName', num2str(pars)); 
% hold on
% 
% plot(repmat(ntimes, size(meandata.hxt1.g1percent,1), 1)', meandata.hxt1.g1percent')
% yyaxis right; plot(ntimes, cy5(n,:))
% title(strjoin([modelName, ': ' ,modelFeatures.paramNames(num)]))

