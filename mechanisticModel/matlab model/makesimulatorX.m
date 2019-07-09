function handle=makesimulatorX(modelName, params )
%make makesimulator 3 is for mechanistic models, 
%just as make simulator 2, but it allows to
%do the fitting only on a subset of parameters while keeping
%others fixed.
%makesimulator 4 is the same as before but this time argsim
%can contain the initial conditions of the model.

input=params.input;
data=params.data;
lsq=params.lsq;
mig1=params.mig1;
%data is a m replicate x n timepoint matrix
%% run specific params
modelFeatures=extractModelFeatures(modelName, [], 0); %making non negative
times= linspace(0,20, 250);
times=times(1:numel(input));
params.times= times;

func=eval(['@' modelName ';']);

handle= @rampSim;
if lsq
    outfun=@(y,d) sum((y-d).^2) ;
else
    outfun=@(y,d) y-d;
end
%establishing initial conditions from data
params.initialconditions= params.datameans(1, :);
initinds=1:numel(params.initialconditions);
params.initialconditions= [params.initialconditions repmat([1,0,0,0], 1, 18)];
function [lsqdiff, tf,yf]=rampSim(pars) 
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

sim=func(pars,params);

[tf, yf]=ode113(sim, times, params.initialconditions.*(params.initialconditions>0), params.opts);

yf=real(yf);
yout=yf(:, params.outvar);

% try
% yf=topUp(real(yf), numel(d), Inf)';
% catch
% yf=Inf(numel(d), 1);
% end
% 
% lsqdiff=outfun(yf,d);

lsqdiff= nansum(nansum((yf(:, initinds) -params.datameans  ).^2));


end

end





% plot(tf, yf, 'Color', cols(c,:), 'LineWidth', 1.5, 'DisplayName', num2str(pars)); 
% hold on
% 
% plot(repmat(ntimes, size(meandata.hxt1.g1percent,1), 1)', meandata.hxt1.g1percent')
% yyaxis right; plot(ntimes, cy5(n,:))
% title(strjoin([modelName, ': ' ,modelFeatures.paramNames(num)]))

