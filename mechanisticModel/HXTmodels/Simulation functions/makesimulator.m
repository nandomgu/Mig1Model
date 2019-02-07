function handle=makesimulator(modelName,params )

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
pars=exp(pars);
d=nanmean(data)';
sim=func(pars,params);
[tf, yf]=ode23(sim, times,d(1) , modelFeatures.options);


d=d(1:numel(times));
lsqdiff=outfun(yf,d);

end

end





% plot(tf, yf, 'Color', cols(c,:), 'LineWidth', 1.5, 'DisplayName', num2str(pars)); 
% hold on
% 
% plot(repmat(ntimes, size(meandata.hxt1.g1percent,1), 1)', meandata.hxt1.g1percent')
% yyaxis right; plot(ntimes, cy5(n,:))
% title(strjoin([modelName, ': ' ,modelFeatures.paramNames(num)]))

