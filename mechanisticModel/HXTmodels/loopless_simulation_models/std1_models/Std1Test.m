function modelh= Std1Test(params, mediaInput)
%%this model is of the negative self regulation of protein Std1 
%%%%% Std1 protein self repressed, glucose-induced degradation and basal degradation too.
%%this model also includes mth1 as a basic protein that gets degraded with glucose.
%%this model includes active mig1 as going in and out of the nucleus, 
%%but Std1 raises the threshold of activation of Std1.
%%at first glance, this model adds some extra parameters to account directly
%%for the Std1 effect on mig1.

modelh= @transcriptModelh;

function dout=transcriptModelh(t, x)

mediaInput=smooth(mediaInput);
Glucose=interp1(0:(numel(mediaInput)-1), mediaInput,t);


Std1=x(1);
Mth1=x(2);
activeMig1=x(3);


VStd1=params(1);
KSelfStd1=params(2);
KDegStd1=params(3);
threshStd1=params(4);
basalDegStd1=params(5);
VMth1= params(6);
basalDegMth1=params(7);
KDegMth1=params(8);
threshMth1=params(9);
KRepMth1=params(10);
thresholdStd1Mig1=params(11); %%threshold concentration of Std1 at which its mig1 effect is noticeable.
hillactiveMig1=params(12);

hillStd1Mig1=2; %%hill factor for the std1-mig1 effect
%%equation that describes how the mig1 activtion threshold gets affected by std1
KactiveMig1= @(Std1,thresholdStd1Mig1, hillStd1Mig1)  Std1^hillStd1/(thresholdStd1Mig1^hillStd1+Std1^hillStd1Mig1);


DStd1= VStd1/(1+(Std1/KSelfStd1)+(Mth1/KRepMth1))- basalDegStd1*Std1 -(KDegStd1*Glucose/ (threshStd1+Glucose)*Std1); 
DMth1=  VMth1- basalDegMth1*Mth1- ((KDegMth1*Glucose)/(threshMth1+Glucose))*Mth1;

%active state of mig1 is determined by glucose and std1
activeMig1= (Glucose^hillactiveMig1)/(KactiveMig1(Std1)^hillactiveMig1+ Glucose^hillactiveMig1)

dout= zeros(2,1);
dout(1)=DStd1;
dout(2)=DMth1;
dout(3)=activeMig1;

end
end


%%% 