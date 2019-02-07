function modelh= HXTmodelStd1L(params, mediaInput)
%%Model FFL1 is simple threshold-activated production of Hxt, including the
%%Hill Hypothesis, and threshold induced repression of Mig1. This function
%%returns a handle for the ODE solver.


modelh= @transcriptModelh;

function dout=transcriptModelh(t, x)

mediaInput=smooth(mediaInput);
Glucose=interp1(0:(numel(mediaInput)-1), mediaInput,t);

mRNA= x(1);
HXT= x(2);
Std1=x(3);


VmRNA=params(1);
KmRNA=params(2);
hillmRNA=params(3);
KdegmRNA=params(4);
alphaHXT=params(5);
KdegHXT= params(6);
KactiveMig1=params(7);
RepThreshMig1=params(8);
hillactiveMig1=params(9);
VStd1=1; %we fix max production of Std1 to 1
KSelfStd1=params(10)




 activeMig1= (Glucose^hillactiveMig1)/(KactiveMig1^hillactiveMig1+ Glucose^hillactiveMig1);
 
 %testing the above equation
%  Glucose=linspace(0,1, 10000);
%  am1= @(Glucose, hillactiveMig1, kactiveMig1) (Glucose^hillactiveMig1)/(KactiveMig1^hillactiveMig1+ Glucose^hillactiveMig1);
%  figure; plot(Glucose, arrayfun(am1, Glucose, repmat(hillactiveMig1, 1, numel(Glucose)), repmat(KactiveMig1, 1, numel(Glucose))))
%  xlabel('Glucose');
%  ylabel('activeMig1');
 
 DStd1= VStd1/ (1+Std1/KSelfStd1)- KDegStd1
 
 DmRNA= (VmRNA*(Glucose)^hillmRNA /(KmRNA^hillmRNA+Glucose^hillmRNA))/(1+(activeMig1/RepThreshMig1))+ -KdegmRNA*mRNA;
 DHXT= alphaHXT*mRNA- KdegHXT*HXT;
 
 
 dout= zeros(2,1);
 dout(1)=DmRNA;
 dout(2)=DHXT;

 
end 
end