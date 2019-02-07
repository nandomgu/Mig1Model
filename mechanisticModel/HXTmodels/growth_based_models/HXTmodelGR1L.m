function modelh= HXTmodelFFL12L(params, mediaInput)
%%Model FFL1 is simple threshold-activated production of Hxt, including the
%%Hill Hypothesis, and threshold induced repression of Mig1. This function
%%returns a handle for the ODE solver.


modelh= @transcriptModelh;

function dout=transcriptModelh(t, x)

mediaInput=smooth(mediaInput);
Glucose=interp1(0:(numel(mediaInput)-1), mediaInput,t);

mRNA= x(1);
HXT= x(2);



VmRNA=params(1);
KmRNA=params(2);
hillmRNA=params(3);
KdegmRNA=params(4);
alphaHXT=params(5);
KdegHXT= params(6);
KactiveMig1=params(7);
RepThreshMig1=params(8);
hillactiveMig1=params(9);



 activeMig1= (Glucose^hillactiveMig1)/(KactiveMig1^hillactiveMig1+ Glucose^hillactiveMig1);
 DmRNA= (VmRNA*(Glucose)^hillmRNA /(KmRNA^hillmRNA+Glucose^hillmRNA))/(1+(activeMig1/RepThreshMig1))+ -KdegmRNA*mRNA;
 DHXT= alphaHXT*mRNA- KdegHXT*HXT;
 
 
 dout= zeros(2,1);
 dout(1)=DmRNA;
 dout(2)=DHXT;

 
end 
end