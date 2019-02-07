
function modelh= HXTmodelFFL32L(params, mediaInput)
%%Model FFL3L is an incoherent type 2 FFL. Glucose destroys Mth1 which is
%%repressing the Hxt. on the other hand glucose directly inhibits the hxt through mig1. 
%%Hill Hypothesis, and threshold induced repression of Mig1. The only
%%difference from FFL2L is that this one implements glucose controled HXT
%%degradation.
%preferedOutVar=3;
%preferedInitialConditions= [0 0 0 1]

modelh= @transcriptModelh;

function dout=transcriptModelh(t, x)

mediaInput=smooth(mediaInput);
Glucose=interp1(0:(numel(mediaInput)-1), mediaInput,t);

mRNA= x(1);
HXT= x(2);


VmRNA=params(1);
hillmRNA=params(2);
KdegmRNA=params(3);
alphaHXT=params(4);
KdegHXT= params(5);
KactiveMig1=params(6);
RepThreshMig1=params(7);
degThreshMth1=params(8);
RepThreshMth1=params(9);
hilldegMth1=params(10);
hillactiveMig1=params(11);
hillRepMth1=params(12);


 Mig1= (Glucose^hillactiveMig1)/(KactiveMig1^hillactiveMig1+ Glucose^hillactiveMig1);
 Mth1= 1/ (1+ (Glucose^hilldegMth1/degThreshMth1^hilldegMth1));
 DmRNA= VmRNA/(1+(Mig1^hillmRNA/RepThreshMig1^hillmRNA)+(Mth1^hillRepMth1/RepThreshMth1^hillRepMth1) ) -KdegmRNA*mRNA;
 DHXT= alphaHXT*mRNA- KdegHXT*HXT;

 dout= zeros(2,1);
 dout(1)=DmRNA;
 dout(2)=DHXT;

 
end 
end