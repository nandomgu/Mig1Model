
function modelh= HXTmodelFFL2L(params, mediaInput)
%%Model FFL2L is an incoherent type 2 FFL. Glucose destroys Mth1 which is
%%repressing the Hxt. on the other hand glucose directly inhibits the hxt through mig1. 
%%Hill Hypothesis, and threshold induced repression of Mig1.


modelh= @transcriptModelh;

function dout=transcriptModelh(t, x)

mediaInput=smooth(mediaInput);
Glucose=interp1(0:(numel(mediaInput)-1), mediaInput,t);

mRNA= x(1);
HXT= x(2);
Mig1=x(3);
Mth1=x(4);


VmRNA=params(1);
hillmRNA=params(2);
KdegmRNA=params(3);
alphaHXT=params(4);
KdegHXT= params(5);
alphaMig1=params(6);
KdegMig1=params(7);
KactiveMig1=params(8);
RepThreshMig1=params(9);
alphaMth1=params(10);
maxdegMth1=params(11);
degThreshMth1=params(12);
RepThreshMth1=params(13);
hilldegMth1=params(14);
K=params(15);

 
 DMig1= alphaMig1-KdegMig1*Mig1;
 activeMig1= (Mig1*Glucose^4)/(KactiveMig1+ Glucose^4);
 DmRNA= VmRNA/(1+(activeMig1/RepThreshMig1)+(Mth1^hillmRNA/RepThreshMth1^hillmRNA) )+ -KdegmRNA*mRNA;
 DHXT= alphaHXT*mRNA- KdegHXT*HXT;
 
 DMth1= alphaMth1 - (K-(maxdegMth1/(1+(Glucose^hilldegMth1/degThreshMth1^hilldegMth1)))) *Mth1;
 
 dout= zeros(4,1);
 dout(1)=DmRNA;
  dout(2)=DMig1;
 dout(3)=DHXT;
 dout(4)=DMth1;

 
end 
end