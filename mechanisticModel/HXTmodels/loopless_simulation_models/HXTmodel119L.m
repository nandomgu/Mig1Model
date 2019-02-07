
function modelh= HXTmodel119L(params, mediaInput)
%model 114. third model (after model 112 and 113) that includes 2 sensors processing the glucose plus
%most of the items from model 111, minus threshold induced degradation of
% HXT
%input. In particular, hthis model also loses threshold/induced activation
%of mRNA.
%Model 115. reduces the above model by temoving the Hill functions and adds
%degradation of mth1 and std1. in this case, std1 and mth1 are produced
%constitutively and degraded by glucose to different extents through snf3
%and rgt2. mth1 and std1 actively repress the mRNA, each of them to
%different extents. 
%This model includes the repression of mth1 and std1 as hill functions.
%sensors are mm kinetics

modelh= @transcriptModelh;

function dout=transcriptModelh(t, x)
mediaInput=smooth(mediaInput);
Glucose=interp1(0:(numel(mediaInput)-1), mediaInput,t);


mRNA= x(1);
HXT= x(2);
MatureHXT= x(3);
Mth1= x(4);
Std1= x(5);

maxVmRNA=params(1);
KdegmRNA=params(2);
alphaHXT=params(3);
maxdegHXT= params(4);
KmatHXT=params(5);
KBleach=params(6); 
degThreshHXTbelow= params(7);
  KSnf3=params(8);
  KRgt2=params(9);
  sensorRatio=params(10);
VMth1= params(11);
maxdegMth1= params(12);
VStd1= params(13);
maxdegStd1= params(14);
RepThreshMth1=params(15);
RepThreshStd1=params(16);
hillMth1=params(17);
hillStd1=params(18);
         
         
 
snf3Gate= 1/(1+(Glucose/KSnf3));
rgt2Gate= (Glucose)/(KRgt2+ Glucose);
DMth1= VMth1-  (sensorRatio*rgt2Gate* maxdegMth1)*Mth1-(1-sensorRatio*snf3Gate* maxdegMth1)*Mth1;
DStd1= VStd1- (sensorRatio*rgt2Gate* maxdegStd1)*Std1-(1-sensorRatio*snf3Gate* maxdegStd1)*Std1;


VmRNASum= maxVmRNA/(1+(Mth1*hillMth1/RepThreshMth1*hillMth1))+ maxVmRNA/(1+(Std1^hillStd1/RepThreshStd1^hillStd1));

 DmRNA= VmRNASum-KdegmRNA*mRNA;
 DHXT= alphaHXT*mRNA-((maxdegHXT)/(degThreshHXTbelow+Glucose))*HXT;

  DMatureHXT= KmatHXT;
  
  
 dout= zeros(5,1);
 dout(1)=DmRNA;
 dout(2)=DHXT;
 dout(3)=DMatureHXT*HXT-(KBleach*MatureHXT);
 dout(4)= DMth1;
 dout(5)= DStd1;
 
end 
end