
function modelh= HXTmodel161L(params, mediaInput)
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
%Model 117. model 116 gets simplified parameter-wise, and as repression can
%only by maximal (vmax=0) then the relative contributions to repression by
%Mth1 and std1 can be defined as relative to each other. 

modelh= @transcriptModelh;

function dout=transcriptModelh(t, x)

mediaInput=smooth(mediaInput);
Glucose=interp1(0:(numel(mediaInput)-1), mediaInput,t);


mRNA= x(1);
HXT= x(2);
MatureHXT= x(3);
Mth1=x(4);
Std1=x(5);

maxVmRNA=params(1);
KdegmRNA=params(2);
alphaHXT=params(3);
maxdegHXT= params(4);
KmatHXT=params(5);
KBleach=params(6); 
degThreshHXTbelow= params(7);
  KdegMth1=params(8);
  KdegStd1=params(9);
VMth1= params(10);
maxdegMth1= params(11);
VStd1= params(12);
maxdegStd1= params(13);
RepThreshMth1=params(14);
RepThreshStd1=params(15);
RepThreshMig1=params(16);
         
         
DMth1= VMth1-  (maxdegMth1/1+(Glucose/KdegMth1))*Mth1;
DStd1= VStd1-  (maxdegStd1*(Glucose/KdegStd1+Glucose))*Std1;


VmRNASum= maxVmRNA/(1+(Mth1/RepThreshMth1))+ maxVmRNA/(1+(Std1/RepThreshStd1))+maxVmRNA/(1+(Glucose/RepThreshMig1));

 DmRNA= VmRNASum-KdegmRNA*mRNA;
DHXT= alphaHXT*mRNA-(maxdegHXT/1+(degThreshHXTbelow/(Glucose)))*HXT;

  DMatureHXT= KmatHXT;
  
  
 dout= zeros(5,1);
 dout(1)=DmRNA;
 dout(2)=DHXT;
 dout(3)=DMatureHXT*HXT-(KBleach*MatureHXT);
 dout(4)=DMth1;
 dout(5)= DStd1;
 
 
end 
end