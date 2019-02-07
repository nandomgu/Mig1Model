
function modelh= HXTmodel400L(params, mediaInput)
%model 114. third model (after model 112 and 113) that includes 2 sensors processing the glucose plus
%most of the items from model 111, minus threshold induced degradation of
% HXT
%input. In particular, hthis model also loses threshold/induced activation
%of mRNA.
%Model 115. reduces the above model by temoving the Hil functions and adds
%degradation of mth1 and std1. in this case, std1 and mth1 are produced
%constitutively and degraded by glucose to different extents through snf3
%and rgt2. mth1 and std1 actively repress the mRNA, each of them to
%different extents. 
%Model 117. model 116 gets simplified parameter-wise, and as repression can
%only by maximal (vmax=0) then the relative contributions to repression by
%Mth1 and std1 can be defined as relative to each other. 

%Model1162 has mth1 and std1 , and their degradation dependent on hil
%functions. 
%model 1163 includes snf3 and rgt2 as gates

%model300 has a correction from model 1163. the repression to the max mrna
%rate occurs in the same equation such that the rate cannot be maximum if
%there is even a little bit of any of the two repressors.


%in model 300 the glucose signal directly causes degradation of mth1 and
%std1, and the magnitude of that degradation is decided by the sensor
%gates. This was perceived as redundant but model 300 is able to achieve
%quite good fits to hxt4 data.

%model 300 corrected corrects the above redundancy, making the sensor gates
%the sole input for the degradation of mth1 and std1. for some reason doing
%so makes the hxt4 phenotype unreachable.

%model 400 implements yck1 and yck2 in the game, such that the yck/sensor
%interaction activates the degradation of hxt4. 



modelh= @transcriptModelh;

function dout=transcriptModelh(t, x)

mediaInput=smooth(mediaInput);
Glucose=interp1(0:(numel(mediaInput)-1), mediaInput,t);

mRNA= x(1);
HXT= x(2);
MatureHXT= x(3);
Mth1=x(4);
Std1=x(5);
Snf3= 1;
Rgt2=1;
Yck1=1;
Yck2=1;
Yck2Snf3=0;
Yck1Rgt2=0;


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
hillDegMth1=params(17);
hillDegStd1=params(18);
hillRepMig1=params(19);
hillSnf3=params(20);
hillRgt2=params(21);
KSnf3=params(22);
KRgt2=params(23);
VYck1=params(24);
VYck2=params(25);
degYck1=params(26);
degYck2=params(27);
VSnf3=params(28);
VRgt2=params(29);
degSnf3=params(30);
degRgt2=params(31);
koffYck1Rgt2=params(32);
koffYck2Snf3=params(33);
degMth1= params(34);
degStd1= params(35);


konYck2Snf3= 1/(1+(Glucose/KSnf3)^hillSnf3);
konYck1Rgt2= (Glucose)/(KRgt2+ Glucose);
DSnf3=VSnf3-degSnf3*Snf3;
DRgt2=VRgt2-degRgt2*Rgt2;

DYck1=VYck1-degYck1-konYck1Rgt2*Yck1*Rgt2+koffYck1Rgt2*Yck1Rgt2;
DYck2=VYck2-degYck2-konYck2Snf3*Yck2*Snf3+koffYck2Snf3*Yck2Snf3;


DYck1Rgt2= konYck1Rgt2*Yck1*Rgt2-koffYck1Rgt2*Yck1Rgt2;
DYck2Snf3= konYck2Snf3*Yck2*Snf3-koffYck2Snf3*Yck2Snf3;

DMth1= VMth1-  ((maxdegMth1/(1+(Yck2Snf3^hillDegMth1/KdegMth1^hillDegMth1)))+degMth1)*Mth1;
DStd1= VStd1-  ((maxdegStd1*(Yck1Rgt2^hillDegStd1)/(KdegStd1^hillDegStd1+Yck1Rgt2^hillDegStd1))+degStd1)*Std1;

VmRNASum= maxVmRNA/(1+(Mth1/RepThreshMth1)+(Std1/RepThreshStd1)+(Glucose^hillRepMig1/RepThreshMig1^hillRepMig1));

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