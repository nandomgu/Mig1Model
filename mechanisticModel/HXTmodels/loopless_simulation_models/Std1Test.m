function modelh= Std1Test(params, mediaInput)
%%this model is of the negative self regulation of protein Std1 
%%%%% Std1 protein self repressed, glucose-induced degradation and basal degradation too.
%%this model also includes mth1 as a basic protein that gets degraded with glucose.
%%this model includes 2 0-1 equations, one of them snf1 and another one for active mig1.
%%snf1 deactivates mig1. the trick is though, that snf1 deactivation of mig1 is enhanced by std1.
%%this happens through additive repression

modelh= @transcriptModelh;

function dout=transcriptModelh(t, x)

%mediaInput=smooth(mediaInput);
Glucose=interp1(0:(numel(mediaInput)-1), mediaInput,t);


Std1=x(1);
Mth1=x(2);
Snf1=x(3);
activeMig1=x(4);
HXT1=x(5);
HXT2=x(6);
HXT3=x(7);
HXT4=x(8);
HXT7=x(9);
Hxt2=x(10);


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
KSnf1= params(11);
KRepSnf1Mig1=params(12);
KEnhance=params(13);
hillactiveMig1=params(14);
KMig1HXT4=params(15);
KStd1HXT4=params(16);
KMth1HXT4=params(17);
KMig1HXT2=params(18);
KMth1HXT2=params(19);
KDegSnf1=params(20);
basalDegSnf1=params(21);
basalDegHXT2=params(22);

%DStd1= VStd1/(1+(Std1/KSelfStd1)+(Mth1/KRepMth1))- basalDegStd1*Std1 -(KDegStd1*Glucose/ (threshStd1+Glucose)*Std1); 
DStd1= VStd1/(1+(Std1/KSelfStd1) +(Mth1/KRepMth1)) - basalDegStd1*Std1 -KDegStd1*Glucose*Std1;%-(KDegStd1*Glucose/ (threshStd1+Glucose)*Std1) ;       %-((KDegStd1*Glucose/(threshStd1+Glucose)^2)*Std1); 
DMth1=  VMth1- basalDegMth1*Mth1- ((KDegMth1*Glucose^4)/(threshMth1^4+Glucose^4))*Mth1;
%active state of mig1 is determined by glucose and Snf1 and Std1

%%playing with the phosphatase FFL. reg1 dephosphorylates snf1
DSnf1= 1-0.05*Snf1-.05*Snf1;


%.05+.015*Snf1+Glucose/(.01+Glucose)  - 0.04*Snf1;
%1/(1+(Glucose/.1))-.1*Snf1; %%we could try  keeping an activation thresh. constant or to have the threshold increase with Std1
%activeMig1= 1/ (1+ (Snf1/KRepSnf1Mig1));

%glucose induces activation of mig1, and is inhibited by Snf1, only in the
%presence of Std1
DactiveMig1= Glucose/(0.01+Glucose)  -.1*activeMig1; %-.1*activeMig1; %-.5*Snf1*activeMig1;  %*Std1;
DHXT1= 1/(1+ Mth1/.0001)- .1*HXT1; %mth1 activation only
DHXT2= 1/(1+((Mth1/5)^4)+(activeMig1/10)^4)-basalDegHXT2*HXT2;
DHXT3= 1/(1+ Mth1/.001)- .1*HXT3;
DHXT4= +1/(1+(Mth1/KMth1HXT4)+(activeMig1/KMig1HXT4))- .05*HXT4;
DHXT7=  1/(1+ activeMig1/.1)- .1*HXT7;

DHxt2= .3*HXT4- .3*HXT4;


dout= zeros(10,1);
dout(1)=DStd1;
dout(2)=DMth1;
dout(3)=DSnf1;
dout(4)=DactiveMig1;
dout(5)=DHXT1;
dout(6)=DHXT2;
dout(7)=DHXT3;
dout(8)=DHXT4;
dout(9)=DHXT7;
dout(10)=DHxt2;
end
end




