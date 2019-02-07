
function modelh= HXTmodel113L(params, mediaInput)
%model 113. First model that includes 2 sensors processing the glucose plus
%most of the items from model 111, minus threshold induced degradation of
% HXT
%input

modelh= @transcriptModelh;

function dout=transcriptModelh(t, x)
mediaInput=smooth(mediaInput);
Glucose=interp1(0:(numel(mediaInput)-1), mediaInput,t);

mRNA= x(1);
HXT= x(2);
MatureHXT= x(3);


maxVmRNA=params(1);
KmRNA=params(2);
hillmRNA=params(3);
KdegmRNA=params(4);
alphaHXT=params(5);
maxdegHXT= params(6);
KmatHXT=params(7);
KBleach=params(8);
    RepThreshmRNA=params(9);
      hillRepmRNA=params(10);
degThreshHXTbelow= params(11);
  hillDegHXTbelow= params(12);
      sensorRatio=params(13);
            KRgt2=params(14);
         hillRgt2=params(15);
            KSnf3=params(16);
         hillSnf3=params(17);


 
snf3Gate= 1/(1+(Glucose/KSnf3)^hillSnf3);
rgt2Gate= (Glucose^hillSnf3)/(KRgt2^hillRgt2+ Glucose^hillRgt2);

VmRNASum= (maxVmRNA* sensorRatio*rgt2Gate)+(maxVmRNA* (1-sensorRatio)*snf3Gate);

 DmRNA= (VmRNASum/(1+(Glucose/RepThreshmRNA)^hillRepmRNA))*(Glucose)^hillmRNA /(KmRNA^hillmRNA+Glucose^hillmRNA)-KdegmRNA*mRNA;
 DHXT= alphaHXT*mRNA-((maxdegHXT^hillDegHXTbelow)/(degThreshHXTbelow^hillDegHXTbelow+Glucose^hillDegHXTbelow))*HXT;

  DMatureHXT= KmatHXT;
  
  
 dout= zeros(3,1);
 
 dout(1)=DmRNA;
 dout(2)=DHXT;
 dout(3)=DMatureHXT*HXT-(KBleach*MatureHXT);
 
 
end 
end