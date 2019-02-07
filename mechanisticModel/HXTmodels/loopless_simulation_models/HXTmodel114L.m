
function modelh= HXTmodel114L(params, mediaInput)
%model 114. third model (after model 112 and 113) that includes 2 sensors processing the glucose plus
%most of the items from model 111, minus threshold induced degradation of
% HXT
%input. In particular, hthis model also loses threshold/induced activation
%of mRNA

modelh= @transcriptModelh;

function dout=transcriptModelh(t, x)
mediaInput=smooth(mediaInput);
Glucose=interp1(0:(numel(mediaInput)-1), mediaInput,t);

mRNA= x(1);
HXT= x(2);
MatureHXT= x(3);


maxVmRNA=params(1);
KdegmRNA=params(2);
alphaHXT=params(3);
maxdegHXT= params(4);
KmatHXT=params(5);
KBleach=params(6);
    RepThreshmRNA=params(7);
      hillRepmRNA=params(8);
degThreshHXTbelow= params(9);
  hillDegHXTbelow= params(10);
      sensorRatio=params(11);
            KRgt2=params(12);
         hillRgt2=params(13);
            KSnf3=params(14);
         hillSnf3=params(15);


 
snf3Gate= 1/(1+(Glucose/KSnf3)^hillSnf3);
rgt2Gate= (Glucose^hillSnf3)/(KRgt2^hillRgt2+ Glucose^hillRgt2);

VmRNASum= (maxVmRNA* sensorRatio*rgt2Gate)+(maxVmRNA* (1-sensorRatio)*snf3Gate);

 DmRNA= (VmRNASum/(1+(Glucose/RepThreshmRNA)^hillRepmRNA))-KdegmRNA*mRNA;
 DHXT= alphaHXT*mRNA-((maxdegHXT^hillDegHXTbelow)/(degThreshHXTbelow^hillDegHXTbelow+Glucose^hillDegHXTbelow))*HXT;

  DMatureHXT= KmatHXT;
  
  
 dout= zeros(3,1);
 
 dout(1)=DmRNA;
 dout(2)=DHXT;
 dout(3)=DMatureHXT*HXT-(KBleach*MatureHXT);
 
 
end 
end