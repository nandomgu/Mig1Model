
function modelh= HXTmodel8L(params, mediaInput)
%Model 8 will include glucose inhibition of Mth1 and Mth1 inhibition of
%Mig1, and Mig1 inhibition of transscription. 

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
degThreshHXT=params(9);
hillDegHXT=params(10);
RepThreshmRNA=params(11);
hillRepmRNA=params(12);

  
 DmRNA= (maxVmRNA/(1+(Glucose/RepThreshmRNA)^hillRepmRNA))*(Glucose)^hillmRNA /(KmRNA^hillmRNA+Glucose^hillmRNA)-KdegmRNA*mRNA;
 DHXT= alphaHXT*mRNA- (maxdegHXT/(1+(Glucose/degThreshHXT)^hillDegHXT))*HXT;

  DMatureHXT= KmatHXT;
  
  
 dout= zeros(3,1);
 dout(1)=DmRNA;
 dout(2)=DHXT;
 dout(3)=DMatureHXT*HXT-(KBleach*MatureHXT);
 
 
end 
end