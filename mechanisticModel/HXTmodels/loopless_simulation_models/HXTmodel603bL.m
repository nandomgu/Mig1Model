function modelh= HXTmodel603L(params, mediaInput)
%%model 603 Loopless.
%%a) there are 2 transcription rates. They are implemented through
%%parammeter v1, which represents one rate as a fraction of the other.
%%This way it is intended to reflect a basal expression rate, never 0.

modelh= @transcriptModelh;

function dout=transcriptModelh(t, x)

mediaInput=smooth(mediaInput);
Glucose=interp1(0:(numel(mediaInput)-1), mediaInput,t);


mRNA= x(1);
HXT= x(2);
MatureHXT= x(3);

VmRNA=params(1);
KmRNA=params(2);
hillmRNA=params(3);
KdegmRNA=params(4);
alphaHXT=params(5);
maxdegHXT= params(6);
KmatHXT=params(7);
KBleach=params(8);
degThreshHXT=params(9);
hillDegHXT=params(10);
K2mRNA=params(11);
hill2mRNA=params(12);


 DmRNA= ((VmRNA*(Glucose)^hillmRNA /(KmRNA^hillmRNA+Glucose^hillmRNA))/(1+(Glucose^hill2mRNA/K2mRNA^hill2mRNA)))- (KdegmRNA*(Glucose)^hillDegHXT /(degThreshHXT^hillDegHXT+Glucose^hillDegHXT))*mRNA;
 
DHXT= alphaHXT*mRNA-  ((maxdegHXT/(1+(Glucose/degThreshHXT)^hillDegHXT))*HXT);

 DMatureHXT= KmatHXT;
 dout= zeros(3,1);
 dout(1)=DmRNA;
 dout(2)=DHXT;
 dout(3)=DMatureHXT*HXT-(KBleach*MatureHXT);
 
 
end 
end