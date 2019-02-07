
function modelh= HXTmodel110L(params, mediaInput)
%model110. model 9 had threshold-induced degradation of protein above a
%threshold. This one also incorporates threshold induced inhibition of protein degradation below a threshold. 

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
degThreshHXTbelow= params(11);
hillDegHXTbelow= params(12);


DmRNA= VmRNA*((Glucose^hillmRNA) /((KmRNA^hillmRNA)+(Glucose^hillmRNA)))- (KdegmRNA*mRNA);

DHXT= alphaHXT*mRNA- ((maxdegHXT^hillDegHXT)/(degThreshHXT^hillDegHXT+Glucose^hillDegHXT))-((maxdegHXT^hillDegHXTbelow)/(degThreshHXTbelow^hillDegHXTbelow+Glucose^hillDegHXTbelow))*HXT;

  DMatureHXT= KmatHXT;
 dout= zeros(3,1);
 dout(1)=DmRNA;
 dout(2)=DHXT;
 dout(3)=DMatureHXT*HXT-(KBleach*MatureHXT);
 
 
end 
end