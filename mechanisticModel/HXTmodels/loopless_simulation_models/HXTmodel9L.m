
function modelh= HXTmodel9L(params, mediaInput)
%model3. model 2 plus threshold-induced degradation of protein. the
%difference from model 6 is that max degradation is induced ABOVE a glucose
%threshold rather than below.

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


%Model equations
 %DmRNA= vmaxRNA*(input)^n /(KmRNA^n+input^n)-KdegmRNA*mRNA
 %DProtein= alphaHXT*mRNA- KdegProt*Protein
 %DGlucose=0;
  
DmRNA= VmRNA*((Glucose^hillmRNA) /((KmRNA^hillmRNA)+(Glucose^hillmRNA)))- (KdegmRNA*mRNA);
%DmRNA=0; 
DHXT= alphaHXT*mRNA- ((maxdegHXT/(degThreshHXT^hillDegHXT+Glucose^hillDegHXT)))*HXT;

  DMatureHXT= KmatHXT;
 dout= zeros(3,1);
 %dout(1)= DGlucose;
 dout(1)=DmRNA;
 dout(2)=DHXT;
 dout(3)=DMatureHXT*HXT-(KBleach*MatureHXT);
 
 
end 
end