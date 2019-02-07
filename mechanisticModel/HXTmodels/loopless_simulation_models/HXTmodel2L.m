
function modelh= HXTmodel2L(params, mediaInput)
%%model2. model 1 adding the GFP maturation reaction 
%preferedOutVar=3

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
KdegHXT= params(6);
KmatHXT=params(7);
KBleach=params(8);




 DmRNA= VmRNA*(Glucose)^hillmRNA /(KmRNA^hillmRNA+Glucose^hillmRNA)-KdegmRNA*mRNA;
 DHXT= alphaHXT*mRNA- KdegHXT*HXT;
 DMatureHXT= KmatHXT;
 
 dout= zeros(3,1);
 dout(1)=DmRNA;
 dout(2)=DHXT;
 dout(3)=DMatureHXT*HXT-(KBleach*MatureHXT);
 
 
end 
end