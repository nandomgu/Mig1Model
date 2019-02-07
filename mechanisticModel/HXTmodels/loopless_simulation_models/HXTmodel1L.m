
function modelh= HXTmodel1L(params, mediaInput)
%%Model 1 is simple threshold-activated production of Hxt, including the
%%hill hypothesis.
%preferedOutVar=2

modelh= @transcriptModelh;

function dout=transcriptModelh(t, x)

mediaInput=smooth(mediaInput);
Glucose=interp1(0:(numel(mediaInput)-1), mediaInput,t);

mRNA= x(1);
HXT= x(2);

VmRNA=params(1);
KmRNA=params(2);
hillmRNA=params(3);
KdegmRNA=params(4);
alphaHXT=params(5);
KdegHXT= params(6);






 DmRNA= VmRNA*(Glucose)^hillmRNA /(KmRNA^hillmRNA+Glucose^hillmRNA)-KdegmRNA*mRNA;
 DHXT= alphaHXT*mRNA- KdegHXT*HXT;
 
 dout= zeros(2,1);
 dout(1)=DmRNA;
 dout(2)=DHXT;
 
 
end 
end