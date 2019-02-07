
function modelh= HXTmodel602L(params, mediaInput)
%model3. model 2 plus threshold-induced degradation of protein

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
V2mRNA=params(13)/VmRNA;
  
 DmRNA= VmRNA*(Glucose)^hillmRNA /(KmRNA^hillmRNA+Glucose^hillmRNA)+ (VmRNA/(1+(Glucose^hillmRNA/K2mRNA^hill2mRNA)))-KdegmRNA*mRNA;
 DHXT= alphaHXT*mRNA- ((maxdegHXT/(1+(Glucose/degThreshHXT)^hillDegHXT))*HXT) - maxdegHXT*(Glucose)^hillmRNA /(KmRNA^hillmRNA+Glucose^hillmRNA)*HXT;

 DMatureHXT= KmatHXT;
 dout= zeros(3,1);
 dout(1)=DmRNA;
 dout(2)=DHXT;
 dout(3)=DMatureHXT*HXT-(KBleach*MatureHXT);
 
 
end 
end