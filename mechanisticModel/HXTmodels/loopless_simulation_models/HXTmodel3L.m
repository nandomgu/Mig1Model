
function modelh= HXTmodel3L(params, mediaInput)
%model3. model 2 plus threshold-induced degradation of mRNA
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
MaxdegmRNA=params(4);
alphaHXT=params(5);
KdegHXT= params(6);
KmatHXT=params(7);
KBleach=params(8);
degThreshmRNA=params(9);
hillDegmRNA=params(10);

  
 DmRNA= VmRNA*(Glucose)^hillmRNA /(KmRNA^hillmRNA+Glucose^hillmRNA)- (MaxdegmRNA/(1+(Glucose/degThreshmRNA)^hillDegmRNA)*mRNA);
 DHXT= alphaHXT*mRNA- KdegHXT*HXT;

  DMatureHXT= KmatHXT;
 dout= zeros(3,1);
 dout(1)=DmRNA;
 dout(2)=DHXT;
 dout(3)=DMatureHXT*HXT-(KBleach*MatureHXT);
 
 
end 
end