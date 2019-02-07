function modelh= HXTmodel703L(params, mediaInput)
%%model 703l. model 603 l but with an mRNA 'accumulator'. that is, above
%%certain level of glucose, mRNA is accumulated, and below it it is
%%released.

modelh= @transcriptModelh;

function dout=transcriptModelh(t, x)

mediaInput=smooth(mediaInput);
Glucose=interp1(0:(numel(mediaInput)-1), mediaInput,t);


mRNA= x(1);
HXT= x(2);
MatureHXT= x(3);
V2mRNA=1;
accumRNA=x(4);

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
kAccumRNA=params(13);
kReleasemRNA=params(14);


DaccumRNA= (kAccumRNA*(Glucose/KmRNA+Glucose) * mRNA )- kReleasemRNA*(accumRNA); %rna accumulation depends on the
DmRNA=  kReleasemRNA*(accumRNA)+ ((VmRNA*(Glucose)^hillmRNA /(KmRNA^hillmRNA+Glucose^hillmRNA))/(1+(Glucose^hill2mRNA/K2mRNA^hill2mRNA)))- (KdegmRNA*(Glucose)^hillDegHXT /(degThreshHXT^hillDegHXT+Glucose^hillDegHXT))*mRNA -kAccumRNA;
 
 


DHXT= alphaHXT*mRNA-  ((maxdegHXT/(1+(Glucose/degThreshHXT)^hillDegHXT))*HXT);

 DMatureHXT= KmatHXT;
 dout= zeros(3,1);
 dout(1)=DmRNA;
 dout(2)=DHXT;
 dout(3)=DMatureHXT*HXT-(KBleach*MatureHXT);
 
 
end 
end