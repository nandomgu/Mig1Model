
function modelh= HXTmodelFFL4L(params, mediaInput)
%%Glucose activates an Rgt2 signaling gate. Glucose also activates Std1
%%(and possibly a repression mechanism).
%Rgt2 degrades Std1
%Std1 activates Hxt4


modelh= @transcriptModelh;

function dout=transcriptModelh(t, x)

mediaInput=smooth(mediaInput);
Glucose=interp1(0:(numel(mediaInput)-1), mediaInput,t);

mRNA= x(1);
HXT= x(2);
Mig1=x(3);
Std1=x(4);


VmRNA=params(1);
%KmRNA=params(2);
hillmRNA=params(2);
KdegmRNA=params(3);
alphaHXT=params(4);
KdegHXT= params(5);
alphaMig1=params(6);
KdegMig1=params(7);
KactiveMig1=params(8);
RepThreshMig1=params(9);
alphaStd1=params(10);
maxdegStd1=params(11);
KactivateStd1=params(12);
RepThreshMth1=params(13);
hilldegSth1=params(14);
KRgt2= params(15);
alphaRgt2=params(16);

%Model equations
 %DmRNA= vmaxRNA*(input)^n /(KmRNA^n+input^n)-KdegmRNA*mRNA
 %DProtein= alphaHXT*mRNA- KdegProt*Protein
 
 Rgt2= Glucose/(KRgt2+Glucose);
 DMig1= alphaMig1*Rgt2/(KactiveMig1+Rgt2) -KdegMig1*Mig1;
 %Drepressor= Rgt2;
 %activeMig1= (Mig1*Glucose^4)/(KactiveMig1+ Glucose^4);
 activeMig1=Mig1;
 DmRNA= (VmRNA*Std1/ (KactivateStd1+Std1))/(1+(activeMig1/RepThreshMig1) )-KdegmRNA*mRNA;  
 
 DHXT= alphaHXT*mRNA- KdegHXT*HXT;
 
 DStd1= alphaStd1*Glucose/(1+Glucose) -  Rgt2*Std1;
%(K-(maxdegStd1/(1+(Glucose^hilldegStd1/degThreshStd1^hilldegStd1))))
 dout= zeros(3,1);
 %dout(1)= DGlucose;
 dout(1)=DmRNA;
  dout(2)=DMig1;
 dout(3)=DHXT;
 dout(4)=DStd1;

 
end 
end