
function modelh= HXTmodelFFL1L(params, mediaInput)
%%Model FFL1 is simple threshold-activated production of Hxt, including the
%%Hill Hypothesis, and threshold induced repression of Mig1.


modelh= @transcriptModelh;

function dout=transcriptModelh(t, x)

mediaInput=smooth(mediaInput);
Glucose=interp1(0:(numel(mediaInput)-1), mediaInput,t);

mRNA= x(1);
HXT= x(2);
Mig1=x(3);
Gene=x(4);


VmRNA=params(1);
KmRNA=params(2);
hillmRNA=params(3);
KdegmRNA=params(4);
alphaHXT=params(5);
KdegHXT= params(6);
alphaMig1=params(7);
KdegMig1=params(8);
KactiveMig1=params(9);
RepThreshMig1=params(10);
hillactiveMig1=params(11);
alphaGene=params(12);
KGene= params(13);
KdegGene= params(14);
hillGene=params(15);



 
 DMig1= alphaMig1-KdegMig1*Mig1;
 activeMig1= (Mig1*Glucose^hillactiveMig1)/(KactiveMig1^hillactiveMig1+ Glucose^hillactiveMig1);
 DmRNA= (VmRNA*(Glucose)^hillmRNA /(KmRNA^hillmRNA+Glucose^hillmRNA))/(1+(activeMig1/RepThreshMig1))+ -KdegmRNA*mRNA;
 DHXT= alphaHXT*mRNA- KdegHXT*HXT;
 DGene=(alphaGene*HXT^hillGene)/(KGene^hillGene+HXT^hillGene) -KdegGene*Gene;
 
 dout= zeros(4,1);
 dout(1)=DmRNA;
  dout(2)=DMig1;
 dout(3)=DHXT;
 dout(4)=DGene;

 
end 
end