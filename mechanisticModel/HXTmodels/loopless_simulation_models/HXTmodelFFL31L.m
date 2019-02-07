function modelh= HXTmodelFFL31L(params, mediaInput)
%%an incoherent type 3 ffl , where glucose directly inhibits the gene, and also activates
%%B which activates the gene.


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
KactiveB=params(7);
KB=params(8);
hillactiveB=params(9);




 activeB= (Glucose^hillactiveB)/(KactiveB^hillactiveB+ Glucose^hillactiveB);
 DmRNA= (VmRNA*activeB /(KB+activeB))/(1+((Glucose^hillmRNA)/(KmRNA^hillmRNA)))+ -KdegmRNA*mRNA;
 DHXT= alphaHXT*mRNA- KdegHXT*HXT;
 
 
 dout= zeros(2,1);
 dout(1)=DmRNA;
 dout(2)=DHXT;

 
end 
end