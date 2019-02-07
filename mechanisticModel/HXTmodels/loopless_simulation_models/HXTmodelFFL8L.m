
function modelh= HXTmodelFFL8L(params, mediaInput)
%The transcription rate is proportional to the rgt2 gate. 
%translatability of mRNA is (inversely) proportional to the rgt2 gate



modelh= @transcriptModelh;

function dout=transcriptModelh(t, x)

mediaInput=smooth(mediaInput);
Glucose=interp1(0:(numel(mediaInput)-1), mediaInput,t);

mRNA= x(1);
HXT= x(2);
matureHXT=x(3);



KRgt2= params(1);
KdegmRNA=params(2);
KdegHXT=params(3);
hillRgt2=params(4);
KmatHXT=params(5);



 Rgt2= Glucose^hillRgt2/(KRgt2+Glucose^hillRgt2);
 translatability= 1/ (1+ (Glucose/KRgt2)^hillRgt2);
 DmRNA= Rgt2- KdegmRNA*mRNA;
 DHXT= translatability*mRNA  -KdegHXT*HXT;
 DmatureHXT= KmatHXT*HXT- KdegHXT*matureHXT;
  dout= zeros(3,1);
 dout(1)=DmRNA;
 dout(2)=DHXT;
 dout(3)=DmatureHXT;


 
end 
end