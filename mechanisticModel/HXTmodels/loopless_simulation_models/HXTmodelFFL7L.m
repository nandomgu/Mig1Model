
function modelh= HXTmodelFFL7L(params, mediaInput)
%The transcription rate is proportional to the rgt2 gate. 
%production of the inhibitor is also proportional to rgt2
%translation rate is inversely proportional to proportional to the inhibitor.



modelh= @transcriptModelh;

function dout=transcriptModelh(t, x)

mediaInput=smooth(mediaInput);
Glucose=interp1(0:(numel(mediaInput)-1), mediaInput,t);

mRNA= x(1);
Inhibitor=x(2);
HXT= x(3);
matureHXT=x(4);





KRgt2= params(1);
KdegmRNA=params(2);
KdegInhibitor=params(3);
alphaHXT=params(4);
KdegHXT=params(5);
hillRgt2=params(6);
KmatHXT=params(7);



 Rgt2= Glucose^hillRgt2/(KRgt2+Glucose^hillRgt2);
 DmRNA= Rgt2- KdegmRNA*mRNA;
 DInhibitor= Rgt2- KdegInhibitor*Inhibitor;
 DHXT= alphaHXT*mRNA  - Inhibitor*HXT -KdegHXT*HXT;
 DmatureHXT= KmatHXT*HXT- KdegHXT*matureHXT;
  dout= zeros(3,1);
 dout(1)=DmRNA;
  dout(2)=DInhibitor;
 dout(3)=DHXT;
 dout(4)=DmatureHXT;


 
end 
end