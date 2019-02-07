
function modelh= HXTmodel608L(params, mediaInput)
%model with 2 transcription hill functions, this time one activated by
%glucose and another repressed by glucose
%model 605 contains 2 different thresholds of glucose
modelh= @transcriptModelh;

function dout=transcriptModelh(t, x)

mediaInput=smooth(mediaInput);
Glucose=interp1(0:(numel(mediaInput)-1), mediaInput,t);


mRNA= x(1);
HXT= x(2);
MatureHXT= x(3);
rep1= 1;
rep2= 1;

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
V2mRNA=params(13);
kdegrep1=params(14);
kdegrep2=params(15);
Vrep1=.1;
Vrep2=.1;
degThreshBelow=params(16);


Drep1= Vrep1 -kdegrep1*Glucose*rep1;
Drep2= Vrep2 -kdegrep2*Glucose*rep2;


DmRNA=VmRNA/(1+(rep1^hillmRNA/KmRNA^hillmRNA)+(rep2^hill2mRNA/K2mRNA^hill2mRNA))- (KdegmRNA*(Glucose)^hillDegHXT /(degThreshHXT^hillDegHXT+Glucose^hillDegHXT))*mRNA;
degHXT= maxdegHXT*(Glucose)^hillDegHXT /(degThreshHXT^hillDegHXT+Glucose^hillDegHXT);
degBelowHXT= VmRNA/(1+(Glucose^hillDegHXT/degThreshBelow^hillDegHXT));
 
 
  
 DHXT= (alphaHXT/(1+(Glucose^hillDegHXT/degThreshHXT^hillDegHXT)))*mRNA -(degHXT +degBelowHXT)*HXT ;
 

 DMatureHXT= KmatHXT;
 dout= zeros(3,1);
 dout(1)=DmRNA;
 dout(2)=DHXT;
 dout(3)=DMatureHXT*HXT-(KBleach*MatureHXT)-degBelowHXT -degHXT*MatureHXT;
 
 
end 
end