
function modelh= simpleModel6(params, args)

description='simpleModel6 is a threshold induced glucose induction and repression model but wih a hill function on induction too, similar to model 1';
mediaInput= args.input;
modelh= @transcriptModelh;
times=linspace(0,20,250);
times=times(1:numel(mediaInput));
function dout=transcriptModelh(t, x)

%mediaInput=smooth(mediaInput);
Glucose=interp1(times, mediaInput,t);

HXT= x(1);

VHXT=params(1);
KHXT=params(2);
hillHXT=params(3);
KdegHXT=params(4);
KrepHXT=params(5);
hillrepHXT=params(6);

DHXT=  (VHXT*Glucose/KHXT)^hillHXT / (1 + (Glucose/KHXT)^hillHXT+ (Glucose/KrepHXT)^hillrepHXT)-KdegHXT*HXT;
 %(1+ (Glucose/KrepHXT)^hillrepHXT))
dout= zeros(1,1);
dout(1)=DHXT;

 
 
end 
end