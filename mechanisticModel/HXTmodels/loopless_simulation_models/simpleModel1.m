
function modelh= simpleModel1(params, args)
%%Model 1 is simple threshold-activated production of Hxt, including the
%%hill hypothesis.
%preferedOutVar=2
mediaInput=args.input;
times=args.times;
modelh= @transcriptModelh;

function dout=transcriptModelh(t, x)

%mediaInput=smooth(mediaInput);
Glucose=interp1(times, mediaInput,t);

HXT= x(1);

VHXT=params(1);
KHXT=params(2);
hillHXT=params(3);
KdegHXT=params(4);


DHXT= (VHXT*(Glucose)^hillHXT /(KHXT^hillHXT+Glucose^hillHXT)) -KdegHXT*HXT;
 %(1+ (Glucose/KrepHXT)^hillrepHXT))
dout= zeros(1,1);
dout(1)=DHXT;

 
 
end 
end