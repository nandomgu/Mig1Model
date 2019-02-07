
function modelh= simpleModel5(params, args)
%%Model 1 is simple threshold-activated production of Hxt, including the
%%hill hypothesis.
%preferedOutVar=2
mediaInput=args.input;
mig1Input=args.mig1;

modelh= @transcriptModelh;
times=linspace(0,20,250);
times=times(1:numel(mediaInput));
mig1Input=mig1Input(1:numel(mediaInput));
function dout=transcriptModelh(t, x)

%mediaInput=smooth(mediaInput);
Glucose=interp1(times, mediaInput,t);

mig1= interp1(times, mig1Input,t);

HXT= x(1);

VHXT=params(1);
KHXT=params(2);
hillHXT=params(3);
KdegHXT=params(4);
KrepHXT=params(5);
hillrepHXT=params(6);

DHXT= (VHXT*Glucose/KHXT)^hillHXT / (1 + (Glucose/KHXT)^hillHXT+ (mig1/KrepHXT)^hillrepHXT)-KdegHXT*HXT;
 %(1+ (Glucose/KrepHXT)^hillrepHXT))
dout= zeros(1,1);
dout(1)=DHXT;

 
 
end 
end