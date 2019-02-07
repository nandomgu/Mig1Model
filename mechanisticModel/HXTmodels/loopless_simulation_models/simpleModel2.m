
function modelh= simpleModel2(params, args)
%%Model 1 is simple threshold-activated production of Hxt, including the
%%hill hypothesis.
%preferedOutVar=2
mediaInput=args.input;
modelh= @transcriptModelh;
times=args.times;
function dout=transcriptModelh(t, x)

%mediaInput=smooth(mediaInput);
Glucose=interp1(times, mediaInput,t);

HXT= x(1);

VHXT=params(1);
KHXT=params(2);
hillHXT=params(3);
VdegHXT=params(4);
Kdeg=params(5);
hilldeg=params(6);
basaldeg=params(7);


deg=VdegHXT*(Glucose)^hilldeg /(Kdeg^hilldeg+Glucose^hilldeg)+basaldeg;
DHXT= (VHXT*(Glucose)^hillHXT /(KHXT^hillHXT+Glucose^hillHXT)) -deg*HXT;
 %(1+ (Glucose/KrepHXT)^hillrepHXT))
dout= zeros(1,1);
dout(1)=DHXT;

 
 
end 
end