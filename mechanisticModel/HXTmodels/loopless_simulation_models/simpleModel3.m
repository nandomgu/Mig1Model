
function modelh= simpleModel3(params, args)
%%Model 3 is glucose repressed degradatio of Hxt
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


deg=VdegHXT/(1+ (Glucose/Kdeg)^hilldeg)+basaldeg;
DHXT= (VHXT*(Glucose)^hillHXT /(KHXT^hillHXT+Glucose^hillHXT)) -deg*HXT;
 %(1+ (Glucose/KrepHXT)^hillrepHXT))
dout= zeros(1,1);
dout(1)=DHXT;

 
 
end 
end