
function modelh= simpleModel8(params, args)

description='simpleModel8- Mth1 is degraded glucose';
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
basaldeg=params(4);
KrepHXT=params(5);
hillrepHXT=params(6);
VdegHXT=params(7);
threshdegHXT=params(8);
hilldeg=params(9);

degMT=maxdegMT*(Glucose)^hilldegMT /(KdegMT^hilldegMT+Glucose^hilldegMT);

DMTH1= 1 / (1 + (Glucose/KrepMT)^hillrepMT)  -   (degMT+ basaldegMT)*MTH1 




deg=VdegHXT/(1+ (Glucose/threshdegHXT)^hilldeg)+basaldeg;



DHXT=  (VHXT*Glucose/KHXT)^hillHXT / (1 + (Glucose/KHXT)^hillHXT+ (Glucose/KrepHXT)^hillrepHXT)-deg*HXT;
 %(1+ (Glucose/KrepHXT)^hillrepHXT))
dout= zeros(1,1);
dout(1)=DHXT;

 
 
end 
end