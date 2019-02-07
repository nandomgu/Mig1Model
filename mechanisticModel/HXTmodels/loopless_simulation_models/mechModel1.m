
function modelh= mechModel1(params, args)

description='mechModel1- Mth1 is degraded glucose, Mig1 ';
mediaInput= args.input;
modelh= @transcriptModelh;
times=linspace(0,20,250);
times=times(1:numel(mediaInput));
function dout=transcriptModelh(t, x)

%mediaInput=smooth(mediaInput);
Glucose=interp1(times, mediaInput,t);

HXT= x(1);
MTH1=x(2);
MIG1=x(3);
HXTMTH1=x(4);
HXTMIG1=x(5);
MTH1MIG1=x(6);

VHXT=params(1);
Vloc=params(2);
basaldeg=params(3);
maxdegMT=params(4);
KdegMT=params(5);
hilldegMT=params(6);
basaldegMT=params(7);
deloc=params(8);
%repression of MTH1 by Mig1
KrepMT=params(9);
hillrepMT=params(10);
%Repression of HXT by Mig1
KrepHMG=params(11);
hillHMG=params(12);
%repression of HXT by MTH1
KrepHMT=params(13);
hillHMT=params(14);
%Mig1 localisation by glucose
KMG=params(15);
hillMG=params(16);


degMT=maxdegMT*(Glucose)^hilldegMT /(KdegMT^hilldegMT+Glucose^hilldegMT);

%fixing mth1 production to be 1
DMTH1= 1 / (1 + (MIG1/KrepMT)^hillrepMT) - (degMT+ basaldegMT)*MTH1 ;

%the equation of Mth1 when there is no mig1 to repress it
DMTH1MIG1= 1 / (1 + (0/KrepMT)^hillrepMT) - (degMT+ basaldegMT)*MTH1 ;



cytMIG1=1-MIG1;
DMIG1= Vloc*cytMIG1*(Glucose/KMG)^hillMG/(1 + (Glucose/KMG)^hillMG) -deloc*MIG1;

%start with basal degradation only
%deg=VdegHXT/(1+ (Glucose/threshdegHXT)^hilldeg)+basaldeg;

DHXT=  VHXT/ (1+ (MIG1/KrepHMG)^hillHMG +(MTH1/KrepHMT)^hillHMT) - basaldeg*HXT;
DHXTMTH1=  VHXT/ (1+ (MIG1/KrepHMG)^hillHMG +(0/KrepHMT)^hillHMT) - basaldeg*HXT;
DHXTMIG1=  VHXT/ (1+ (0/KrepHMG)^hillHMG +(MTH1/KrepHMT)^hillHMT) - basaldeg*HXT;


%(VHXT*Glucose/KHXT)^hillHXT / (1 + (Glucose/KHXT)^hillHXT+ (Glucose/KrepHXT)^hillrepHXT)-basaldeg*HXT;
 %(1+ (Glucose/KrepHXT)^hillrepHXT))
dout= zeros(3,1);
dout(1)=DHXT;
dout(2)=DMTH1;
dout(3)=DMIG1;
dout(4)=DHXTMTH1;
dout(5)=DHXTMIG1;
dout(6)=DMTH1MIG1;
 
end 
end