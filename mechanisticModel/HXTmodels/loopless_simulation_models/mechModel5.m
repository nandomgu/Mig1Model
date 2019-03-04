
function modelh= mechModel5(params, args)

description='mechModel3- Mth1 is degraded glucose, glucose -> nuclear mig1; Mth1 is NOT repressed by Mig1Std1 is degraded by glucoseStd1 induces delocalisation of Mig1 mig1, mth1 repress th hxt, hxt is degraded in the absence of glucose ';
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
STD1=x(7);
HXTP2=x(8);
HXTP4=x(9);

VHXT=params(1);
Vloc=params(2);
basaldeg=params(3);
maxdegMT=params(4);
%gluc induced deg threshold of Std1
KdegMT=params(5);
hilldegMT=params(6);
basaldegMT=params(7);
maxdegST=params(8);
%gluc induced deg threshold of Std1
KdegST=params(9);
hilldegST=params(10);
basaldegST=params(11);
%with KinhSMG, std1 inhibits the action of mig1.
%first we try a proportional one. 
KinhSMG=params(12);
deloc=params(13);
%Repression of HXT by Mig1
KrepHMG=params(14);
hillHMG=params(15);
%repression of HXT by MTH1
KrepHMT=params(16);
hillHMT=params(17);
%Mig1 localisation by glucose
KMG=params(18);
hillMG=params(19);
VdegHXT=params(20);
threshdegHXT=params(21);
hilldegHXT=params(22);

degMT=maxdegMT*(Glucose)^hilldegMT /(KdegMT^hilldegMT+Glucose^hilldegMT);
degST=maxdegST*(Glucose)^hilldegST /(KdegST^hilldegST+Glucose^hilldegST);



%fixing mth1 production to be 1
DMTH1= 1  - (degMT+ basaldegMT)*MTH1 ;

%fixing std1 production to be 1
DSTD1= 1 - (degST+ basaldegST)*STD1 ;



%the equation of Mth1 when there is no mig1 to repress it
DMTH1MIG1= 1  - (degMT+ basaldegMT)*MTH1 ;



cytMIG1=1-MIG1;
DMIG1= Vloc*cytMIG1*(Glucose/KMG)^hillMG/(1 + (Glucose/KMG)^hillMG) -deloc*MIG1 -KinhSMG*STD1*MIG1;

%start with basal degradation only
deg=VdegHXT/(1+ (Glucose/threshdegHXT)^hilldegHXT)+basaldeg;

DHXT=  VHXT/ (1+ (MIG1/KrepHMG)^hillHMG +(MTH1/KrepHMT)^hillHMT) - deg*HXT;
DHXTMTH1=  VHXT/ (1+ (MIG1/KrepHMG)^hillHMG +(0/KrepHMT)^hillHMT) - deg*HXT;
DHXTMIG1=  VHXT/ (1+ (0/KrepHMG)^hillHMG +(MTH1/KrepHMT)^hillHMT) - deg*HXT;


%(VHXT*Glucose/KHXT)^hillHXT / (1 + (Glucose/KHXT)^hillHXT+ (Glucose/KrepHXT)^hillrepHXT)-basaldeg*HXT;
 %(1+ (Glucose/KrepHXT)^hillrepHXT))
dout= zeros(3,1);
dout(1)=DHXT;
dout(2)=DMTH1;
dout(3)=DMIG1;
dout(4)=DHXTMTH1;
dout(5)=DHXTMIG1;
dout(6)=DMTH1MIG1;
dout(7)=DSTD1 ;
dout(8)=DHXTP2 ;
dout(9)=DHXTP4 ;
end 
end