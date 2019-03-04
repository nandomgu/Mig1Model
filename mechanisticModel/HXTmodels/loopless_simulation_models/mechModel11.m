
function modelh= mechModel11(params, args)

description=['mechModel9- glucose deg mth1'  ...
             'glucose deg std1'...
             'mth1 rep Std1' ...
             'std1 inhibits mig1' ...
             'Mig1, Mth1, AND Std1 repress hxt'...
             'Mig2 is controled by Mth1'
             'Mig2 represses Mth1. also substituting that for Mig1 here. '
             'This model sumulates also 3 glucose concentrations simultaneously'
             ];
             
             
mediaInput= args.input;
modelh= @transcriptModelh;
times=linspace(0,20,250);
times=times(1:numel(mediaInput));
function dout=transcriptModelh(t, x)

%mediaInput=smooth(mediaInput);
Glucose=interp1(times, mediaInput,t);
Glucose2=Glucose *.2;
Glucose4=Glucose*.4;


HXT= x(1);
MTH1=x(2);
MIG1=x(3);
HXTMTH1= x(4);
HXTMIG1= x(5);
MTH1MIG1=x(6);
STD1=x(7);
MIG2=x(8);


HXT2= x(9);
MTH12=x(10);
MIG12=x(11);
HXTMTH12= x(12);
HXTMIG12= x(13);
MTH1MIG12=x(14);
STD12=x(15);
MIG22=x(16);

HXT4= x(17);
MTH14=x(18);
MIG14=x(19);
HXTMTH14= x(20);
HXTMIG14= x(21);
MTH1MIG14=x(22);
STD14=x(23);
MIG24=x(24);





VHXT=params(1);
Vloc=params(2);
basaldeg=params(3);
maxdegMT=params(4);
%gluc induced deg threshold of Mth1
KdegMT=params(5);
hilldegMT=params(6);
basaldegMT=params(7);
maxdegST=params(8);
%mth1 repression of Std1
KrepMTST=params(9);
hillrepMTST=params(10);
%gluc induced deg threshold of Std1
KdegST=params(11);
hilldegST=params(12);
basaldegST=params(13);
%with KinhSMG, std1 inhibits the action of mig1.
%first we try a proportional one. 
KinhSMG=params(14);
deloc=params(15);
%Repression of HXT by Mig1
KrepHMG=params(16);
hillHMG=params(17);
%repression of HXT by MTH1
KrepHMT=params(18);
hillHMT=params(19);
%repression of HXT by STD1
KrepHST=params(20);
hillrepHST=params(21);
%Mig1 localisation by glucose
KMG=params(22);
hillMG=params(23);
VdegHXT=params(24);
threshdegHXT=params(25);
hilldegHXT=params(26);
basaldegM2=params(27);
KrepSTM2=params(28);
hillrepSTM2=params(29);
KrepHM2=params(30);
hillrepHM2=params(31);
KinhSM2=params(32);
KrepM2MT=params(33);
hillrepM2MT=params(34);

degMT=maxdegMT*(Glucose)^hilldegMT /(KdegMT^hilldegMT+Glucose^hilldegMT);
degST=maxdegST*(Glucose)^hilldegST /(KdegST^hilldegST+Glucose^hilldegST);
DMIG2= 1/(1+ (MTH1/KrepSTM2)^hillrepSTM2) -KinhSM2*STD1*MIG2 -basaldegM2*MIG2;
%fixing mth1 production to be 1
DMTH1= 1/ (1 + (MIG2/KrepM2MT)^hillrepM2MT)  - (degMT+ basaldegMT)*MTH1 ;
%fixing std1 production to be 1
DSTD1= 1/(1+(MTH1/KrepMTST)^hillrepMTST) - (degST+ basaldegST)*STD1 ;
%the equation of Mth1 when there is no mig1 to repress it
DMTH1MIG1= 1  - (degMT+ basaldegMT)*MTH1 ;
cytMIG1=1-MIG1;
DMIG1= Vloc*cytMIG1*(Glucose/KMG)^hillMG/(1 + (Glucose/KMG)^hillMG) -deloc*MIG1 -KinhSMG*STD1*MIG1;
%start with basal degradation only
deg=VdegHXT/(1+ (Glucose/threshdegHXT)^hilldegHXT)+basaldeg;
DHXT=  VHXT/ (1+ (MIG1/KrepHMG)^hillHMG +(MTH1/KrepHMT)^hillHMT +(MIG2/KrepHM2)^hillrepHM2 +(STD1/KrepHST)^hillrepHST) - deg*HXT;
DHXTMTH1=  VHXT/ (1+ (MIG1/KrepHMG)^hillHMG +(0/KrepHMT)^hillHMT) - deg*HXT;
DHXTMIG1=  VHXT/ (1+ (0/KrepHMG)^hillHMG +(MTH1/KrepHMT)^hillHMT) - deg*HXT;


degMT2=maxdegMT*(Glucose2)^hilldegMT /(KdegMT^hilldegMT+Glucose2^hilldegMT);
degST2=maxdegST*(Glucose2)^hilldegST /(KdegST^hilldegST+Glucose2^hilldegST);
DMIG22= 1/(1+ (MTH12/KrepSTM2)^hillrepSTM2) -KinhSM2*STD12*MIG22 -basaldegM2*MIG22;
%fixing mth1 production to be 1
DMTH12= 1/ (1 + (MIG22/KrepM2MT)^hillrepM2MT)  - (degMT2+ basaldegMT)*MTH12 ;
%fixing std1 production to be 1
DSTD12= 1/(1+(MTH12/KrepMTST)^hillrepMTST) - (degST2+ basaldegST)*STD1 ;
%the equation of Mth1 when there is no mig1 to repress it
DMTH1MIG12= 1  - (degMT+ basaldegMT)*MTH1 ;
cytMIG12=1-MIG12;
DMIG12= Vloc*cytMIG12*(Glucose2/KMG)^hillMG/(1 + (Glucose2/KMG)^hillMG) -deloc*MIG12 -KinhSMG*STD12*MIG12;
%start with basal degradation only
deg2=VdegHXT/(1+ (Glucose2/threshdegHXT)^hilldegHXT)+basaldeg;
DHXT2=  VHXT2/ (1+ (MIG12/KrepHMG)^hillHMG +(MTH12/KrepHMT)^hillHMT +(MIG22/KrepHM2)^hillrepHM2 +(STD12/KrepHST)^hillrepHST) - deg*HXT;
DHXTMTH12=  VHXT/ (1+ (MIG12/KrepHMG)^hillHMG +(0/KrepHMT)^hillHMT) - deg*HXT2;
DHXTMIG12=  VHXT/ (1+ (0/KrepHMG)^hillHMG +(MTH12/KrepHMT)^hillHMT) - deg*HXT2;


degMT4=maxdegMT*(Glucose4)^hilldegMT /(KdegMT^hilldegMT+Glucose4^hilldegMT);
degST4=maxdegST*(Glucose4)^hilldegST /(KdegST^hilldegST+Glucose4^hilldegST);
DMIG24= 1/(1+ (MTH14/KrepSTM2)^hillrepSTM2) -KinhSM2*STD14*MIG24 -basaldegM2*MIG24;
%fixing mth1 production to be 1
DMTH14= 1/ (1 + (MIG24/KrepM2MT)^hillrepM2MT)  - (degMT4+ basaldegMT)*MTH14 ;
%fixing std1 production to be 1
DSTD14= 1/(1+(MTH14/KrepMTST)^hillrepMTST) - (degST4+ basaldegST)*STD14 ;
%the equation of Mth1 when there is no mig1 to repress it
DMTH1MIG14= 1  - (degMT4+ basaldegMT)*MTH14 ;
cytMIG14=1-MIG14;
DMIG14= Vloc*cytMIG14*(Glucose4/KMG)^hillMG/(1 + (Glucose4/KMG)^hillMG) -deloc*MIG14 -KinhSMG*STD14*MIG14;
%start with basal degradation only
deg=VdegHXT/(1+ (Glucose4/threshdegHXT)^hilldegHXT)+basaldeg;
DHXT4=  VHXT/ (1+ (MIG14/KrepHMG)^hillHMG +(MTH14/KrepHMT)^hillHMT +(MIG24/KrepHM2)^hillrepHM2 +(STD14/KrepHST)^hillrepHST) - deg*HXT4;
DHXTMTH14=  VHXT/ (1+ (MIG14/KrepHMG)^hillHMG +(0/KrepHMT)^hillHMT) - deg*HXT4;
DHXTMIG14=  VHXT/ (1+ (0/KrepHMG)^hillHMG +(MTH14/KrepHMT)^hillHMT) - deg*HXT4;




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
dout(8)=DMIG2 ;
end 
end