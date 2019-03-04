
function modelh= mechModel10(params, args)

description=['mechModel10- glucose deg mth1'  ...
             'glucose deg std1'...
             'mth1 rep Std1' ...
             'std1 inhibits mig1' ...
             'Mig1, Mth1, AND Std1 repress hxt'...
             'Mig2 is controled by Mth1'...
             'Mig2 represses Mth1. also substituting that for Mig1 here. '...
             'Mig2 leaves the nucleus by means independent of Std1'
             ];
             
             
mediaInput= args.input;
modelh= @transcriptModelh;
times=linspace(0,20,250);
times=times(1:numel(mediaInput));
function dout=transcriptModelh(t, x)

%mediaInput=smooth(mediaInput);
Glucose=interp1(times, mediaInput,t);

% if isfield(params, 'mutant') && ~isempty(params.mutant)  && params.mutant~=0 %if we'd like the simulation to have a mutant
%     
%     x(params.mutant)=0;
% end

HXT= x(1);
MTH1=x(2);
MIG1=x(3);
HXTMTH1=x(4);
HXTMIG1=x(5);
MTH1MIG1=x(6);
STD1=x(7);
MIG2=x(8);



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
VinactGM2=params(35);
KinactGM2=params(36);
hillinactGM2=params(37);


 %if we'd like the simulation to have a mutant

degMT=maxdegMT*(Glucose)^hilldegMT /(KdegMT^hilldegMT+Glucose^hilldegMT);
degST=maxdegST*(Glucose)^hilldegST /(KdegST^hilldegST+Glucose^hilldegST);

if isfield(params, 'mutant') && ~isempty(params.mutant)  && params.mutant==8
    DMIG2=0;
else
DMIG2= 1/(1+ (MTH1/KrepSTM2)^hillrepSTM2) -KinhSM2*STD1*MIG2 -basaldegM2*MIG2   -  (VinactGM2 / (1 + (Glucose/KinactGM2)^hillinactGM2))*MIG2;
end
if isfield(params, 'mutant') && ~isempty(params.mutant)  && params.mutant==2
    DMTH1=0;
else
%fixing mth1 production to be 1
DMTH1= 1/ (1 + (MIG2/KrepM2MT)^hillrepM2MT)  - (degMT+ basaldegMT)*MTH1 ;
end
%fixing std1 production to be 1
if isfield(params, 'mutant') && ~isempty(params.mutant)  && params.mutant==7
    DSTD1=0;
else
DSTD1= 1/(1+(MTH1/KrepMTST)^hillrepMTST) - (degST+ basaldegST)*STD1 ;
end


%the equation of Mth1 when there is no mig1 to repress it
DMTH1MIG1= 1  - (degMT+ basaldegMT)*MTH1 ;


if isfield(params, 'mutant') && ~isempty(params.mutant)  && params.mutant==3
    DMIG1=0;
else
cytMIG1=1-MIG1;
DMIG1= Vloc*cytMIG1*(Glucose/KMG)^hillMG/(1 + (Glucose/KMG)^hillMG) -deloc*MIG1 -KinhSMG*STD1*MIG1;
end
%start with basal degradation only
deg=VdegHXT/(1+ (Glucose/threshdegHXT)^hilldegHXT)+basaldeg;

DHXT=  VHXT/ (1+ (MIG1/KrepHMG)^hillHMG +(MTH1/KrepHMT)^hillHMT +(MIG2/KrepHM2)^hillrepHM2 +(STD1/KrepHST)^hillrepHST) - deg*HXT;
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
dout(8)=DMIG2 ;
end 
end