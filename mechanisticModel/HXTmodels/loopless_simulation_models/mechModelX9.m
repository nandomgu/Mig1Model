
function modelh= mechModelX9(params, args)

description=['simulation of Hxt in all mutants and 3 glucose concentrations (0.2 0.4 1)'...
              'mth1 is not repressed at all'...
              'std1 does not inhibit Mig2'...
              'make sure to simulate without params 33 and 34'
];
             
             
mediaInput= args.input;
modelh= @transcriptModelh;
times=linspace(0,20,250);
times=times(1:numel(mediaInput));


function dout=transcriptModelh(t, x)

%mediaInput=smooth(mediaInput);
Glucose1=interp1(times, mediaInput,t);
% if isfield(params, 'mutant') && ~isempty(params.mutant)  && params.mutant~=0 %if we'd like the simulation to have a mutant
%     
%     x(params.mutant)=0;
% end
regulators=[];
HXT_2= x(1);
HXT_4= x(2);
HXT_1= x(3);
HXTMIG1_2=x(4);
HXTMIG1_4=x(5);
HXTMIG1_1=x(6);
HXTMTH1_2=x(7);
HXTMTH1_4=x(8);
HXTMTH1_1=x(9);
HXTSTD1_2=x(10);
HXTSTD1_4=x(11);
HXTSTD1_1=x(12);
HXTRGT2_2=x(13);
HXTRGT2_4=x(14);
HXTRGT2_1=x(15);
HXTSNF3_2=x(16);
HXTSNF3_4=x(17);
HXTSNF3_1=x(18);

j=19;
k=1;
while j< (18*4)
regulators(j: j+3)=[x(j) x(j+2) x(j+2) x(j+3)];
j=j+5;

end

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

function [DHXT, DMTH1, DMIG1, DSTD1, DMIG2] =system(Glucose, MTH1, MIG1, STD1, MIG2, HXT) 
cytMIG1=1-MIG1;
degMT=maxdegMT*(Glucose)^hilldegMT /(KdegMT^hilldegMT+Glucose^hilldegMT);
degST=maxdegST*(Glucose)^hilldegST /(KdegST^hilldegST+Glucose^hilldegST);
deg=VdegHXT/(1+ (Glucose/threshdegHXT)^hilldegHXT)+basaldeg;
DMTH1= 1  - (degMT+ basaldegMT)*MTH1 ;
DMIG1= Vloc*cytMIG1*(Glucose/KMG)^hillMG/(1 + (Glucose/KMG)^hillMG) -deloc*MIG1 -KinhSMG*STD1*MIG1;
DSTD1= 1/(1+(MTH1/KrepMTST)^hillrepMTST) - (degST+ basaldegST)*STD1 ;
DMIG2= 1/(1+ (MTH1/KrepSTM2)^hillrepSTM2) -basaldegM2*MIG2   -  (VinactGM2 / (1 + (Glucose/KinactGM2)^hillinactGM2))*MIG2; 
DHXT=VHXT/ (1+ (MIG1/KrepHMG)^hillHMG +(MTH1/KrepHMT)^hillHMT +(MIG2/KrepHM2)^hillrepHM2 +(STD1/KrepHST)^hillrepHST) - deg*HXT;

end
 

function [DHXT, DMTH1, DMIG1, DSTD1, DMIG2] =systemSNF3(Glucose, MTH1, MIG1, STD1, MIG2, HXT) 
cytMIG1=1-MIG1;
degMT=maxdegMT*(Glucose)^hilldegMT /(.4^hilldegMT+Glucose^hilldegMT);
degST=maxdegST*(Glucose)^hilldegST /(KdegST^hilldegST+Glucose^hilldegST);
deg=VdegHXT/(1+ (Glucose/threshdegHXT)^hilldegHXT)+basaldeg;
DMTH1= 1  - (degMT+ basaldegMT)*MTH1 ;
DMIG1= Vloc*cytMIG1*(Glucose/KMG)^hillMG/(1 + (Glucose/KMG)^hillMG) -deloc*MIG1 -KinhSMG*STD1*MIG1;
DSTD1= 1/(1+(MTH1/KrepMTST)^hillrepMTST) - (degST+ basaldegST)*STD1 ;
DMIG2= 1/(1+ (MTH1/KrepSTM2)^hillrepSTM2)  -basaldegM2*MIG2   -  (VinactGM2 / (1 + (Glucose/KinactGM2)^hillinactGM2))*MIG2; 
DHXT=VHXT/ (1+ (MIG1/KrepHMG)^hillHMG +(MTH1/KrepHMT)^hillHMT +(MIG2/KrepHM2)^hillrepHM2 +(STD1/KrepHST)^hillrepHST) - deg*HXT;

 
 end

function [DHXT, DMTH1, DMIG1, DSTD1, DMIG2] =systemRGT2(Glucose, MTH1, MIG1, STD1, MIG2, HXT) 
%RGT2 system has a Kdeg for STD1 of 10
    cytMIG1=1-MIG1;
degMT=maxdegMT*(Glucose)^hilldegMT /(KdegMT^hilldegMT+Glucose^hilldegMT);
degST=maxdegST*(Glucose)^hilldegST /(10^hilldegST+Glucose^hilldegST);
deg=VdegHXT/(1+ (Glucose/threshdegHXT)^hilldegHXT)+basaldeg;
DMTH1= 1  - (degMT+ basaldegMT)*MTH1 ;
DMIG1= Vloc*cytMIG1*(Glucose/KMG)^hillMG/(1 + (Glucose/KMG)^hillMG) -deloc*MIG1 -KinhSMG*STD1*MIG1;
DSTD1= 1/(1+(MTH1/KrepMTST)^hillrepMTST) - (degST+ basaldegST)*STD1 ;
DMIG2= 1/(1+ (MTH1/KrepSTM2)^hillrepSTM2) -basaldegM2*MIG2   -  (VinactGM2 / (1 + (Glucose/KinactGM2)^hillinactGM2))*MIG2; 
DHXT=VHXT/ (1+ (MIG1/KrepHMG)^hillHMG +(MTH1/KrepHMT)^hillHMT +(MIG2/KrepHM2)^hillrepHM2 +(STD1/KrepHST)^hillrepHST) - deg*HXT;

 
 end

dout=zeros((18+18*4), 1);
Glucose=Glucose1*.2;

f=19; %variable to store all the regulators
[DHXTMTH1_2,dout(f),dout(f+1),dout(f+2),dout(f+3)] =system(Glucose,0,x(f+1),x(f+2),x(f+3),HXTMTH1_2) ;
dout(f)=0;
f=f+4;
[DHXTMIG1_2,dout(f),dout(f+1),dout(f+2),dout(f+3)] =system(Glucose,x(f),0,x(f+2),x(f+3), HXTMIG1_2) ;
dout(f+1)=0;
f=f+4;
[DHXTSTD1_2,dout(f),dout(f+1),dout(f+2),dout(f+3)] =system(Glucose,x(f),x(f+1),0,x(f+3),HXTSTD1_2); 
dout(f+2)=0;
f=f+4;
[DHXTRGT2_2,dout(f),dout(f+1),dout(f+2),dout(f+3)] =systemRGT2(Glucose,x(f),x(f+1),x(f+2),x(f+3), HXTRGT2_2) ;
f=f+4;
[DHXTSNF3_2,dout(f),dout(f+1),dout(f+2),dout(f+3)] =systemSNF3(Glucose,x(f),x(f+1),x(f+2),x(f+3), HXTSNF3_2) ;
f=f+4;
[DHXT_2,dout(f),dout(f+1),dout(f+2),dout(f+3)] =system(Glucose,x(f),x(f+1),x(f+2),x(f+3), HXT_2) ;


Glucose=Glucose1*.4;

f=f+4;
[DHXTMTH1_4,dout(f),dout(f+1),dout(f+2),dout(f+3)] =system(Glucose,0,x(f+1),x(f+2),x(f+3), HXTMTH1_4) ;
dout(f)=0;
f=f+4;
[DHXTMIG1_4,dout(f),dout(f+1),dout(f+2),dout(f+3)] =system(Glucose,x(f),0,x(f+2),x(f+3), HXTMIG1_4) ;
dout(f+1)=0;
f=f+4;
[DHXTSTD1_4,dout(f),dout(f+1),dout(f+2),dout(f+3)] =system(Glucose,x(f),x(f+1),0,x(f+3), HXTSTD1_4); 
dout(f+2)=0;
f=f+4;
[DHXTRGT2_4,dout(f),dout(f+1),dout(f+2),dout(f+3)] =systemRGT2(Glucose,x(f),x(f+1),x(f+2),x(f+3),HXTRGT2_4) ;
f=f+4;
[DHXTSNF3_4,dout(f),dout(f+1),dout(f+2),dout(f+3)] =systemSNF3(Glucose,x(f),x(f+1),x(f+2),x(f+3),HXTSNF3_4) ;
f=f+4;
[DHXT_4, dout(f),dout(f+1),dout(f+2),dout(f+3)] =system(Glucose,x(f),x(f+1),x(f+2),x(f+3), HXT_4) ;


Glucose=Glucose1;
f=f+4;
[DHXTMTH1_1,dout(f),dout(f+1),dout(f+2),dout(f+3)] =system(Glucose, 0,x(f+1),x(f+2),x(f+3), HXTMTH1_1) ;
dout(f)=0;
f=f+4;
[DHXTMIG1_1,dout(f),dout(f+1),dout(f+2),dout(f+3)] =system(Glucose,x(f),0,x(f+2),x(f+3),HXTMIG1_1) ;
dout(f+1)=0;
f=f+4;
[DHXTSTD1_1,dout(f),dout(f+1),dout(f+2),dout(f+3)] =system(Glucose,x(f),x(f+1),0,x(f+3),HXTSTD1_1); 
dout(f+2)=0;
f=f+4;
[DHXTRGT2_1,dout(f),dout(f+1),dout(f+2),dout(f+3)] =systemRGT2(Glucose,x(f),x(f+1),x(f+2),x(f+3), HXTRGT2_1) ;
f=f+4;
[DHXTSNF3_1,dout(f),dout(f+1),dout(f+2),dout(f+3)] =systemSNF3(Glucose,x(f),x(f+1),x(f+2),x(f+3),HXTSNF3_1) ;
f=f+4;
[DHXT_1,dout(f),dout(f+1),dout(f+2),dout(f+3)] =system(Glucose,x(f),x(f+1),x(f+2),x(f+3), HXT_1) ;


dout(1)=DHXT_2;
dout(2)=DHXT_4;
dout(3)=DHXT_1;
dout(4)=DHXTMIG1_2;
dout(5)=DHXTMIG1_4;
dout(6)=DHXTMIG1_1;
dout(7)=DHXTMTH1_2;
dout(8)=DHXTMTH1_4;
dout(9)=DHXTMTH1_1;
dout(10)=DHXTSTD1_2;
dout(11)=DHXTSTD1_4;
dout(12)=DHXTSTD1_1;
dout(13)=DHXTRGT2_2;
dout(14)=DHXTRGT2_4;
dout(15)=DHXTRGT2_1;
dout(16)=DHXTSNF3_2;
dout(17)=DHXTSNF3_4;
dout(18)=DHXTSNF3_1;




end 
end