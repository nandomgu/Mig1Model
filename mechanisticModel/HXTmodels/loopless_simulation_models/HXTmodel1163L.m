
function modelh= HXTmodel1163L(params, mediaInput)
%model 114. third model (after model 112 and 113) that includes 2 sensors processing the glucose plus
%most of the items from model 111, minus threshold induced degradation of
% HXT
%input. In particular, hthis model also loses threshold/induced activation
%of mRNA.
%Model 115. reduces the above model by temoving the Hill functions and adds
%degradation of mth1 and std1. in this case, std1 and mth1 are produced
%constitutively and degraded by glucose to different extents through snf3
%and rgt2. mth1 and std1 actively repress the mRNA, each of them to
%different extents. 
%Model 117. model 116 gets simplified parameter-wise, and as repression can
%only by maximal (vmax=0) then the relative contributions to repression by
%Mth1 and std1 can be defined as relative to each other. 

%Model1162 has mth1 and std1 , and their degradation dependent on hill
%functions. 
%model 1163 includes snf3 and rgt2 as gates

modelh= @transcriptModelh;

function dout=transcriptModelh(t, x)

%mediaInput=[0;0.0223032063571172;0.0405047314534221;0.0577574155317936;0.0708017171371265;0.074879900620598;0.0758356059661576;0.0771340928584555;0.0770603734213315;0.0766017074310701;0.0825233317190962;0.0909744554418606;0.093993550248502;0.102173022231747;0.116563929130705;0.143203882900438;0.176368881684685;0.223713693914748;0.281052257734873;0.346466307100618;0.407476155958698;0.471069387757387;0.533242734001389;0.593344344903071;0.646896154803748;0.701085931884189;0.738849711878534;0.778092173550209;0.811575003034639;0.840202031942638;0.856290962143873;0.882968428322724;0.897828449866704;0.904189374798246;0.904315649063733;0.906949983540048;0.904302812442446;0.89571826831523;0.884935318073253;0.875509014907067;0.86321464541776;0.846263132839432;0.83026705376725;0.813688756593054;0.79768698434489;0.778942343617509;0.756122945858676;0.728587294009643;0.695907127396961;0.653422645982901;0.598590983814381;0.539549562479458;0.479394232916031;0.414454994302869;0.348601794470689;0.285308136848455;0.224213405996943;0.167090244832356;0.122325009547493;0.0901200542782974;0.0726337968136678;0.0647512577549924;0.0600383055243802;0.0550657572829314;0.0492174819494744;0.0441327496532019;0.0404784360648094;0.0386288373598618;0.0392827393567328;0.0416060245174376;0.04385917090531;0.0492989753395697;0.0573014894896179;0.0739614838095577;0.0994063835142636;0.131546699085503;0.168006119875606;0.20781938996431;0.248857320083113;0.290178889871922;0.335607409038797;0.388906458123698;0.449891255385504;0.511048599375835;0.575290050983856;0.639074832665356;0.701137058482281;0.764278106678239;0.82464351652881;0.87528951712772;0.91648686201527;0.945361227023334;0.961844754355451;0.972233559512459;0.981067057337478;0.98776798745466;0.989072093065479;0.983223292248401;0.975825441497439;0.959369224507547;0.941174416237702;0.927473435570653;0.914943640537118;0.902157327430209;0.893530447953199;0.882762771839404;0.869213756591971;0.858531802932826;0.844644637638763;0.820074164859752;0.790494804327956;0.752437898339734;0.701541856613617;0.635115335626438;0.57337245957311;0.506332572458254;0.437997050295204;0.374797729083169;0.318603618410695;0.263920417707352;0.219355412800032;0.181480990485946;0.147997849252685;0.123091198645273;0.104051523676733;0.0870487941621835;0.0738718773828963;0.0639971956947056;0.0562793281769372;0.0495754666431618;0.0448937810251198;0.0412568388348441;0.0384279360909205;0.0358014676960468;0.0335895169856507;0.0316324820712432;0.0297321354517882;0.0275405856334663;0.0256602563979294;0.0234919909693372;0.0217603007531544;0.0201354702808326;0.0186647892246733;0.017383650443999;0.016511174268384;0.0152562395520462;0.0142146133505696;0.013299036415849;0.0122787686065637;0.0113007857746887;0.0104571487391209;0.00966251216427188;0.00902375582022215;0.00839836325494823;0.00778158911118866;0.00718619652525696];
mediaInput=smooth(mediaInput);
Glucose=interp1(0:5:(numel(mediaInput)-1)*5, mediaInput,t);

mRNA= x(1);
HXT= x(2);
MatureHXT= x(3);
Mth1=x(4);
Std1=x(5);

maxVmRNA=params(1);
KdegmRNA=params(2);
alphaHXT=params(3);
maxdegHXT= params(4);
KmatHXT=params(5);
KBleach=params(6); 
degThreshHXTbelow= params(7);
  KdegMth1=params(8);
  KdegStd1=params(9); 
VMth1= params(10);
maxdegMth1= params(11);
VStd1= params(12);
maxdegStd1= params(13);
RepThreshMth1=params(14);
RepThreshStd1=params(15);
RepThreshMig1=params(16);
hillDegMth1=params(17);
hillDegStd1=params(18);
hillRepMig1=params(19);
hillSnf3=params(20);
hillRgt2=params(21);
KSnf3=params(22);
KRgt2=params(23);
         

 
snf3Gate= 1/(1+(Glucose/KSnf3)^hillSnf3);
rgt2Gate= (Glucose)/(KRgt2+ Glucose);
DMth1= VMth1-  snf3Gate*(maxdegMth1/(1+(Glucose^hillDegMth1/KdegMth1^hillDegMth1)))*Mth1;
DStd1= VStd1-  rgt2Gate*(maxdegStd1*(Glucose^hillDegStd1/KdegStd1^hillDegStd1+Glucose^hillDegStd1))*Std1;


VmRNASum= maxVmRNA/(1+(Mth1/RepThreshMth1))+ maxVmRNA/(1+(Std1/RepThreshStd1))+maxVmRNA/(1+(Glucose^hillRepMig1/RepThreshMig1^hillRepMig1));

DmRNA= VmRNASum-KdegmRNA*mRNA;
DHXT= alphaHXT*mRNA-(maxdegHXT/1+(degThreshHXTbelow/(Glucose)))*HXT;

  DMatureHXT= KmatHXT;
  
  
 dout= zeros(5,1);

 dout(1)=DmRNA;
 dout(2)=DHXT;
 dout(3)=DMatureHXT*HXT-(KBleach*MatureHXT);
 dout(4)=DMth1;
 dout(5)= DStd1;
 
 
end 
end