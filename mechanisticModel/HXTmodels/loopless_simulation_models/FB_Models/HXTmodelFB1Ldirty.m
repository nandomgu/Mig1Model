function modelh= HXTmodelFB1Ldirty(params, mediaInput)
%%This new family of models inserts a feedback into the network


modelh= @transcriptModelh;

function dout=transcriptModelh(t, x)

mediaInput=smooth(mediaInput);
Glucose=interp1(0:(numel(mediaInput)-1), mediaInput,t);

A= x(1);
B= x(2);



% VmRNA=params(1);
% KmRNA=params(2);
% hillmRNA=params(3);
% KdegmRNA=params(4);
% alphaHXT=params(5);
% KdegHXT= params(6);
% KactiveMig1=params(7);
% RepThreshMig1=params(8);
% hillactiveMig1=params(9);

KGA=params(1);
KBA=params(2);
KAB=params(3);
KdegB=params(4);


%%% traditional equations
%A receives the glucose input  and the action of B.
%DA= A^hillAB/(KAB +A^hillAB)+ Glucose^hillGA/(KGA + Glucose^hillGA)-KdegA*A;
%B receives A, either in an activation or repression fashion
%DB= A^hillAB/(KAB +A^hillAB)-KdegB*B;
%repression     DB= 1/ (1 + A^hillAB/KAB^hillAB)-KdegB*B;
%C is ideally the mRNA
%DC= A^hillAC/(KAC + A^hillAC)-KdegC*C;

%equations from ma et al



DA= KGA*Glucose*(1-A)-KBA*A*B;

DB= KAB*A* ((1-B)/(KAB+1-B))- KdegB* (B/(KdegB+B));





%  activeMig1= (Glucose^hillactiveMig1)/(KactiveMig1^hillactiveMig1+ Glucose^hillactiveMig1);
%  DmRNA= (VmRNA*(Glucose)^hillmRNA /(KmRNA^hillmRNA+Glucose^hillmRNA))/(1+(activeMig1/RepThreshMig1))+ -KdegmRNA*mRNA;
%  DHXT= alphaHXT*mRNA- KdegHXT*HXT;
 
 
 dout= zeros(2,1);
 dout(1)=DB;
 dout(2)=DA;

 
end 
end