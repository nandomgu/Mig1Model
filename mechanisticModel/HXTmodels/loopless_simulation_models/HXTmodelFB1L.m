function modelh= HXTmodelFB1L(params, mediaInput)
%%This new family of models inserts a feedback into the network


modelh= @transcriptModelh;

function dout=transcriptModelh(t, x)

mediaInput=smooth(mediaInput);
Glucose=interp1(0:(numel(mediaInput)-1), mediaInput,t);

A= x(1);
B= x(2);


KGA=params(1);
KBA=params(2);
KAB=params(3);
KdegB=params(4);


DA= KGA*Glucose*(1-A)-KBA*A*B;
DB= KAB*A* ((1-B)/(KAB+1-B))- KdegB* (B/(KdegB+B));

 
 dout= zeros(2,1);
 dout(1)=DB;
 dout(2)=DA;

 
end 
end