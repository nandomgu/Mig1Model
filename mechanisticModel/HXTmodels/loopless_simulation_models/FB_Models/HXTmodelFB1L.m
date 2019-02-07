function modelh= HXTmodelFB1L(params, mediaInput)
%%This new family of models inserts a feedback into the network


modelh= @transcriptModelh;

function dout=transcriptModelh(t, x)

mediaInput=smooth(mediaInput);
Glucose=interp1(0:(numel(mediaInput)-1), mediaInput,t);

A= x(1);
B= x(2);


KactivGA=params(1);
KdeactivBA=params(2);
KactivAB=params(3);
KdeactivBA=params(4);
KThresh3=params(5)
KThresh4=params(6);


DA= KactivGA*Glucose*(1-A)-KdeactivBA*A*B;
DB= KactivAB*A* ((1-B)/(KThresh3+1-B))- KdegB* (B/(KThresh4+B));

 
 dout= zeros(2,1);
 dout(1)=DB;
 dout(2)=DA;

 
end 
end