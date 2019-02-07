
function modelh= simpleModel1(params, args)
%%Model 1 is simple threshold-activated production of Hxt, including the
%%hill hypothesis.
%preferedOutVar=2
mediaInput=args.input;
times=args.times;
modelh= @transcriptModelh;

function dout=transcriptModelh(t, x)

%mediaInput=smooth(mediaInput);
Glucose=interp1(times, mediaInput,t);

mRNA= x(1);

VmRNA=params(1);
KmRNA=params(2);
hillmRNA=params(3);
KdegmRNA=params(4);


DmRNA= (VmRNA*(Glucose)^hillmRNA /(KmRNA^hillmRNA+Glucose^hillmRNA)) -KdegmRNA*mRNA;
 %(1+ (Glucose/KrepmRNA)^hillrepmRNA))
dout= zeros(1,1);
dout(1)=DmRNA;

 
 
end 
end