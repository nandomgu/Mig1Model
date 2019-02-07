function modelh= Std1Test(params, mediaInput)
%%this model is of the negative self regulation of protein Std1 
%%%%% Std1 protein self repressed, glucose-induced degradation and basal degradation too.


modelh= @transcriptModelh;

function dout=transcriptModelh(t, x)

mediaInput=smooth(mediaInput);
Glucose=interp1(0:(numel(mediaInput)-1), mediaInput,t);


Std1=x(1);



VStd1=params(1);
KSelfStd1=params(2);
KDegStd1=params(3);
threshStd1=params(4);

DStd1= VStd1- basalDegStd1*Std1 -(KDegStd1*Glucose/ (threshStd1+Glucose)*Std1); 



dout= [];

dout(1)=DStd1;


end
end