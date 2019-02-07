


% mediaInput=zeros(156,1);
%mediaInput(100:end)=1;
mediaInput(1:50)=linspace(0,1,50);
mediaInput(50:100)=1;
mediaInput(100:156)=linspace(1,0,57);
KI= 0.7;
	
 KA = 3;

 N_TAR =7;
	
 ALFA_TAR_MEASP = 2;
	
 M0_TAR_MEASP = 1;
	% parameter in equation (3)
KR = 0.005;
	% 1/s;  linear rate for methylation process
 KB = KR;
	% 1/s;  linear rate for demethylation process
L=mediaInput;
m=1;
M=length(mediaInput);
for i=1:length(mediaInput)
    
 a(i) = 1 ./ (1 + exp(N_TAR .* ((ALFA_TAR_MEASP .* (M0_TAR_MEASP - m))...
        - (log((1 + L(i)./KA)./(1 + L(i)./KI))))));



    m = m + (KR .* (1 - a(i)) - KB .* a(i)) .* 1;
end
%for i=1:length(mediaInput)
 %   if a(i)<0
  %      a(i)=0;
   % end
%end
figure(10)
plot(mediaInput)
figure(11)
plot(a)


rsim603L=rampSimHandle3(@HXTmodel603L, a, [0 0 0 0 ], options, 3, gfpOutput')
%%now this is a function that you can use.


[x,y]=rsim603L(NiceParamsHxt4.model603L(1,:));  %x is the least square difference, y is the actual simulated time series.
figure(12); plot(y)

