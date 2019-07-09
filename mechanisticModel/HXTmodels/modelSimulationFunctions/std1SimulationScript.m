

% VStd1=params(1);
% KSelfStd1=params(2);
% KDegStd1=params(3);
% threshStd1=params(4);
% basalDegStd1=params(5);
% VMth1= (6);
% basalDegMth1=params(7);
% KDegMth1=params(8);
% threshMth1=params(9);
% KRepMth1= params(10);
mf=extractModelFeatures('Std1Test')
options=mf.options
% one way of storing the parameters could be a structure like below


% params.VStd1=1;
% params.KSelfStd1=.01;
% params.KDegStd1=8;
% params.threshStd1=.01;
% params.basalDegStd1=.05;
% params.VMth1=1;
% params.basalDegMth1=.05;
% params.KDegMth1=.3;
% params.threshMth1=.01;
% params.KRepMth1=.1;
%then we retrieve the numers as follows.
%a=struct2cell(params); ps=log([a{:}]); 
%%the disadvantage of doing this is that you are stuck
%%with using all the parameters. therefore a second way to store the params
%would be in a cell array. 

params={};
params(1,:)= {'VStd1', 1, [-4,4]};%maximum production rate of Std1
params(2,:)= {'KSelfStd1', .4 , [-4, 4]};% threshold at which Std1 represses itself.
params(3,:)= {'KDegStd1', 1, [-4,4]}; % maximum glucose induced degradation rate of Std1
params(4,:)= {'threshStd1', .01, [-4, 4]}; % threshold at which Std1 is degraded by glucose
params(5,:)= {'basalDegStd1', .05, [-4, 4]}; %basal degradation rate for Std1
params(6,:)= {'VMth1', 1, [0, 1]}; %Maximum production rate for Mth1
params(7,:)= {'basalDegMth1', .05, [-4, 4]}; %basal degradation rate for Mth1
params(8,:)= {'KDegMth1', 1, [-4 4]};% maximum glucose induced deg. rate for Mth1
params(9,:)= {'threshMth1', .01, [-4, 4]}; %% threshold for Mth1 to be degraded by glucose
params(10,:)= {'KRepMth1', .01, [-4, 4]}; %% Threshold for Mth1 to repress Std1
params(11,:)= {'thresholdStd1Mig1', .1, [-4, 4]}; %%threshold for std to inhibit mig1/ raise the activation threshold of mig1
params(12,:)={'hillactiveMig1', 4, [ 1, 4]}; %hill coefficient for mig1 to become active upon glucose
params(13,:)={'KSnf1', .01, [-6 0]}; %%threshold at which glucose inhibits snf1
params(14,:)={'KRepSnf1Mig1', .9, [-6, 6]}; %%threshold of snf1 needed to inhibit Mig1
params(15,:)={'KEnhance', .1, [-6,6]}; %threshold at which Std1 enhances Snf1
params(16,:)={'KMig1HXT4', .005,[-6,6]};
params(17,:)={'KStd1HXT4', .2,[-6,6]};
params(18,:)={'KMth1HXT4', .2,[-6,6]};
params(19,:)={'KMig1HXT2', 5,[-6,6]};
params(20,:)={'KMth1HXT2', 15,[-6,6]};
params(21,:)={'KDegSnf1', .1, [-6, 0]};
params(22,:)={'basalDegSnf1', 0.1, [-6, 0]};
params(23,:)={'basalDegHXT2', 0.05, [-6 0]} 

%for j=1:size(params, 1)
%params{j,2}= randRange(params{j,3})
%end

[r, inds]= ismember(mf.paramNames, params(:,1)')

ps= log([params{inds, 2}]);
incons=zeros(1, mf.numInit); %%get all possible initial conditions
incons(1:4)=[0,20, 0, 0];
outvarnames={'HXT2'};
[~, outvar]= ismember( outvarnames, mf.varnames) %%insinde the first brackets add the variable names you are interested in. this looks for std1in variable names, and the second output is the index which is.
%[~, whichMig1]= ismember( 'activeMig1', outvarnames)



figure; 
expand=50;
totrow=5;
rampcolors= cmapFromTo([0 0 0],[0, .8, .6],totrow);
%responsecolors=cmapFromTo([],[],10);

for j=1:totrow
%for k=1:10
%input=fullTransition(300,0,50+expand*k, expand*j)
input=fullTransition(300,expand*j,200, expand*j)
rsimStd1= rampSimHandle2(@Std1Test, input, incons,options, outvar)
[~,y,t]=rsimStd1(ps)
%if whichMig1~=0
%df=diff(y(:, whichMig1));
%y(:, whichMig1)= [df;df(end)];
%end
%subplot(,10, (j-1)*10+k);
axes();
 plot(t, y, 'Color', rampcolors(j,:));addHLine(0); legend([outvarnames]); legend('hide')
%addHLine(12.53)
excedent{j}=trapz(y(y>12.53))
yyaxis right; plot(input, 'Color', rampcolors(j,:));addVLine(300+expand*j+expand+expand*j);
yyaxis left;
pause(.5)
%end
end

stackPlots(gcf)

%.2*Snf1* Glucose- 0.05*Snf1

%%%%
%%%%
%%%% changing 2 parameters in paralel

incons=zeros(1, mf.numInit); %%get all possible initial conditions
incons(1:4)=[0,20, 0, 1];
outvarnames={'Snf1'};
[~, outvar]= ismember( outvarnames, mf.varnames) %%insinde the first brackets add the variable names you are interested in. this looks for std1in variable names, and the second output is the index which is.
[~, whichMig1]= ismember( 'activeMig1', outvarnames)


figure; 
axes(); %generate only one axes
expand=50;
totrow=50;
rampcolors= cmapFromTo([0 0 0],[0, .8, .6],totrow);
%responsecolors=cmapFromTo([],[],10);
p1Name='KDegSnf1'; 
p2Name='KRepMth1'; 
p3Name='KMig1HXT4';
p4Name='KDegMth1';
changeParamNames={p1Name, p2Name, p3Name, p4Name};

[is, p1]=ismember(p1Name, params(:,1));
[is, p2]=ismember(p2Name, params(:,1)); %%fetching the number of the parameters
[is, p3]=ismember(p3Name, params(:,1)); %%fetching the number of the parameters
[is, p4]=ismember(p4Name, params(:,1)); %%fetching the number of the parameters


changeParams= [p1];


vals=linspace(-10,0, totrow); %values to try for these parameters

for j=1:totrow
%for k=1:10
%input=fullTransition(300,0,50+expand*k, expand*j)
input=fullTransition(300,150,150, 150)
rsimStd1= rampSimHandle2(@Std1Test, input, incons,options, outvar)

[r, inds]= ismember(mf.paramNames, params(:,1)');
ps= log([params{inds, 2}]);
ps(changeParams)= vals(j);
[~,y,t]=rsimStd1(ps)
if whichMig1~=0
df=diff(y(:, whichMig1));
y(:, whichMig1)= [df;df(end)];
end
%subplot(,10, (j-1)*10+k);
plot(t, y, 'Color', rampcolors(j,:));
hold on
addHLine(0); legend(mf.varnames{outvar}); legend('hide')
%addHLine(12.53)
excedent{j}=trapz(y(y>12.53))
%pause(.5)
end
yyaxis right; plot(input, 'Color', rampcolors(1,:));
yyaxis left;
%textbp( strjoin(changeParamNames{changeParams}, 'Delimiter', ' '))
%for which only parameter 1 is in a range
textbp([ p1Name '[' num2str(exp(p1Range(1))) ',' num2str(exp(p1Range(end))) ']'  ]) ;

%%%
%%%
%%%  modeling several concentration steps

incons=zeros(1, mf.numInit); %%get all possible initial conditions
incons(1:4)=[0,20, 0, 1];
outvarnames={'HXT2'};
[~, outvar]= ismember( outvarnames, mf.varnames) %%insinde the first brackets add the variable names you are interested in. this looks for std1in variable names, and the second output is the index which is.
[~, whichMig1]= ismember('activeMig1', outvarnames)

%stepsFigure=figure;
%figure; 
%axes(); %generate only one axes
expand=50;
totrow=10; %%total concentrations
totvals=1;
vals=linspace(-6,3, totrow); %values to try for these parameters and concentrations
pvals=linspace(0,3, totvals);

rampcolors= cmapFromTo([0 0 1],[1, 0, 1],totrow); %blue to magenta for concentration
%responsecolors=cmapFromTo([],[],10);
p1Name='basalDegSnf1'; 
p2Name='threshStd1'; 
p3Name='KMth1HXT2';
p4Name='KMig1HXT2';
changeParamNames={p4Name};

[is, p1]=ismember(p1Name, params(:,1));
[is, p2]=ismember(p2Name, params(:,1)); %%fetching the number of the parameters
[is, p3]=ismember(p3Name, params(:,1)); %%fetching the number of the parameters
[is, p4]=ismember(p4Name, params(:,1)); %%fetching the number of the parameters

changeParamNumbers=[p4];
paramsOfInterest=[1];
changeParams= [changeParamNumbers(paramsOfInterest)];

figure;
for j=1:totvals %looping over param values
axes();

for k=1:totrow %%looping over concentrations

%input=fullTransition(300,0,50+expand*k, expand*j)
input=fullTransition(600,0,600,0)* k/totrow %%steps of 300 minutes. maximum concentration for the sugar is one
rsimStd1= rampSimHandle2(@Std1Test, input, incons,options, outvar)
[r, inds]= ismember(mf.paramNames, params(:,1)');
ps= log([params{inds, 2}]);
ps(changeParams)= pvals(j);
[~,y,t]=rsimStd1(ps)
%if whichMig1~=0
%df=diff(y(:, whichMig1));
%y(:, whichMig1)= [df;df(end)];
%end
%subplot(,10, (j-1)*10+k);
plot(t, y, 'Color', rampcolors(k,:));
hold on
addHLine(0); %legend(mf.varnames{outvar}); legend('hide')
%addHLine(12.53)
%excedent{j}=trapz(y(y>12.53))
%pause(.5)


end
textbp( [[changeParamNames{paramsOfInterest}], ' =' num2str(exp(pvals(j)))]);
end
stackPlots(gcf)



yyaxis right; plot( input, 'Color', rampcolors(1,:));
yyaxis left;
textbp( [outvarnames{:}])




figure; plot(diff(y)); yyaxis right; plot(input)

%%%
%%% parameter variation subscript
%%%
%%looping through all params
allResponses=cell(numel(params(:,1)),numel(params(:,1)))
allRatios=cell(numel(params(:,1)),numel(params(:,1)))
%for b=1:numel(params(:,1))
%for d=2:numel(params(:,1))


%%assigning parameters to work with
p1Name='KMig1HXT2'; 
p2Name='KMth1HXT2'; 


%% making sure you only do stuff once
%if b==d || ~isempty(allResponses{b,d}) || ~isempty(allResponses{d,b})
%continue;
%end

%%assigning the parameters from looping
%p1Name=params(b,1);
%p2Name=params(d, 1); 

%%finding what numbers these parameters are. redundant for looping, but anyway.
[is, p1]=ismember(p1Name, params(:,1))

[is, p2]=ismember(p2Name, params(:,1))


%%assigning num of columns and rows
totrow=10;
totcol=10;

%%assigning ranges for those parmeters

p1Range=linspace(params{p1, 3}(1), params{p1, 3}(2), totrow)
p2Range=linspace(params{p2, 3}(1), params{p2, 3}(2), totcol)


%assinging starting values to parameters
[r, inds]= ismember(mf.paramNames, params(:,1)')
ps= [params{inds, 2}];

%% initial conditions
incons=zeros(1, mf.numInit); %%get all possible initial conditions
incons(1:4)=[0,20, 0, 1];
%%we look for the variable we want to show

outvarnames={'HXT2'};
[~, outvar]= ismember( outvarnames, mf.varnames) %%insinde the first brackets add the variable names you are interested in. this looks for std1in variable names, and the second output is the index which is.
[~, whichMig1]= ismember('activeMig1', outvarnames)
[~, whichSnf1]= ismember('Snf1', outvarnames)
h=figure; 
responses={}
upDownRatios=[];
c=.95;
for j=1:totrow
for k=1:totcol
%input=fullTransition(300,0,50+expand*k, expand*j) %in case we want to change ramp high glucose time
%input=fullTransition(300,0,50, expand*j) % in case we want to change the downshift duration
input=fullTransition(300,150,150, 150) %%fixed ramping down time
 rsimStd1= rampSimHandle2(@Std1Test, input, incons,options, outvar)
 [r, inds]= ismember(mf.paramNames, params(:,1)')
 ps= log([params{inds, 2}]); %temporary line to get the log of the numbers
 ps(p1)= p1Range(j); %substitution params are already in log scale
 ps(p2)= p2Range(k);
[~,y,t]=rsimStd1(ps)

%try
%df=diff(y(:, whichMig1));
%y(:, whichMig1)= [df;df(end)];
%df=diff(y(:, whichSnf1));
%y(:, whichSnf1)= [df;df(end)];
%end

responses{j,k}=y;
[~, ~, ~, s]=getAllTransitionIndices(input, .001, .999); %getting the ramping sections

upSection=y((s==2));
downSection=y((s==4));
upTime=numel(t(s==2));
downTime=numel(t(s==4));

upDownRatios(j,k)= log((trapz(downSection)/upTime) / (trapz(upSection)/downTime));
%subplot(totrow,totcol, (j-1)*totcol+k);
axes(h);
 plot(t, y);addHLine(0); legend(mf.varnames{outvar}); legend('hide')
 set(gca,'Position', [(1/totcol)*(k-1), (j/totrow)-(1/totrow), (1/totcol)*c, (1/totrow)*c]) 
 xlabel([p1Name,'= ' num2str(exp(ps(p1)))])
 ylabel([p2Name, '= ' num2str(exp(ps(p2)))])
 %text([num2str(p2Range(k)) ',' num2str(p1Range(j))])
yyaxis right; plot(input);
yyaxis left;
end
end

allResponses{p1, p2}=responses;
allRatios{p1,p2}= upDownRatios;
allp1Ranges{p1,p2}= p1Range;
allp2Ranges{p1,p2}=p2Range;

end
end


%axes(l); 
%set(gca,'Position', [(1/numel(ps))*(p2-1), (p1/numel(ps))-(1/numel(ps)), (1/numel(ps))*c, (1/numel(ps))*c]) 
%imshow(imresize(upDownRatios,200,'nearest'), [],'colormap', colorPalette([.5 .3 0; 0 0 0; 0 .8 .5],20))
%xlabel(p1Name)
%ylabel(p2Name)
%zlabel('log(Downshift/Upshift)')
%hideAllTicks(gcf)


%%plotting all lines of one parameter for one value of the other parameter

figure; 
lgnds={};
m= allResponses{p1, p2};
totrow=size(m, 1);
totcol=size(m,2)

 cm=cmapFromTo([0 .1 .6],[1 0 1],  totcol);%% colormap used for std1
for j=1:totrow
axes();
for k=1:totcol
	hold on; 
	
 	plot(t, m{j,k}, 'Color', cm(k,:));addHLine(0);lgnds{k}=[p2Name '=' num2str(exp(p2Range(k)))];
 	%legend('show'); lg=legend(gca); lg.String{k}=[p2Name '=' num2str(exp(p2Range(j)))]; legend('hide') %pain in the ass line to manipulate legends
	set(gca,'Position', [0, (j/totrow)-(1/totrow), 1, (1/totrow)*c]) %%plotting rows from bottom (j=1) to top
	
end
yyaxis right; plot(input);
%legend(gca, lgnds{:}, cm)
textbp([p1Name '=' num2str(exp(p1Range(j)))]);
textbp([ p2Name '[' num2str(exp(p2Range(1))) ',' num2str(exp(p2Range(end))) ']'  ]) ;
end
hideAllTicks()



%%only plotting the grid a particular combination of parameters

figure; 

m= allResponses{p1, p2};
totrow=size(m, 1);
totcol=size(m,2)
for j=1:totrow
for k=1:totcol
	axes();
 	plot(t, m{j,k});addHLine(0); legend(mf.varnames{outvar}); legend('hide')
	set(gca,'Position', [(1/totcol)*(k-1), (j/totrow)-(1/totrow), (1/totcol)*c, (1/totrow)*c]) %%plotting rows from bottom (j=1) to top
	
end
end
hideAllTicks()






%%testing mig activation functions 
%%we could alternatively 
Snf1=linspace(0,1, 1000);
activeMig1= @(Snf1, Std1, KRepSnf1Mig1, KEnhance) 1/ (1+ (Snf1/KRepSnf1Mig1)^4+(Std1/KEnhance));
Std1=linspace(0,1,1000);
sn= repmat(Snf1, 1000,1); %Snf1 goes 01 from left to right
st= flipud(repmat(Std1, 1000,1)'); %makes std go 0-1 from the bottom up
figure;
totrow=4;
totcol=10;
c=.9;
colRange=linspace(-4,4,totcol);
rowRange=linspace(-4,0,totrow);
for j=1:totrow
for k=1:totcol
krsm= repmat(exp(rowRange(j)), 1000, 1000); %KRepSnf1Mig1
ke= repmat(exp(colRange(k)), 1000, 1000); %KEnhance
disp(j)
axes()
set(gca,'Position', [(1/totcol)*(k-1), (j/totrow)-(1/totrow), (1/totcol)*c, (1/totrow)*c]) 
imshow(arrayfun(activeMig1, sn, st, krsm, ke), [0,1], 'colormap', parula);
xlabel(exp(colRange(k)))
ylabel(exp(rowRange(j)))
end 

end

%%%checking hxt activation functions



activeMig1= @(Glucose) Glucose/(.1 + Glucose);

Mth1=linspace(0,1, 1000);
Glucose=linspace(0,1,1000);
%hill factor for the mth1 term
%hxt2=@(Mth1, KMth1HXT2, activeMig1, KMig1HXT2) 1/(1+((Mth1/KMth1HXT2)^4)+(activeMig1/KMig1HXT2))
%%hill factor for both mig and mth1
%hxt2=@(Mth1, KMth1HXT2, activeMig1, KMig1HXT2) 1/(1+((Mth1/KMth1HXT2)^4)+(activeMig1/KMig1HXT2)^4)
%% no hill factor at all
hxt2=@(Mth1, KMth1HXT2, activeMig1, KMig1HXT2) 1/(1+((Mth1/KMth1HXT2))+(activeMig1/KMig1HXT2))


%%%it was not obvious to me that if you are below the KMth1HXT2 then the expression is proportional to
%%% glucose but inaccessible to mig1 regulation. if you are above KMth1HXT2  Mth1 may be completely gone
%%% but Mig1 may severely limit the expression. that is why the curve looks nice when you need a lot of mig1
%%%to repress Hxt2 and KMig1 be higher than Mth1. 





gl= repmat(Mth1, 1000,1); %Glucose goes 01 from left to right
mt= flipud(repmat(Glucose, 1000,1)'); %makes std go 0-1 from the bottom up

mig1= arrayfun(activeMig1, gl);

black2green=cmapFromTo([0, 0, 0],[0,1,.6], 100);
black2magenta=cmapFromTo([0, 0, 0],[1,0,1], 100);
figure; 
subplot(1,3,1)
imshow(gl, [0,1])
title('glucose')
subplot(1,3,2)
imshow(mt, [0,1])
title('Mth1')
subplot(1,3,3)

imshow(mig1, [0,1], 'colormap', black2green)
title('active Mig1=F(glucose)')

figure;
totrow=4;
totcol=10;
c=.9;
colRange=linspace(-4,4,totcol);
rowRange=linspace(-4,0,totrow);
for j=1:totrow
for k=1:totcol
kmth1= repmat(exp(rowRange(j)), 1000, 1000); %KRepSnf1Mig1
kmig1= repmat(exp(colRange(k)), 1000, 1000); %KEnhance
disp(j)
axes()
set(gca,'Position', [(1/totcol)*(k-1), (j/totrow)-(1/totrow), (1/totcol)*c, (1/totrow)*c]) 
imshow(arrayfun(hxt2, mt, kmth1, mig1, kmig1), [0,1], 'colormap', black2magenta);
xlabel(exp(colRange(k)))
ylabel(exp(rowRange(j)))
end 

end



%%fitting one particular gene's parameters
[r, inds]= ismember(mf.paramNames, params(:,1)') %%getting parameters present in the model
parsPresent= params(inds, :);
outvarnames={'HXT2'};
dataName={'hxt2'};
[means, nms]=processMean(multichamber20170221) %%obtaining the means for this experiment
[~, datainds]= ismember( dataName, nms); %%which is the index that corresponds to the variables of interest?

[~, outvar]= ismember( outvarnames, mf.varnames); %%insinde the first brackets add the variable names you are interested in. this looks for std1in variable names, and the second output is the index which is.
params2fit= {'KMth1HXT2', 'KMig1HXT2','basalDegHXT2'};
[is, parnumbers]=ismember(params2fit, parsPresent(:,1)); %%which are the numbers of the parameters in the model
fitStartPoint= [parsPresent{parnumbers,2}]

input=fullTransition(180,0,480,0) %%steps of 480 minutes. maximum concentration for the sugar is one
rsimStd1= rampSimHandle2(@Std1Test, input, incons,options, outvar, means(datainds,:)')

%%getting means to fit and names of the variables in experiment

ps= log([parsPresent{:, 2}]);

%%minimisation function
simParamSubset= @(newvalues)  rsimStd1(paramReplace(ps, parnumbers, newvalues));
simParamSubset(fitStartPoint)


