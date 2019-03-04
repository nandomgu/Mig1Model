function plotall(plotargs)
%plotting mechanistic simulations for all mutants and concentrations
yall= plotargs.y.(plotargs.modelname);
concolors=plotargs.concolors;
t=plotargs.t;
input=plotargs.input;
if ~isfield(plotargs, 'lw')
    lw=1.5;
else
    lw=plotargs.lw;
end

%%
strs={'WT', 'Mig1?', 'Mth1?','std1?','rgt2?','snf3'};
figure;
k=1;
c=1;
while k<18
axes(); 
j=0; 
plot(t, real(yall(:, k+j)), 'color', concolors(j+1,:), 'LineWidth', lw);
hold on;
j=j+1
plot(t,real(yall(:, k+j)), 'color', concolors(j+1,:), 'LineWidth', lw);
hold on;
j=j+1
plot(t,real(yall(:, k+j)), 'color', concolors(j+1,:), 'LineWidth', lw);
hold on;
title([plotargs.modelname ' ' strs{c}])
k=k+3;
c=c+1;
end
stackPlots(gcf)


        disp('plotting sugar')
        allaxes=get(gcf, 'children')
        for j=1:numel(allaxes)
      
            axes(allaxes(j))
            yyaxis right
        area( t, input, 'FaceAlpha', 0.1, 'FaceColor', [1, .6, 0], 'EdgeColor', 'none'); 
        ylim([0, 1])
        yyaxis left
        
        end
       

%%
%axes(); k=k+3; plot(real(y(:, k:(k+2)))); axes();k=k+3; plot(real(y(:, k:(k+2)))); axes();k=k+3; plot(real(y(:, k:(k+2)))); axes();k=k+3; plot(real(y(:, k:(k+2))));k=k+3;axes(); plot(real(y(:, k:(k+2)))); stackPlots(gcf)





