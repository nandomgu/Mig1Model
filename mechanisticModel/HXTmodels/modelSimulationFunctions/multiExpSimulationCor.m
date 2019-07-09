function handle=multiExpSimulationCor( modelname,outvar, initialConditions, figs,scaleFactors, normalize, varargin)
%% this function generates an array of inputs and outputs for each cExperiment included, simulates the model and calculates the costs with the outputs.
%% the outputs are by default a) assumed to be the same protein, b)normalized to the same minimum and maximum, c) assumed to be in the same 5min interval imaging.


modelFeatures= extractModelFeatures(modelname); 

model=eval(strjoin({'@' modelname}, ''));
if nargin<4
    initialConditions= zeros(1, modelFeatures.numInit);
end

srs=multiInputPanelNoFig([1:180], varargin{:});


if isempty(scaleFactors)
    scaleFactors= ones(1, numel(varargin));
end



inputs= srs.cy5;

for j=1:numel(inputs)
  inputs{j}= inputs{j}* scaleFactors(j);
    
end


outputs=srs.mean;
globmax=[];
globmin=[];
%normalize the TS conserving the relative maxima among experiments
%normoutputs=(outputs-min(min(outputs)))./max(max((outputs-min(min(outputs)))));
% for j= 1:numel(outputs)
%     globmax= [globmax max(outputs{j})];
%     globmin= [globmin min(outputs{j})];
% end
% globmax=max(globmax);
% globmin=min(globmin);
% 
% for j=1:numel(outputs)
%     
%     outputs{j}= (outputs{j}-globmin)/(globmax-globalmin);
% end

%normalize the TS regardless of the relative maxima (shape is the most
%important factor

for j=1:numel(outputs)
    
    if ~isempty(normalize)
    outputs{j}= normalizeTS(outputs{j});
    end
end

handle=@expSim;

    function [corvar, y]= expSim(params)
    


lsqdiff=0;
corvar=0;
    for j= 1:numel(inputs)
        
        
        func=rampSimHandle2cor(model, inputs{j},initialConditions, modelFeatures.options, outvar, outputs{j});
        [x{j},y{j}]=func(params);
        
        if figs==1
        subplot(2, numel(inputs), j);
       
        %func=rampSimHandle2(model, mediaInput,initialConditions, modelFeatures.options, outvar, gfpOutput);
        
        plot(inputs{j}, 'Color', [0 0 1]);
       subplot(2, numel(inputs), j+numel(inputs));
       
           
       
       
       if isempty(normalize)
           plot(1:numel([outputs{j}]), [outputs{j}] , 'LineWidth', 2, 'Color', 'blue');
       else
          plot(1:numel([outputs{j}]), normalizeTS([outputs{j}]) , 'LineWidth', 2, 'Color', 'blue');
       end
       
           hold on;
           if 1/x{j}>0
               costColor=0;
           else
           costColor=  1-1/x{j};
           end
           
       if isempty(normalize)
        plot(1:numel([y{j}]), y{j}' );  
           
       else
          plot(1:numel([y{j}]), normalizeTS([y{j}]') );
       end
       
        hold on;
        
        end
        
    end
    
    
    
    
    
    
    
    
    for j=1:numel(x)
        
        corvar= corvar+ x{j};
        lsqdiff= lsqdiff+sum((y{j}'-outputs{j}).^2);
        
    end

    
    
   
    
    system(sprintf(['echo ''' modelname '\t' mat2str(initialConditions,2)  '\t' mat2str(exp(params),3) '\t' mat2str(lsqdiff,3) '\t' mat2str(corvar,3) '''>> ' reverseDate '_modelcosts.txt']));



    end




end








