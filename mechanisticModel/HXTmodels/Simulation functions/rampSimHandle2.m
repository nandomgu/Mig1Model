function handle=rampSimHandle2(modelHandleProducer, inputRamp, initialConditions, options, outvar,data)
%%ramp simulation 2 receives handles to loopless models (ending in L)
%%which, instead of conventionally receiving just the list of parameters,
%%they also receive the media input. in the model handle producer function.
%%this should make computation easier but as otf 20151018 the glucose input
%%cannot be seen.


handle= @rampSim;
%%the output of this function is a handle for a simulation function. that
%%simulation function has the function parameters as an input and returns the simulated time series for the
%%mature hxt as well as the element-wise square difference from the data.
%% this simuation can be used on its own   to test and plot different parameter values manually or to insert it into an optimization function that minimizes
%%lsqdiff
    function [lsqdiff,y,t]=rampSim(params)
       %% these params are assumed to be, actually, logarithms of the params. we immediately
       %%bring these out so we can input them into our model. we do this
       %%because we don't want the fitting algorithms to make our
       %%parameters negative.
       
       params=exp(params);
        
        
        func=modelHandleProducer(params, inputRamp);
        
        
        
%         output=[];
%         stincond=initialConditions;
        
%         dbg=0;
%         for j=1:length(inputRamp)
            
            %disp(initialConditions)
%             initialConditions(1)= inputRamp(j);
            
            
            
            [t, y]=ode23(func, [0 numel(0:5:((numel(inputRamp)*5)-1))], initialConditions, options);
           %[t, y]=ode23(func, [0:numel(inputRamp)-1], initialConditions, options);
            
%             output(j, :)= y(end,:);
%             
% %             i
% %             if (dbg>5)
% %                 
% %             dbg='hello world';
% %             
% %             end
% %             dbg=dbg+1
% %             
%             
%             initialConditions=y(end,:);
%             
%         end
        %no scaling performed
        %y=output(:,outvar);
        y= y(:,outvar);
% y=interp1(t, y, 0:numel(data)-1);
 %y=normalizeTS(y);
        %the line below scales the sim between 0 and 1
        %y=(output(:,outvar)-min(output(:,outvar)))/max(output(:,outvar)-min(output(:,outvar)));
        %data=(data-min(data))/max(data);
        %y=(y-data).^2;
        %output=[];
       % initialConditions= stincond;
       %computing the sum of squares for fminsearch
       %lsqdiff=sum((y-interp1(0:5:(numel(data)*5)-1, data, t)).^2);
       %computing the element wise difference for lsqnonlin to work with. 
       if isempty(data)
      lsqdiff=NaN;
       else
            lsqdiff=(data-y';
       end
        %figure; plot(y)
        %initialConditions= stincond;
    end
 

end