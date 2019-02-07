        function x2=exploreparamFMS(gene, argsim)
        %% routine to optimise a certain gene with all the argsim parameter structure.
        if~isfield(argsim, 'numstarts')
            argsim.numstarts=100;
        end
        
        if~isfield(argsim, 'numtop')
            argsim.numtop=5;
        end
        if~isfield(argsim, 'numsteps')
            argsim.numsteps=30;
        end
        numtop=argsim.numtop;
        numstarts=argsim.numstarts;
        numsteps=argsim.numsteps;
        if ~isfield(argsim, 'simfunction')
        argsim.simfunction= @makesimulator4;
        end
        if ~isfield(argsim, 'model')
        argsim.model= 'mechModel3';
        end
        if ~isfield(argsim, 'opts') 
            argsim.opts= optimset('fminsearch')
            argsim.opts.MaxIter=5;
        end
        if ~isfield(argsim, 'opts2') 
            argsim.opts2= optimset('fminsearch')
            argsim.opts2.MaxIter=10;
        end
        if ~isfield(argsim, 'stepsize') 
        argsim.opts2.MaxIter=10;
        else
          argsim.opts2.MaxIter=argsim.stepsize;  
            
        end
       
        argsim.mf= extractModelFeatures(argsim.model);
        data=argsim.meandata.(gene).g1percent;
        argsim.data=data;
        datastart=nanmean(data(:,1));
        argsim.initialconditions= [datastart, 1,    0,   argsim.mth1ko(1), argsim.mig1ko(1), 1, 1];
        argsim.initialconditions(argsim.initialconditions<0)=0;
        argsim.initialconditions=argsim.initialconditions(1:numel(argsim.mf.varnames));
        %simulating and scoring the model with 100 startpoints
        startparams= repmat(argsim.defaultparams, numstarts, 1);
        randomvals=randn(numstarts, numel(argsim.mf.paramNames)); %creating random matrix of starting values
        startparams(:,argsim.onlyparams)= randomvals(:, argsim.onlyparams); %only modifying the params specified in onlyparams
        simulator=argsim.simfunction(argsim.model, argsim);
        
        l=[];
        disp('preliminary simulations');
        for n=1:size(startparams,1)
            disp(num2str(n))
            try
        [l(n),t,y,d]=simulator(startparams(n,:)); %plot(t,y, 'DisplayName',num2str(n)); hold on;
            catch
                l(n)=Inf;
            end
        end
        [a,b]= sort(l);
        startparams=startparams(b(1:argsim.numtop), :) %best params
        
        %startparams(:, 15)= rand(10, 1).*randi(4, 10,1);
        %startparams(:, 1)= startparams(:, 1)*randi(20);
        
        for k=1:size(startparams, 1)  
            %lastk.(genes{j}).(models{m})=k;
            try
                disp('trying fminsearch with new random param set')
        [x1(k, :), arr(k)]=fminsearch(simulator, startparams(k, :), argsim.opts); %we run fit twice just in case we run out of steps.
            catch
                arr(k)=Inf;
            end
        end
        fval.(gene).(argsim.model)=arr;
        trialparams.(gene).(argsim.model)=x1;
        [a,b]=min(arr); %looking for the smallest cost in all the simulations. 
        x2=x1(b,:);
        
        for o =1:argsim.numsteps
                tic;
            x2=fminsearch(simulator, x2, argsim.opts2);
            elapsed=toc;
            %refitparams2.(genes{j}).(models{m}).fminsearch=x2;
            if elapsed> 300 %if the duration of this bit is more than 5 mins then we reduce the number of steps by half
                argsim.numsteps= round(argsim.numsteps/2);
            end
            end
            end
               % refitparams2.(genes{j}).(models{m}).fminsearch=x2;
        %[l,t,y,d]=simulator(x2); plot(t,y, 'DisplayName', [genes{j} ' ' models{m}]); hold on;