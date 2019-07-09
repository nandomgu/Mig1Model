        function [x2, simulator]=exploreparmutants(gene, argsim, pars)
        %% routine to optimise a certain gene with all the argsim parameter structure.
        if~isfield(argsim, 'numstarts')
            argsim.numstarts=1000;
        end
        
        if~isfield(argsim, 'numtop')
            argsim.numtop=5;
        end
        if~isfield(argsim, 'numsteps') || isempty(argsim.opts)
            argsim.numsteps=10;
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

       % if ~isfield(argsim, 'opts2') || isempty(argsim.opts2)
            switch argsim.optimtype
                case {'fminsearch', 'fmincon'} 
                    argsim.opts= optimset(argsim.optimtype);
                    argsim.opts2= optimset(argsim.optimtype);
                otherwise
                    argsim.opts= optimoptions(argsim.optimtype);
                    argsim.opts2= optimoptions(argsim.optimtype);
            end
            argsim.opts.MaxIter=50;
            argsim.opts2.MaxIter=10;
            argsim.opts.Display='iter';
            argsim.opts2.Display='iter';
        %end
        if ~isfield(argsim, 'stepsize') 
        argsim.stepsize=200;
        end
          argsim.opts2.MaxIter=argsim.stepsize;  
            
        
       
        argsim.mf= extractModelFeatures(argsim.model);
        argsim.defaultparams=defaultparams(argsim.model);
        data=argsim.meandata.(gene).g1percent;
        argsim.data=data;
        datastart=nanmean(data(:,1));
        argsim.initialconditions= [datastart, 1,    0,   argsim.mth1ko(1), argsim.mig1ko(1), 1, 0, 0,0];
        argsim.initialconditions(argsim.initialconditions<0)=0;
        argsim.initialconditions=argsim.initialconditions(1:numel(argsim.mf.varnames));
        %simulating and scoring the model with 100 startpoints
        startparams= repmat(argsim.defaultparams, numstarts, 1);
        randomvals=randn(numstarts, numel(argsim.mf.paramNames)); %creating random matrix of starting values
        startparams(:,argsim.onlyparams)= randomvals(:, argsim.onlyparams); %only modifying the params specified in onlyparams
        simulator=argsim.simfunction(argsim.model, argsim);
        if nargin<3 || isempty(pars)
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
                switch argsim.optimtype
                    case 'fminsearch'
                        disp('optimisation algorithm: fminsearch')
        [x1(k, :), arr(k)]=fminsearch(simulator, startparams(k, :), argsim.opts); %we run fit twice just in case we run out of steps.
                    case 'simulannealbnd'
                        disp('optimisation algorithm: simmulannealbnd (simulated annealing)')
         [x1(k, :), arr(k)]=simulannealbnd(simulator, startparams(k, :), [], [], argsim.opts); 
                    otherwise
                        error('optimisation type not recognised');
                end
        
        catch
                arr(k)=Inf;
            end
        end
        fval.(gene).(argsim.model)=arr;
        trialparams.(gene).(argsim.model)=x1;
        [a,b]=min(arr); %looking for the smallest cost in all the simulations. 
        x2=x1(b,:);
        else
            
            x2=pars;
            
        end
        figure;
        for o =1:argsim.numsteps
                tic;
                argsim.opts2.MaxIter=argsim.stepsize;
                switch argsim.optimtype
                    case 'fminsearch'
            x2=fminsearch(simulator, x2, argsim.opts2)
                end
            elapsed=toc;
            [l,t,y,d]=simulator(x2); plot(t,y, 'DisplayName', [gene ' ' argsim.model]); hold on;
            if o==1
                plot(t, d(1:numel(t))); hold on;
            end
            xlim([0,20]);
            ylim([0,15]);
            pause(.5);
            %refitparams2.(genes{j}).(models{m}).fminsearch=x2;
            if elapsed> 300 %if the duration of this bit is more than 5 mins then we reduce the number of steps by half
                argsim.numsteps= round(argsim.numsteps/2);
                
            end
        end
            
        
        
            end
               % refitparams2.(genes{j}).(models{m}).fminsearch=x2;
        %[l,t,y,d]=simulator(x2); plot(t,y, 'DisplayName', [genes{j} ' ' models{m}]); hold on;