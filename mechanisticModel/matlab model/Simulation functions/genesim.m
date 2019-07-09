        function [l, y]=genesim(gene, argsim, pars)
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
        if nargin>2 &&  ~isempty(pars)
        argsim.defaultparams=defaultparams(argsim.model, paramstruct(argsim.model, pars));
%         
        end
         %if nargin>2 &&  ~isempty(pars) && isfield(argsim, 'onlyparams') && ~isempty(argsim.onlyparams)
%             argsim.defaultparams(argsim.onlyparams)=pars(argsim.onlyparams);
%         end
        data=argsim.meandata.(gene).g1percent;
        argsim.data=data;
        datastart=nanmean(data(:,1));
        if ~isfield(argsim, 'initialconditions') || isempty(argsim.initialconditions)
        argsim.initialconditions= [datastart, 1,    0,   argsim.mth1ko(1), argsim.mig1ko(1), 0, 0, 0,0];
        argsim.initialconditions(argsim.initialconditions<0)=0;
        argsim.initialconditions=argsim.initialconditions(1:numel(argsim.mf.varnames));
        end
        %simulating and scoring the model with 100 startpoints
        %startparams= repmat(argsim.defaultparams, numstarts, 1);
        %randomvals=randn(numstarts, numel(argsim.mf.paramNames)); %creating random matrix of starting values
        %startparams(:,argsim.onlyparams)= randomvals(:, argsim.onlyparams); %only modifying the params specified in onlyparams
        simulator=argsim.simfunction(argsim.model, argsim);
        
            
            [l,t,y,d]=simulator(pars);% plot(t,y, 'DisplayName', [gene ' ' argsim.model]); hold on;
            if ~isfield(argsim, 'plotdata') || argsim.plotdata==1
                plot(t, d); hold on;
            end
        %[l,t,y,d]=simulator(x2); plot(t,y, 'DisplayName', [genes{j} ' ' models{m}]); hold on;