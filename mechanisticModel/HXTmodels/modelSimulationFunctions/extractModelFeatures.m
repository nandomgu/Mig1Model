function out= extractModelFeatures(modelname, modelpath, nonnegative)
%%this function extracts all important features for a model simulation by
%%parsing the model function file. This, of course, implies that the model
%%file should be considerably clean so that the parsing can proceed
%%flawlessly. specifically:
%number of parameters is determined by searching lines with 'params('
%number of model variables is retrieved by searching lines with 'x(' .
%number of outputs is determined by the number of lines with 'dout('.

if nargin<2 || isempty(modelpath)
model=eval(strjoin({'@' modelname}, ''));
localdir= pwd;
modelpath=strjoin({localdir '/HXTmodels/loopless_simulation_models/'},'');
end

modelpath=strjoin({modelpath '/' modelname '.m'},'');


    [s1,  s3]=system(['grep  ''params('' ' modelpath  ' | wc -l']);
    
   numpar= str2double(s3);
   
    [s1,  s7]=system(['grep  ''params('' ' modelpath  ' | cut -d ''='' -f 1 ']);
    
   parnames= strsplit(s7);
   parnames=parnames(1:end-1);
    %retrieve number of initial conditions in the model
   [s1,  s4]=system(['grep  ''x('' ' modelpath  ' | wc -l']);
   numInit=str2double(s4) ;
   [s1,  s8]=system(['grep  ''x('' ' modelpath '|cut -d "=" -f1' ]);
       namecell=strsplit(s8, '\n');
    varnames=namecell(1:end-1);
    %retrieve which parameters in the model are Hill coefficients. to run
    %this make sure that all the models have all parameters clustered
    %together and their lines contain the string params.
   % [s1,  sc]=system(['grep  ''params('' ' fnames{j}  ' | grep -n ''hill'' ']);
    [s1,  s5]=system(['grep  ''params('' ' modelpath  ' | grep -n ''hill'' | cut -d '':'' -f1 ']);

    
    hillIndices= str2double(strtrim(strsplit(s5, '\n')))';
  
    hillIndices=hillIndices(find(~isnan(hillIndices)));
    
    
    [s1,  s6]=system(['grep  ''dout('' ' modelpath  ' | wc -l']);
    numOutputs=str2double(s6) ;
    if nargin>2 && nonnegative==1
    options.NonNegative= 1:numOutputs;
    else
        options.NonNegative=[];
    end
    
    options.RelTol= .0001;

    
    out.options=options;
    out.varnames=varnames;
    out.numInit=numInit;
    out.paramNum= 1:numpar;
    out.numOutputs=numOutputs;
    out.paramNames=parnames;
    out.hillIndices=hillIndices;
    
end
