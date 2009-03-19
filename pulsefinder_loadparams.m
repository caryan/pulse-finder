function params = pulsefinder_loadparams(paramsfile)

%Helper function to load the params file and to do some error checking
%implement some default values

params = struct();

%Load it and evaluate line by line because eval(file) does not work in compiled code
paramsFID = fopen([paramsfile '.m'],'r');
if(paramsFID ~= -1)
    paramslines = textscan(paramsFID,'%s','delimiter','\n');
    paramslines = paramslines{1};
    for linect = 1:1:size(paramslines,1)
        eval(paramslines{linect});
    end
else
    error('paramsfile does not exist.');
end
fclose(paramsFID);

%Check for a molecule file and that the file exists
if(isfield(params,'nucleifile'))
    if(~exist(params.nucleifile,'file'))
        error('nuclei file specified does not exist.');
    end
else
    error('nucleifile field in params not specified');
end

%Check for the plength field
if(~isfield(params,'plength'))
    warning('No plength field in params.  Using default of 100');
    params.timesteps = 100;
end

%Check for timestep length
if(~isfield(params,'timestep'))
    warning('No timestep field in params.  Using default of 2us');
    params.timestep = 2e-6;
end

%Check for stepsize 
if(~isfield(params,'stepsize'))
    warning('No stepsize field in params.  Using default of 1e-1');
    params.stepsize = 1e-1;
end

%Check for subsystems
if(~isfield(params,'subsystem'))
    warning('No subsystem field in params.  Using default of all the spins.');
    params.subsystem{1} = [];
end

if(~isfield(params,'subsys_weight'))
    warning('No subsys_weight field in params.  Using default of equal weight');
    params.subsys_weight = ones(1,length(params.subsystem));
end







        
    
