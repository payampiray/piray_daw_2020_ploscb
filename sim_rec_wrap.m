function sim_rec_wrap(nsim,simcat,simstr)

tempdir  = getdefaults('tempdir');

pipedir  = getdefaults('pipedir');
makedir(fullfile(pipedir,simcat));
fsimfit  = fullfile(pipedir,simcat,sprintf('fit_%s.mat',simstr));

bigpipedir  = getdefaults('bigpipedir');
makedir(fullfile(bigpipedir,simcat));
fsimfull = fullfile(bigpipedir,simcat,sprintf('fit_%s.mat',simstr));

exfit = exist(fsimfit,'file');

if ~exfit
    simcatdir = fullfile(tempdir,simcat,simstr);        
    ok = zeros(1,nsim);
    while ~all(ok)            
        try %#ok<TRYNC>
            for j=1:nsim
                simdir = fullfile(simcatdir,sprintf('sim%04d',j));
                ok(j) = exist(fullfile(simdir,'hbi.mat'),'file');
            end
        end
    end
    for j=1:nsim
        simdir = fullfile(simcatdir,sprintf('sim%04d',j ));
        fdata = fullfile(simdir,'data.mat');
        fdata = load(fdata);
        y = fdata.data';
        y = cell2mat(y);
        
        
        fhbi = fullfile(simdir,'hbi.mat');
        fhbi = load(fhbi); 
        flap = fullfile(simdir,'lap.mat');
        flap = load(flap);
        fsim = fullfile(simdir,'sim.mat');
        fsim = load(fsim);
        sim  = fsim.sim;
        ful(j) = struct('sim',sim,'hbi',fhbi.cbm,'lap',flap.cbm,'data',y);
        fit(j) = struct('sim',sim,'hbi',fhbi.cbm,'lap',flap.cbm);
    end
    fconfig = fullfile(simcatdir,'config.mat');    
    config = load(fconfig); config = config.config;
    
    save(fsimfit,'fit', 'config');
    save(fsimfull,'ful', 'config');
end
end
