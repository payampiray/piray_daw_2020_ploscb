function fit_fit_pf(data, nn, fitcat, mname, pf_model, response_model, sampling_range, sampling_init, fitting_d, ndatapoints, nsamples, nparticles)


%--------------
pipedir   = getdefaults('pipedir');
cmdir     = fullfile(pipedir,fitcat); makedir(cmdir);
tempdir   = getdefaults('tempdir');
lapdir    = fullfile(tempdir,fitcat,'laps'); makedir(lapdir);

if nargin<1
    nn = 1:size(data,1);
end

N  = length(data);
ffit = fullfile(cmdir,sprintf('fit_%s.mat',mname));
ffit_pf = fullfile(cmdir,sprintf('fit_%s_pf.mat',mname));
fname = fullfile(cmdir,sprintf('lap_%s.mat',mname));

if ~exist(ffit,'file')
    for n=nn        
        rng(0); % seed of random generator
        
        flap = fullfile(lapdir,sprintf('fit_%s_%04d.mat',mname,n));
        fpf = fullfile(lapdir,sprintf('fit_%s_%04d_pf.mat',mname,n));
        if ~exist(flap,'file')
            [cbm, pf] = fit(sampling_init, sampling_range, zeros(1,fitting_d), nsamples, nparticles, data(n), pf_model, response_model); %#ok<ASGLU>
            save(flap,'cbm')
            save(fpf,'pf')
        end

        [fnames,ok_laps] = getfileordered(lapdir,sprintf('fit_%s_%s.mat',mname,'%04d'),1:N);
        if ok_laps
            cbm = cbm_lap_aggregate(fnames); %#ok<NASGU>        
            save(ffit,'cbm');            
        end                
    end
end
if ~exist(fname,'file')
    [pfs,ok_pfs] = getfileordered(lapdir,sprintf('fit_%s_%s_pf.mat',mname,'%04d'),1:N);
    if ok_pfs
        pf_aggregate(pfs, ffit_pf, ffit, ndatapoints, fname);
    end
end
end

function [cbm, pf] = fit(pf_init, pf_range, resp_init, nsample, nparticles, data, pf_model, response_model)
v0 = 6.25;

d = length(resp_init);
M_config = struct('range',2*[-ones(1,d);ones(1,d)], 'numinit',1, 'verbose', 0);
M_prior = struct('mean',zeros(d,1),'variance',v0);

d = length(pf_init);
rng = pf_range(:,2)-pf_range(:,1);
mxQ = 0;

m_signals = cell(nsample,1);
v_signals = cell(nsample,1);
xx = nan(nsample,d);
lme = nan(1,nsample);
nverbose = round(nsample/10);

for n=1:nsample
    x = rand(d,1).*rng + pf_range(:,1);
    [xQ, vol] = pf_model(nparticles, x, data{1});
    m_signals{n} = xQ;
    v_signals{n} = vol;
    
    model = @(x,dat)response_model(xQ,x,dat);
    cbm  = cbm_lap(data, model, M_prior, [], M_config);    
    lme(n) = cbm.output.log_evidence;
    cbms(n) = cbm;
       
    xx(n,:) = x;
    mxQ = mxQ + xQ/nsample;
    
    if mod(n,nverbose)==0
        fprintf('%04d ', n);
    end    
end
fprintf('\n');
pf = struct('x_samples',xx,'m_signals',{m_signals}, 'v_signals',{v_signals},'lme',lme);

[~,i] = max(lme);
cbm = cbms(i);
end

function cbm = pf_aggregate(fnames, fpf, ffit, ndatapoints, fname)
cbm = load(ffit); cbm = cbm.cbm;
output = cbm.output;
log_evidence = output.log_evidence;

N = length(fnames);
for n=1:N
    pfn = load(fnames{n}); pfn = pfn.pf;
    xx = pfn.x_samples;
    lme = pfn.lme;
    [~,i] = max(lme);
    x = xx(i,:);

    x_pf(n,:) = x;
    log_evidence(n) = log_evidence(n) -.5*length(x)*log(ndatapoints(n));
    
    pf(n,1) = pfn;
end
output.log_evidence = log_evidence;
output.x_pf = x_pf;

cbm = struct('output',output);
save(fname,'cbm');

save(fpf,'pf');
end