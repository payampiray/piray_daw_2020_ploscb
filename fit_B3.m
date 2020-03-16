function M = fit_B3(inp)
% fit particle-filter version of VKF

if nargin<1, inp = 0; end
nn = inp(1);
fitcat = 'fit_B1';
mname = 'pfvkf';

bmc = inp==0;
%--------------
pipedir = getdefaults('pipedir');

sampling_range = [0 .2;1 10];
sampling_init = [.1 1];
fitting_d = 1;

%--------------
nsamples = 1000;
nparticles = 10000;

if ~bmc
    ndatapoints = 160*ones(size(data));    
    data      = load(fullfile(pipedir,'data','dataset_matt2_noneg.mat')); data = data.data;    
    fit_fit_pf(data, nn, fitcat, mname, @pf_model, @response_model, sampling_range, sampling_init, fitting_d, ndatapoints, nsamples, nparticles);
end

% model comparison
if bmc
    mnames    = {'vkf','hgf','rl','kf','pfvkf'};
    M =length(mnames);
    lme = cell(1,M);
    for i = 1:length(mnames)
        fname = fullfile(pipedir, fitcat, sprintf('lap_%s.mat', mnames{i}) );
        cbm   = load(fname); 
        cbm = cbm.cbm;
        lme{i}  = cbm.output.log_evidence;     
    end   
    lme = cell2mat(lme);
    M = fit_NHI(lme,mnames);
end

end

function [loglik] = response_model(dv, params, data)
p        = 1./(1+exp(-dv));

choice   = data.choice;
Y        = choice==1;
[loglik] = fit_B_response(p,data.trialvalue,Y,params);
end

function [dv, vol, pf] = pf_model(nparticles, params, data)
lambda = params(1);
v0 = params(2);

image    = data.image;
outcome  = data.outcome;

nt = size(outcome,1);
dv = nan(nt,1);
vol= nan(nt,1);
for i=1:2
    tt = image==i;
    [x, v, pf(i)] = pf_vkf_bin(outcome(tt),lambda,v0,nparticles); %#ok<AGROW>
    dv(tt) = x(1:end-1);
    vol(tt) = v(1:end-1,:);
    
end

if any(isnan(dv)) || any(~isfinite(dv))
    dv = zeros(size(dv));
end
end
