function M = fit_A3(inp)
% fit particle-filter version of VKF

if nargin<1, inp = 0; end
nn = inp(1);
fitcat = 'fit_A1';
mname = 'pfvkf';

bmc = inp==0;
%--------------
pipedir   = getdefaults('pipedir');

sampling_range = [0 .2;1 10];
sampling_init = [.1 1];
fitting_d = 4;

%--------------
nsamples = 1000;
nparticles = 10000;


if ~bmc
    ndatapoints = 480*ones(size(data));    
    data      = load(fullfile(pipedir,'data','dataset1.mat')); data = data.data;    
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
choice   = data.choice;

if any(any(isnan(dv))) || any(any(~isfinite(dv)))
    dv = zeros(size(dv));
end

nt       = size(choice,1);
X        = dv(1:nt,:);
Y        = choice==1;
loglik   = fit_A_response(X,Y,params);

end

function [xQ, vol, pf] = pf_model(nparticles, params, data)
lambda = params(1);
v0 = params(2);

choice   = data.choice;
outcome  = data.outcome;

outcome(:,1:2) = 2*outcome(:,1:2) - 1;
outcome(:,3:4) = 2*outcome(:,3:4) + 1;
outcome(choice==2) = -outcome(choice==2);

outcome(outcome==-1) = 0;

N = size(outcome,1);
xQ = nan(N,4);
vol = nan(N,4);
for i=1:size(xQ,2)
    [x, v,pf(i)] = pf_vkf_bin(outcome(:,i),lambda,v0, nparticles); %#ok<AGROW>
    xQ(:,i) = x(1:end-1);
    vol(:,i) = v(1:end-1,:);
end
if any(any(isnan(xQ))) || any(any(~isfinite(xQ)))
    xQ = zeros(size(xQ));
end
end
