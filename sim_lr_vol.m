function [tx_vkf, tx_hgf, mc, nc, nc_supp] = sim_lr_vol
[mc, nc, tx] = main;
nc_supp = supp;

tx_vkf = tx(2,:);
tx_hgf = tx(1,:);
end

function [mc, nc, tx] = main
models = {'hgf','vkf'};

for i=1:length(models)
    [ckv(:,i), tx(i,:)] = run(models{i}); %#ok<AGROW>
end

ckv = fisher(ckv);
mc = invfisher(mean(ckv));
nc = mean(ckv<0);
end

function [ckv, tx] = run(modelname)
simcat = 'comp_switch';
nsim = 500;


switch modelname
    case 'hgf'
        model = @model_hgf;
        d = 3;
    case 'vkf'
        model = @model_vkf;
        d = 3;
    otherwise
        error('!');
end

pipedir = getdefaults('pipedir');
fsim1    = fullfile(pipedir,simcat, sprintf('data_train.mat') );
fsim2    = fullfile(pipedir,simcat, sprintf('data_test.mat') );
ffit  = fullfile(pipedir,simcat, sprintf('fit_%s.mat',modelname) );
makedir(fullfile(pipedir,simcat));

rng(0);
Ys = sim_gen(fsim1,nsim);
sim_gen(fsim2,nsim);

do_fit = ~exist(ffit,'file');


vfit = 15.23;
if do_fit
    config = struct('range',[-5*ones(1,d);5*ones(1,d)],'numinit',7*d,'prior_for_bads',1);
    prior = struct('mean',zeros(d,1),'variance',vfit);
    cbm  = cbm_lap({Ys}, model, prior, [], config); %#ok<NASGU>
    save(ffit, 'cbm');    
end

cbm = load(ffit);
cbm = cbm.cbm;
[ckv, tx] = eval_model(fsim1, model, cbm);

end

function [ckv, tx] = eval_model(fsim, model, cbm)
data = load(fsim);
Ys = data.y;

nsim = size(Ys,2);

fx = cbm.output.parameters;

ckv = nan(nsim,1);
for n=1:nsim
    y = Ys(:,n);
    [~, tx, cc]= model(fx,y);                      
    ckv(n) = cc;        
end

end

function [loglik, tx, c] = model_hgf(params,data)
ux = @(x)(1./(1+exp(-x)));

y = data;
nu = ux(params(1));
kappa = ux(params(2));        
omega = params(3) -1;
tx = [nu kappa omega];

[~,bad_traj,m, mu3, sigma2] = hgf_bin(y,nu,kappa,omega);
m = 1./(1+exp(-m));

mu3 = mu3(1:end-1,:);
sigma2 = sigma2(2:end,:);
c = corr(mu3,sigma2,'type','spearman'); 
    
mu1 = m(1:end-1,:);
p = mu1.*y + (1-mu1).*(1-y);        
loglik = sum(sum(log(p+eps)));
if bad_traj>0
    loglik = -10^16;
end                
end

function [loglik, tx, c] = model_vkf(params,data)
ux = @(x)(1./(1+exp(-x)));

y = data;
lambda = ux(params(1));
v0 = 10*ux(params(2));
omega = exp(params(3));        
tx = [lambda, v0, omega];

[~,k,v,m] = vkf_bin(y,lambda,v0,omega);
m = [.5*ones(1,size(m,2)); m];

c = corr(v,k,'type','spearman');

mu = m(1:end-1,:);
p = mu.*y + (1-mu).*(1-y);
loglik = sum(sum(log(p+eps)));
end

function [Ys, Xs] = sim_gen(fsim,nsim)
if exist(fsim,'file')
    data = load(fsim); 
    Ys = data.y;
    Xs = data.x;    
    return;
end

n  = 20;
p0 = .8;

p  = [p0*ones(1,2) repmat([1-p0 p0],1,2) (1-p0)*ones(1,4)];
p  = [p p];

nb = length(p);  
N  = nb*n;

Xs = nan(N/2,nsim);
Ys = nan(N/2,nsim);
for j=1:nsim
    x  = nan(N,1);    
    y  = zeros(N,1);

    t0 = 0;
    for i=1:nb
        ii = t0 + (1:n);
        x(ii) = p(i);
        ni = randperm(n);
        ni = ni(1: round(p(i)*n));
        y(ii(ni))  = 1;
        t0 = t0 + n;
    end


    y = y(1:N/2,:);
    x = x(1:N/2,:);
    
    Xs(:,j) = x;
    Ys(:,j) = y;
end

end
%--------------------------------------------------------------------------

function nckv = supp
simcat = 'basic';
pipedir = getdefaults('pipedir');
fname = fullfile(pipedir,simcat,'bin.mat');
data = load(fname);
o = data.o';

rng(0);
nsim = 1000;
ckv = nan(nsim,2);

i = 1;
while i<=nsim
    lambda = rand;
    v0 = exp(randn);
    omega = exp(randn);

    [~, k, v] = vkf_bin(o,lambda,v0,omega);
    ckv(i,1) = corr(v,k,'type','spearman');

    nu = rand;
    kappa = rand;
    omega = randn - 3;
    
    [~, bad_traj, ~, mu3, sigma2] = hgf_bin(o,nu,kappa,omega);    
    v = mu3(1:end-1);
    sigma2 = sigma2(2:end);
    ckv(i,2) = corr(sigma2,v,'type','spearman');
    
    i = i+1;
    if bad_traj, i = i-1; end
end

nckv = mean(ckv<0);

end