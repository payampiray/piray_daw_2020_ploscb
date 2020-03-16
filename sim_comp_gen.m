function [me, se, CC, nans, dxstat]=sim_comp_gen

[me(1,:), se(1,:), CC(1,:), nans(1,:), dxstat(1,:)] = run('hgf',1);
[me(2,:), se(2,:), CC(2,:), nans(2,:), dxstat(2,:)] = run('vkf',1);

[me(3,:), se(3,:), CC(3,:), nans(3,:), dxstat(3,:)] = run('hgf',2);
[me(4,:), se(4,:), CC(4,:), nans(4,:), dxstat(4,:)] = run('vkf',2);

end

function [mE, sE, CC, nans, dxstat] = run(modelname,simset,ii)
simcat = 'comp_gen';
np   = 10000;
nsim = 1000;
nchunk = 50;
if nargin<3, ii = 1:nchunk; end

if strcmp(modelname,'hgf')
    switch simset
        case 1
            nu = .5;
            kappa = 1;
            omega = -3;
            sigma = 1;
        case 2
            nu = .25;
            kappa = 1;
            omega = -3;
            sigma = 1;
    end
    tx_sim = struct('sigma',sigma,'nu',nu,'kappa',kappa,'omega',omega);    
    rbpf_model = @(y)rbpf_hgf_lin(y,nu,kappa,omega,sigma,np);
end
if strcmp(modelname,'vkf')
    switch simset
        case 1
            lambda = .15;
            v0 = 1;
            sigma = 1;
        case 2
            lambda = .1;
            v0 = 1;
            sigma = 1;
    end
    tx_sim = struct('lambda',lambda,'v0',v0,'sigma',sigma);
    rbpf_model = @(y)rbpf_vkf_lin(y,lambda,v0,sigma,np);
end


% ------------
pipedir = getdefaults('pipedir');
tempdir = getdefaults('tempdir');
fsim    = fullfile(pipedir,simcat, sprintf('data_%s_set%d.mat',modelname,simset) );
ffit0   = sprintf('rbpf_%s_set%d_np%d',modelname,simset,np);
ffit    = fullfile(pipedir,simcat,sprintf('%s.mat',ffit0));
frandnum= fullfile(pipedir,simcat,sprintf('randnum_nsim%d.mat',nsim));
makedir(fullfile(pipedir,simcat));

do_fit = 1;
if exist(ffit,'file')    
   fpf = load(ffit); pf = fpf.pf;
   do_fit = 0;
end

rng(0);
Ys = sim_gen(fsim,nsim,modelname,tx_sim);

if ~exist(frandnum,'file')
    rng(0);
    randnum = randperm(nsim);
    save(frandnum,'randnum');
else
    frandnum = load(frandnum);
    randnum = frandnum.randnum;
end

if do_fit
    makedir(fullfile(tempdir,simcat));
    for i = ii
        nn = (i-1)*(nsim/nchunk) + (1:(nsim/nchunk));
        ffit_i = fullfile(tempdir,simcat,sprintf('%s_%02d.mat',ffit0,i));
        if ~exist(ffit_i,'file')    
            for j = 1:length(nn)                
                n = nn(j);
                y = Ys(:,n);
                rng(randnum(n));
                [mpf,vpf,cpf] = rbpf_model(y);
                pf(j) = struct('m',mpf,'v',vpf,'c',cpf);
            end
            save(ffit_i,'pf');
            fprintf('sim-chunk %02d\n', i);
        end    
    end
    [ffits,ok_all] = getfileordered(fullfile(tempdir,simcat),sprintf('%s_%s.mat',ffit0,'%02d'),1:nchunk);
    if ok_all
        pf(1:nsim) = deal(struct('m',[],'v',[],'c',[]));
        for i=1:nchunk
            nn = (i-1)*(nsim/nchunk) + (1:(nsim/nchunk));
            pfi = load(ffits{i});
            pf(nn) = pfi.pf;       
        end
        save(ffit,'pf');
    end
end

if exist(ffit,'file')
    [mE, sE, CC, nans, dxstat] = comp_pf(modelname, fsim, pf);
end

end

function [mE, sE, CC, nans, dxstat] = comp_pf(modelname, fsim, pf)
data = load(fsim);
Ys = data.y;
Xs = data.x;
Zs = data.z;

switch modelname
    case 'hgf'
        sigma = data.tx_sim.sigma;
        nu    = data.tx_sim.nu;
        kappa = data.tx_sim.kappa;
        omega = data.tx_sim.omega;
        
        Vs = exp(kappa*Zs+omega);
        model = @(y)hgf_lin(y,sigma,nu,kappa,omega);
    case 'vkf'
        lambda = data.tx_sim.lambda;
        v0     = data.tx_sim.v0;
        sigma  = data.tx_sim.sigma;        
        
        Vs = Zs.^-1;
        model = @(y)vkf_lin(y,lambda,v0,sigma);
end

nsim = size(Ys,2);
bad_traj = zeros(1,nsim);
ef      = nan(nsim,2);
dxs      = nan(nsim,1);
for n=1:nsim
    y = Ys(:,n);        
    x = Xs(2:end,n);
    v = Vs(2:end,n);
    
    mpf = pf(n).m(2:end,:);
    vpf = pf(n).v(2:end,:);    
    
    [mm, temp,vm] = model(y);
    if strcmp(modelname,'hgf')        
        mm(1,:) = [];                
        bad_traj(n) = temp;
        
        vm(1,:) = [];
        vm = exp(kappa*vm+omega);        
        vpf = exp(kappa*vpf+omega);
    end
    
    [mi] = kalman(y,v,sigma);
    
    e1 = median(abs(mm-x));
    e2 = median(abs(mpf-x));
    ef(n,1)  = (e1/e2-1)*100;    
    
    e1 = median(abs(mm-x));
    e2 = median(abs(mi-x));
    ef(n,2)  = (e1/e2-1)*100;
    
    ef(n,3)  = fisher(corr(mm,mpf));
    ef(n,4)  = fisher(corr(vm,vpf,'type','spearman'));        
    
    dx = abs(diff(x)); 
    dxs(n,:) = median(dx);
end
dxs(bad_traj==1,:) = nan;

dxstat = nanmean(dxs);
mE = nanmean(ef);
sE = nanserr(ef);
CC = invfisher(mE(3:4)); 
mE(3:4) = [];
sE(3:4) = [];

nans = sum(bad_traj);

dxstat  = round(dxstat*100)/100;
mE = round(mE*100)/100;
sE = round(sE*100)/100;
CC = round(CC*100)/100;
end

function [m] = kalman(y,v,omega)
[nt]    = size(y,1);
nq      = size(v,2);
y       = repmat(y,1,nq);

m       = zeros(1,nq);
w       = v(1).*ones(1,nq);
sigma   = v;

dv      = nan(nt,nq);
for t  = 1:nt      
    dv(t,:)     = m;
    v = sigma(t);
        
    k           = (w+v)./(w+v + omega);
    delta       = y(t,:) - m;
    m           = m + k.*delta;
    w           = (1-k).*(w+v);
end
dv(t+1,:) = m;

m = dv(2:end);
end

%--------------------------
function [y, x, z] = sim_gen(fsim,nsim,modelname,tx_sim)

switch modelname
    case 'hgf'
        gen_model = @gen_hgf_lin;                
    case 'vkf'
        gen_model = @gen_vkf_lin;                    
    otherwise
        error('!');
end

if ~exist(fsim,'file')
    for n=1:nsim
        [y(:,n),x(:,n),z(:,n)] = gen_model(tx_sim); %#ok<AGROW>
    end
    save(fsim,'y','x','z','tx_sim');    
end

data = load(fsim); 
y = data.y;
x = data.x;
z = data.z;

end

function [y, x, z] = gen_hgf_lin(tx)
N = 200;

sigma = tx.sigma;
nu    = tx.nu;
kappa = tx.kappa;
omega = tx.omega;

ok = 0;
while ~ok
    x2   = nan(N+1,1);
    x3 = nan(N+1,1);

    x2(1,:)= 0;
    x3(1,:) = 1;

    V   = nan(N+1,1);
    V(1,:) = exp(kappa*x3(1)+omega);
    y   = nan(N+1,1);

    for n=2:(N+1)
        % generate top-level
        x3(n) = randnormal(x3(n-1),nu);

        % generate low-level
        V(n) = exp(kappa*x3(n)+omega);
        x2(n) = randnormal(x2(n-1),V(n));

        y(n) = randnormal(x2(n),sigma);
    end    
    ok = 1;
    vars = [x3 x2 y]; vars(1,3) = 0;
    if any(isnan(vars(:)) | isinf(vars(:)))
        ok = 0;
    end    
end
z = x3;
x = x2;
y(1,:) = [];

y = y(1:N/2,:);
x = x(1:N/2+1,:);
z = z(1:N/2+1,:);
end

function [y, x, z] = gen_vkf_lin(tx)
N = 200;

lambda = tx.lambda;
v0     = tx.v0;
sigma  = tx.sigma;
        
eta = 1-lambda;
nu  = .5/(1-eta);

x   = nan(N+1,1);    
z   = nan(N+1,1);
y   = nan(N+1,1);

x(1,:) = 0;    
z(1,:) = v0^-1;
for n=2:(N+1)    
    epsil = betarnd(eta*nu,(1-eta)*nu) + eps;
    z(n) = z(n-1)*(eta.^-1)*epsil;    
    x(n) = randnormal(x(n-1),(z(n))^-1);
    y(n) = randnormal(x(n),sigma);
end
y(1,:) = [];

y = y(1:N/2,:);
x = x(1:N/2+1,:);
z = z(1:N/2+1,:);
end

function z = randnormal(m,v)
z = m + sqrt(v).*randn(size(m));
end

