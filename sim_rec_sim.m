function sim_rec_sim(ii,simcat,simstr,N,model,task,pnames,mu,tau,T,normx,ncue)

rng('shuffle');
for i=ii
    simrun(simcat,simstr,i,N,model,task,pnames,mu,tau,T,normx,ncue);
end

end

function simrun(simcat,simstr,id,N,model,task,pnames,mu,tau,T,normx,ncue) %#ok<INUSL>
simdir  = fullfile(getdefaults('tempdir'),simcat,simstr,sprintf('sim%04d',id));


model = str2func(model);
model = @(params, data)model(params, data, normx);
data = sim(simdir,T,N,mu,tau,model,pnames,normx,ncue);

d = length(mu);
run(simdir,data, model, d);

end

function [data] = sim(simdir,T,N,mu,tau,model,pnames,normx,ncue)
makedir(simdir);
fdatasim = fullfile(simdir,'data.mat');
if ~exist(fdatasim,'file')

    D    = length(mu);

    t = bsxfun(@times,((tau).^-.5),randn(D,N));
    h = bsxfun(@plus,t, + mu - mean(t,2));

    data = cell(N,1);
    for n=1:N
        [y, choice, r1, r2] = sim_vkf_bin(h(:,n), [T ncue], normx);
        data{n} = struct('outcome',y,'choice',choice,'rand1',r1,'rand2',r2);        
    end

    if ~exist(simdir,'dir')
        mkdir(simdir);
    end

    sim  = struct('h',h,'tau',tau,'mu',mu); %#ok<NASGU>
    
    fsim  = fullfile(simdir,'sim.mat'); 
    if ~exist(fsim,'file')
        save(fsim,'sim');
    end
    
    save(fdatasim,'data');
end

fconfig  = fullfile(simdir,'..','config.mat'); 
if ~exist(fconfig,'file')
    for i=1:length(normx)
        normx{i} = func2str(normx{i});
    end    
    config = struct('model',func2str(model),'pnames',{pnames},'normx',{normx}); %#ok<NASGU>    
    save(fconfig,'config');
end

fdata = load(fdatasim);
data = fdata.data;

end

function run(simdir,data, model, d)
v0 = 6.25;

flap = fullfile( simdir, sprintf('lap.mat') );
if ~exist(flap,'file')
    config1.numinit_med = 100;
    config1.numinit_up  = 100;
    config1.prior_for_bads = 1;
    
    config1.numinit = 10;

    prior = struct('mean',zeros(d,1),'variance',v0);
    cbm_lap(data, model, prior, flap, config1);    
end

fname_hbi = fullfile(simdir,sprintf('hbi.mat'));
flog_hbi = fullfile(simdir,sprintf('hbi.log'));
hbiconfig  = struct('flog',flog_hbi,'tolx',0.05);
if ~exist(fname_hbi,'file')
    cbm_hbi(data, {model}, {flap}, fname_hbi, hbiconfig,[]);
end

end

%--------------------------------------------------------------------------
function [y, choice, r1, r2] = sim_vkf_bin(params, siz, normx)
lambda   = normx{1}(params(1));
v0       = normx{2}(params(2));
omega    = normx{3}(params(3));
beta     = normx{4}(params(4));


%------
N = siz(1);
Q = siz(2);

r1 = rand(N,Q);
r2 = rand(N,Q);

m       = zeros(1,Q);
w       = omega*ones(1,Q);
v       = v0*ones(1,Q);

xQ      = nan(N,Q);

y       = nan(N,Q);
ux      = @(x)1./(1+exp(-x));
for t  = 1:N      
    xQ(t,:)     = m(1,:);

    wpre        = w;

    y(t,:)      = randbern(m,r1(t,:));

    delta       = sqrt(w+v).*(y(t,:) - ux(m));

    k           = (w+v)./(w+v + omega);
    m           = m + delta;        
    w           = (1-k).*(w+v);

    v           = v +lambda.*(delta.^2 + k.*wpre - k.*v);
end

choice = randbern(xQ*beta,r2);
choice(choice==0) = 2;
end

function loglik = vkf_binary(params, data, normx)
outcome = data.outcome;
choice  = data.choice;
choice(choice==2) = 0;

lambda   = normx{1}(params(1));
v0       = normx{2}(params(2));
omega    = normx{3}(params(3));
beta     = normx{4}(params(4));

m        = vkf_bin(outcome,lambda,v0,omega);

f = 1./(1+exp(-m*beta));
p = f.*choice + (1-f).*(1-choice);
loglik = sum(sum(log(p+eps)));

end

function y = randbern(mu,r)
p = 1./(1+exp(-mu));
y = nan(size(p));
y(r<=p) = 1;
y(r>p) = 0;

end
