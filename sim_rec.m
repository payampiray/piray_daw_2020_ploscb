function [tx_sim, mx] = sim_rec(inp)
if nargin<1
    inp = 1;
end
i = inp;

[tx_sim, fsimfit] = run(i);
qtx = sim_error(fsimfit);

mx = quantile(qtx,[.25 .5 .75]);

end

function [tx_sim, fsimfit] = run(ii)

simcat = 'rec';
simstr = 'bin';
nsim = 500;

lambda = 0.2;
v0 = 5;
omega = 1;
beta = 1;
tx_sim = [lambda v0 omega beta];

logit = @(x)log(x./(1-x));
mu = [logit(lambda) logit(v0/10) log(omega) log(beta)]';
tau = 4*ones(size(mu));

normx = {@(x)1./(1+exp(-x)), @(x)10./(1+exp(-x)), @exp, @exp};
pnames = {'\lambda', 'v_0','\omega','\beta'};

task = 'sim_vkf_bin';
model = 'vkf_binary';        
N = 50;
T = 120;
ncue = 4;  

pipedir  = getdefaults('pipedir');
fsimfit  = fullfile(pipedir,simcat,sprintf('fit_%s.mat',simstr));
if ~exist(fsimfit,'file')
    sim_rec_sim(ii,simcat,simstr,N,model,task,pnames,mu,tau,T,normx,ncue);
    if ii(end)==nsim
        sim_rec_wrap(nsim,simcat,simstr);
    end
end
end

function [qtx]=sim_error(fsimfit)

simfit = load(fsimfit);
config = simfit.config;
fit = simfit.fit;

normx = config.normx;
for i=1:length(normx)
    normx{i} = str2func(normx{i});
end

qtx = nan(length(fit),4);
for i= 1:length(fit)
    x = fit(i).hbi.output.parameters{1};

    tx = nan(size(x));
    for j=1:size(x,2)
        ftx = (normx{j});
        tx(:,j) = ftx(x(:,j));
    end
    
    qtxi = mean(tx);
    
    qtx(i,:) = qtxi;
end
end
