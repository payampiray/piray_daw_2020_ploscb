function [tt,lr,mlr,m] = fit_A_post_analysis(data,fname,hbi)
model = @(params,data)fit_A_vkf(params,data);

cbm     = load(fname); cbm = cbm.cbm;
[fdir,filename] = fileparts(fname);
flr = fullfile(fdir,sprintf('%s_lr.mat',filename));

N = length(data); 
if hbi
    k = 1;
    fx      = cbm.output.parameters{k};
    
    gfx     = cbm.output.group_mean{k};
else
    fx      = cbm.output.parameters;
    gfx = mean(fx,1);
end

tt = -5:5;

m = nan(N,length(tt));
lr = nan(N,length(tt));
mlr = nan(N,2);
for n=1:N
    [~,out]=model(fx(n,:),data{n});
    outcome = out.outcome;    
    pseq = data{n}.pseq;
    [~, k, ~] = vkf_bin(outcome,out.lambda,out.v0,out.omega);
    m(n,:) = signal_locked_changepoint(tt,pseq,k);    

    
    ihalf = floor(length(tt)/2);
    mlr(n,:) = [mean(m(n,1:ihalf)) mean(m(n,(ihalf+1):end))];
    
    
    [~,out]=model(gfx,data{n});
    outcome = out.outcome;    
    pseq = data{n}.pseq;
    [~, k, ~] = vkf_bin(outcome,out.lambda,out.v0,out.omega);
    lr(n,:) = signal_locked_changepoint(tt,pseq,k);        
end    
save(flr,'tt','lr','mlr','m');

end

function m = signal_locked_changepoint(tt,pseq,signal)
nq = size(pseq,2);
mi = zeros(0,length(tt));
for q=1:nq
    u = pseq(:,q)/100;
    ix = u_changepoint(u);
    miq = x_changepoint(signal(:,q), ix, tt);                
    mi = [mi; miq];
end
m = mean(mi,1);
end

function change_point = u_changepoint(u)

dx = [0; diff(u)~=0];
ix = (find(dx));
ix = [ix; length(u)];

dix = diff(ix);
ix1 = ix(dix<10); % beginning of a change period
ix2 = ix(dix>10); ix2(1) = [];

change_point = ix1;
end

function m = x_changepoint(x, ix, tt)
K = size(ix,1);
ix = [ix; size(x,1)];
xt = nan(K,length(tt));
for i=1:K
    t = ix(i)+ tt;
    xt(i,:) = (x(t));
end
% m = mean(xt,1);
m = xt;
end