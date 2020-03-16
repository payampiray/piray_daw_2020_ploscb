function [tt,lr,mlr,m] = fit_B_post_analysis(data,fname,hbi)

model = @(params,data)fit_B_vkf(params,data);
cbm     = load(fname); cbm = cbm.cbm;
[fdir,filename] = fileparts(fname);
flr = fullfile(fdir,sprintf('%s_lr.mat',filename));


N = length(data);

if hbi
    k = 1;
    fx      = cbm.output.parameters{k};
    gfx      = cbm.output.group_mean{k};
else
    fx      = cbm.output.parameters;
    gfx = mean(fx);
end


tt = -5:5;
m = nan(N,length(tt));
lr = nan(N,length(tt));
mlr = nan(N,2);
for n=1:N
    [~,out]=model(fx(n,:),data{n});
    outcome = out.outcome;
    image = out.image;
    pseq = data{n}.rewardprob;
    
    for i=1:2
        o = outcome(image==i);
        [~, sig(:,i)] = vkf_bin(o,out.lambda,out.v0,out.omega); %#ok<AGROW>
    end
    m(n,:) = signal_locked_changepoint(tt,image,pseq,outcome,sig);
    
    ihalf = floor(length(tt)/2);
    mlr(n,:) = [mean(m(n,1:ihalf)) mean(m(n,(ihalf+1):end))];
    

    [~,out]=model(gfx,data{n});
    for i=1:2
        o = outcome(image==i);
        [~, gsig(:,i)] = vkf_bin(o,out.lambda,out.v0,out.omega);
    end
    lr(n,:) = signal_locked_changepoint(tt,image,pseq,outcome,gsig);        
end

ilr = nanmean(mlr,2);
ibad = false(size(ilr));

mlr(ibad,:) = [];
m(ibad,:) = [];

save(flr,'tt','lr','mlr','m');
end


function m = signal_locked_changepoint(tt,image,pseq,outcome,signal)
xt = zeros(0,length(tt));
for i=1:2
    o = outcome(image==i);
    p = pseq(image==i, i);
    [ix] = u_changepoint(p,o);   
    mi = x_changepoint(signal(:,i), ix, tt);
    xt = [xt ; mi]; %#ok<AGROW>
end
m = mean(xt,1);
end

function [ix, do] = u_changepoint(pseq, outcome)

dx = [0; diff(pseq)~=0];
ix = (find(dx));
ix = [ix; length(pseq)];

dix = diff(ix);
ix = ix(dix>10);

% too early to count as a change
ix(ix<=10) = [];

% make sure that outcome has actually changed
tt0 = -10:9;
thalf = length(tt0)/2;
mo = nan(length(ix),length(tt0));
for i = 1:length(ix)
    tt = ix(i) + tt0;
    mo(i,:) = outcome(tt);
end
mmo1 = mean(mo(:,1:thalf),2);
mmo2 = mean(mo(:,(thalf+1):end),2);
do = [mmo1 mmo2];
ok = (abs(mmo2-mmo1))>=.5;

ix = ix(ok);
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