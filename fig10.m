function fig10
fitdir    = 'fit_B1';
file_hbi  = 'vkf_hgf_rl_kf';
pipedir   = getdefaults('pipedir');
fname_hbi = fullfile(pipedir,fitdir,sprintf('hbi_%s.mat',file_hbi));
fname_lap = fullfile(pipedir,fitdir,'lap_vkf.mat');


k = 1;

pnames = {'\lambda','v_0','\omega','\beta'};
modelnames = {'VKF','HGF','RW','KF'};

cbm     = load(fname_hbi); cbm = cbm.cbm;


%-------
freq = cbm.output.model_frequency;
pxp  = cbm.output.protected_exceedance_prob;

x  = cbm.output.group_mean{k};
xh = cbm.output.group_mean{k}+ cbm.output.group_hierarchical_errorbar{k};
xl = cbm.output.group_mean{k}- cbm.output.group_hierarchical_errorbar{k};

xs = [xl; x; xh];
txs = nan(size(xs));
for i=1:size(xs,1)
    [out] = fit_fx_vkf(xs(i,:));
    lambda = out.lambda;
    v0 = out.v0;
    omega = out.omega;        
    beta = out.beta;
    
    txs(i,:) = [lambda v0 omega beta];
    
    
    lambda = out.lambda;
    v0 = out.v0;
    omega = out.omega;        
    beta = out.beta;
    
    txs(i,:) = [lambda v0 omega beta];
end


for i=1:length(pnames)
    if ~strcmp(pnames{i}(1),'$'), pnames{i} = sprintf('$%s',pnames{i}); end
    if ~strcmp(pnames{i}(end),'$'), pnames{i} = sprintf('%s$',pnames{i}); end
end

txl = txs(1,:);
tx  = txs(2,:);
txh = txs(3,:);

%-------
[fdir,filename] = fileparts(fname_lap);
flr = fullfile(fdir,sprintf('%s_lr.mat',filename));
lr = load(flr);

tt_sig = lr.tt;
lr_sig = lr.lr;


fitting_plot(pxp, freq, modelnames, tx, txl, txh, pnames, tt_sig, lr_sig);
end

function fitting_plot(pxp, freq, modelnames, tx, txl, txh, pnames, tt_sig, lr_sig)
%-------
nr = 1;
nc = 3;
sub_plots = {1,2,3};

fn = getdefaults('fn');
fs = getdefaults('fs');
fsy = getdefaults('fsy');
fnt = getdefaults('fnt');
fst = getdefaults('fst');
fsA = getdefaults('fsA');
xsA = getdefaults('xsA');
ysA = 1.05; getdefaults('ysA');
abc = getdefaults('abc');

xst = .5;
yst = 1.05;getdefaults('yst');

cmaphbi = getdefaults('cmaphbi');
%-------------

fsalpha = 18;
fsxt = 14; % xticklabel


bwp = .25 ;

fpos0 = [0    0.3800    .85*1.0000    .35*0.8133];
alf = .6;
%---------

figure;
set(gcf,'units','normalized');
set(gcf,'position',fpos0);

%---------------
i = 1;
h(i) = subplot(nr,nc,sub_plots{i});
pos = get(h(i),'position');
text(xst,yst,sprintf('Bayesian model comparison'),'fontsize',fst,'fontname',fnt,...
    'Unit','normalized','fontweight','bold','Parent',h(i),'HorizontalAlignment','center'); hold on;
text(xsA,ysA,abc(i),'fontsize',fsA,'Unit','normalized','fontname',fnt); hold on;
set(gca,'visible','off');

% pxp(pxp<.01) = 0.01;
xm = {pxp, freq};
xe = {pxp*0, freq*0};
ylabels = {'PXP','Model frequency'};

np = 2;
dl = .05;
x0 = pos(1);    
lp = 1/np* (pos(3) -(np-1)*dl);
    
for i=1:np
    pos1 = pos;
    pos1(1) = x0+(i-1)*(lp+dl);
    pos1(3) = lp;

    axes('Position',pos1);
    
    
    errorbar1xN(xm{i},xe{i},modelnames,0,cmaphbi); hold on;     

    alpha(gca,alf);
    set(gca,'fontsize',fs,'fontname',fn);
    
    ax = ancestor(gca, 'axes');
    xaxes = get(ax,'XAxis');
    set(xaxes,'Fontsize',fsxt);
    
    ylabel(ylabels{i},'fontsize',fsy);
    set(gca,'ylim',[0 1.05]);
end  
%---------------
i = 2;
h(i) = subplot(nr,nc,sub_plots{i});
pos = get(gca,'Position');

text(xst,yst,sprintf('Parameters of VKF'),'fontsize',fst,'fontname',fnt,...
    'Unit','normalized','fontweight','bold','Parent',h(i),'HorizontalAlignment','center'); hold on;
text(xsA,ysA,abc(i),'fontsize',fsA,'Unit','normalized','fontname',fnt); hold on;
set(gca,'visible','off');

np = length(tx);
dl = .03;
x0 = pos(1);    
lp = 1/np* (pos(3) -(np-1)*dl);
    
for i=1:np
    pos1 = pos;
    pos1(1) = x0+(i-1)*(lp+dl);
    pos1(3) = lp;

    h = axes('Position',pos1);    
    errorbarKxN(tx(i),[txl(i); txh(i)],pnames{i},'',cmaphbi,0,bwp); hold on;        
    alpha(gca,alf);
    set(gca,'fontsize',fs,'fontname',fn);
    
    ax = ancestor(gca, 'axes');
    xaxes = get(ax,'XAxis');
    set(xaxes,'Fontsize',fsalpha,'TickLabelInterpreter','latex','fontweight','bold');
    
    ytick = get(h,'ytick');
    ytick(end) = [];
    set(h,'ytick', ytick);
    
    if i==1
        ylabel('Group parameters','fontsize',fsy);
    end
end  

%---------------
i = 3;

h(i) = subplot(nr,nc,sub_plots{i});
text(xst,yst,sprintf('Changes in learning rate'),'fontsize',fst,'fontname',fnt,...
    'Unit','normalized','fontweight','bold','Parent',h(i),'HorizontalAlignment','center'); hold on;
text(xsA,ysA,abc(i),'fontsize',fsA,'Unit','normalized','fontname','Calibri'); hold on;

mmg = nanmean(lr_sig,1);
N = size(lr_sig,1);

for n=1:N
    plot(tt_sig,lr_sig(n,:),'color',.7*ones(1,3)); hold on;
end
plot(tt_sig,mmg,'color','r','linewidth',2); hold on;
set(gca,'fontsize',fs,'fontname',fn);
set(gca,'xtick',tt_sig);
ylabel('Learning rate','fontsize',fsy);
xlabel('Time (locked to change-point)','fontsize',fsy)
end
