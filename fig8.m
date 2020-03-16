function fig8

fitdir    = 'fit_A1';
file_hbi  = 'vkf_hgf_rl_kf';
pipedir   = getdefaults('pipedir');
fname_hbi = fullfile(pipedir,fitdir,sprintf('hbi_%s.mat',file_hbi));
k = 1;

pnames = {'\lambda','v_0','\omega','\beta'};

cbm     = load(fname_hbi); cbm = cbm.cbm;


%-------

x  = cbm.output.group_mean{k};
xh = cbm.output.group_mean{k}+ cbm.output.group_hierarchical_errorbar{k};
xl = cbm.output.group_mean{k}- cbm.output.group_hierarchical_errorbar{k};

xs = [xl; x; xh];
txs = nan(size(xs,1),length(pnames));
for i=1:size(xs,1)
    [out] = fit_fx_vkf(xs(i,:));    
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

%-------
fname = fullfile(pipedir,'fit_A1','example.mat');

example = load(fname,'x','y','fx');

[output] = fit_fx_vkf(example.fx);
lambda = output.lambda;
v0 = output.v0;
omega = output.omega;

m = vkf_bin(example.y,lambda,v0,omega);
m = 1./(1+exp(-m));

task_plot(m, example.x, example.y);
end

function task_plot(m, x, y)
%-------
nr = 1;
nc = 2;

fn = getdefaults('fn');
fsy = getdefaults('fsy');
fnt = getdefaults('fnt');
fsl = getdefaults('fsl');
fsA = getdefaults('fsA');
xsA = getdefaults('xsA');
ysA = 1.05; getdefaults('ysA');
abc = getdefaults('abc');

cols = getdefaults('colmap');

fpos0 = [0    0.0800    .65*1.0000    .35*0.8133];

%---------

figure;
set(gcf,'units','normalized');
set(gcf,'position',fpos0);
% set(gcf,'name',figname);

%---------
i = 1;
h(i) = subplot(nr,nc,1);
text(xsA,ysA,abc(i),'fontsize',fsA,'Unit','normalized','fontname',fnt); hold on;
set(gca,'visible','off');

%---------------
i = 2;
h(i) = subplot(nr,nc,2);
text(xsA,ysA,abc(i),'fontsize',fsA,'Unit','normalized','fontname',fnt); hold on;

plot(m,'color',cols(1,:),'linewidth',2); hold on;
plot(x,'color',cols(2,:),'linewidth',1); hold on;
if all(~isnan(y))
    plot(y,'.','color','k');
end
set(gca,'fontname',fn);
ylabel(sprintf('Probability\n of outcome being correct'),'fontsize',fsy); hold on;
set(gca,'ticklength', [0 0]);
xlabel('Trial','fontsize',fsy);
ylim([-.19 1.19]);

hlg = legend(h(i),{'VKF','True'},'fontsize',fsl,'location','east');
pos = get(hlg,'Position');
pos(2) = 1.3*pos(2);
set(hlg,'Position',pos);

end
