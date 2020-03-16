function fig6_tapas
addpath(fullfile(pwd,'tapas_HGF'));

simcat = 'basic';
pipedir = getdefaults('pipedir');
fname = fullfile(pipedir,simcat,'bin.mat');
data = load(fname);
o = data.o';
x = data.x';

[tx_vkf, tx_hgf] = sim_lr_vol;

lambda = tx_vkf(1);
v0 = tx_vkf(2);
omega = tx_vkf(3);

[m, k, v] = vkf_bin(o,lambda,v0,omega);
val1 = m;
vol1 = v;
lr1 = k;

nu = tx_hgf(1);
kappa = tx_hgf(2);
omega = tx_hgf(3);
[~, ~, mu2, mu3, sigma2] = hgf_tapas_bin(o,nu,kappa,omega);    
m = mu2(1:end-1);
v = (mu3(1:end-1));
sigma2 = sigma2(2:end);    
val2 = m;
vol2 = v;
lr2 = sigma2;

%------------------------------------------------------
close all;
fig_plot(x,o,val1,vol1,lr1,val2,vol2,lr2);
end

function fig_plot(x,y,val1,vol1,lr1,val2,vol2,lr2)
nr = 3;
nc = 2;
fpos0 = [0.2    0.0800    .55*1.0000    .65*0.8133];


fn = getdefaults('fn');
fnt = getdefaults('fnt');
fst = getdefaults('fst');
% fsl = getdefaults('fsl');
fsA = getdefaults('fsA');
xsA = getdefaults('xsA');
ysA = getdefaults('ysA');
abc = getdefaults('abc');
yst = getdefaults('yst');

%---------
figure;
set(gcf,'units','normalized');
set(gcf,'position',fpos0);

yl = [0 1.5; .2 1.2;-2.5 2.5];


sub_plts = [1 3 5];
hl = sim_C_plot(nr,nc,sub_plts,x,y,vol1,lr1,val1,yl);

text(.5,yst,'VKF','fontsize',fst,'Unit','normalized','fontname',fnt,'parent',hl(1),'HorizontalAlignment','Center','fontweight','bold');


for i=1:2
    text(xsA,ysA,abc(i),'fontsize',fsA,'Unit','normalized','fontname',fn,'parent',hl(i));
end
%---------
sub_plts = [2 4 6];

hr = sim_C_plot(nr,nc,sub_plts,x,y,vol2,lr2,val2,yl);
text(.5,yst,'HGF','fontsize',fst,'Unit','normalized','fontname',fnt,'parent',hr(1),'HorizontalAlignment','Center','fontweight','bold');

for i=1:2
    text(xsA,ysA,abc(i+2),'fontsize',fsA,'Unit','normalized','fontname',fn,'parent',hr(i));
end
end
