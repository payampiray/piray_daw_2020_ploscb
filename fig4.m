function [e, c] = fig4
simcat = 'basic';
pipedir = getdefaults('pipedir');
fname = fullfile(pipedir,simcat,'lin.mat');
data = load(fname);
o = data.o';
x = data.x';


lambda = .1;
v0 = .1;
sigma = .1;
nparticles = 10000;
rng(0);

[m1, ~, v1] = vkf_lin(o,lambda,v0,sigma);

[m2,v2] = rbpf_vkf_lin(o,lambda,v0,sigma,nparticles);
m2(1) = [];
v2(1) = [];

e1 = mean(abs(x - m1));
e2 = mean(abs(x - m2));
e = (e1./e2 - 1)*100;

c(1) = corr(m1,m2);
c(2) = corr(v1,v2,'type','spearman');
%------------------------------------------------------
close all;
fig_plot(x,o,m1,v1,m2,v2);
end

function fig_plot(x,y,val1,vol1,val2,vol2)
nr = 2;
nc = 2;
fpos0 = [0.2    0.0800    .55*1.0000    .5*0.8133];

fn = getdefaults('fn');
fnt = getdefaults('fnt');
fst = getdefaults('fst');
fsl = getdefaults('fsl');
fsA = getdefaults('fsA');
xsA = getdefaults('xsA');
ysA = getdefaults('ysA');
abc = getdefaults('abc');
yst = getdefaults('yst');

%---------
figure;
set(gcf,'units','normalized');
set(gcf,'position',fpos0);

yl = [0 .39; -1.5 1.5];

sub_plts = [1 3];
hl = sim_B_plot(nr,nc,sub_plts,x,y,vol1,val1,yl);

text(.5,yst,'VKF','fontsize',fst,'Unit','normalized','fontname',fnt,'parent',hl(1),'HorizontalAlignment','Center','fontweight','bold');

for i=1:2
    text(xsA,ysA,abc(i),'fontsize',fsA,'Unit','normalized','fontname',fn,'parent',hl(i));
end

%---------
sub_plts = [2 4];
hr = sim_B_plot(nr,nc,sub_plts,x,y,vol2,val2,yl);

text(.5,yst,'Benchmark','fontsize',fst,'Unit','normalized','fontname',fnt,'parent',hr(1),'HorizontalAlignment','Center','fontweight','bold');

legend(hr(2),{'Predicted','True'},'fontsize',fsl,'location','northeast');

for i=1:2
    text(xsA,ysA,abc(i+2),'fontsize',fsA,'Unit','normalized','fontname',fn,'parent',hr(i));
end
end
