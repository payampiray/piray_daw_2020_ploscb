function [e, c] = fig5
simcat = 'basic';
pipedir = getdefaults('pipedir');
fname = fullfile(pipedir,simcat,'bin.mat');
data = load(fname);
o = data.o';
x = data.x';

lambda = .1;
v0 = .1;
omega = .1;
nparticles = 10000;
rng(0);

[m, ~, v] = vkf_bin(o,lambda,v0,omega);
m = m(2:end);
val1 = 1./(1+exp(-m));
vol1 = v(2:end);

[m,v] = pf_vkf_bin(o,lambda,v0,nparticles);
m = m(2:end);
v = v(2:end);
val2 = 1./(1+exp(-m));
vol2 = v;

x(end) = [];
e1 = mean(abs(x - val1));
e2 = mean(abs(x - val2));

e = ((e1./e2) - 1)*100;

c(1) = corr(val1,val2,'type','spearman');
c(2) = corr(vol1,vol2,'type','spearman');
%------------------------------------------------------
close all;
fig_plot(x,o,val1,vol1,val2,vol2);
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

yl = [0.06 .22; -.2 1.2];

sub_plts = [1 3];
hl = sim_B_plot(nr,nc,sub_plts,x,y,vol1,val1,yl);

text(.5,yst,'VKF','fontsize',fst,'Unit','normalized','fontname',fnt,'parent',hl(1),'HorizontalAlignment','Center','fontweight','bold');


for i=1:2
    text(xsA,ysA,abc(i),'fontsize',fsA,'Unit','normalized','fontname',fn,'parent',hl(i));
end
%---------
yl = [0 .34; -.2 1.2];
sub_plts = [2 4];

hr = sim_B_plot(nr,nc,sub_plts,x,y,vol2,val2,yl);

ytick = [0 .1 .2 .3];
set(hr(1),'ytick',ytick);

text(.5,yst,'Benchmark','fontsize',fst,'Unit','normalized','fontname',fnt,'parent',hr(1),'HorizontalAlignment','Center','fontweight','bold');

hlg = legend(hr(2),{'Predicted','True'},'fontsize',fsl,'location','east');
pos = get(hlg,'Position');
pos(2) = 1.2*pos(2);
set(hlg,'Position',pos);


for i=1:2
    text(xsA,ysA,abc(i+2),'fontsize',fsA,'Unit','normalized','fontname',fn,'parent',hr(i));
end
end
