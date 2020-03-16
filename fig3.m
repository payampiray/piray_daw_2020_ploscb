function fig3
[y1, x1] = timeseries_lin;
[y2, x2] = timeseries_bin;


lambda = .1;
v0 = .1;
omega = .1;

[m1, k1, v1] = vkf_lin(y1,lambda,v0,omega);
val(:,1) = m1;
vol(:,1) = v1;
kal(:,1) = k1;

y(:,1) = y1;
x(:,1) = x1;

% -------------------
[m1, k1, v1] = vkf_bin(y2,lambda,v0,omega);
m1 = 1./(1+exp(-m1));
val(:,2) = m1;
vol(:,2) = v1;
kal(:,2) = k1;

y(:,2) = y2;
x(:,2) = x2;

fig_plot(x,vol,kal,val,y);
end

function fig_plot(x,v,lr,m,y)

fpos0 = [0.2    0.0800    .55*1.0000    .7*0.8133];


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

%------------
nr = 3;
nc = 2;
sub_plts = [1 3 5];

yl = [0 .6;.4 .9;-1.4 1.4];

j = 1;
h = sim_A_plot(nr,nc,sub_plts,x(:,j),y(:,j)+nan,v(:,j),lr(:,j),m(:,j),yl);

for i=1:3
    text(xsA,ysA,abc(i),'fontsize',fsA,'Unit','normalized','fontname',fn,'parent',h(i));
%     set(h(i),'ylim',yl(i,:));
end


i =1;
text(.5,yst,'Linear','fontsize',fst,'Unit','normalized','fontname',fnt,'parent',h(i),'HorizontalAlignment','Center','fontweight','bold');

%------------
sub_plts = [2 4 6];

yl = [.08 .16;.4 .5;0 1.2];

j = 2;
h = sim_A_plot(nr,nc,sub_plts,x(:,j),y(:,j)+nan,v(:,j),lr(:,j),m(:,j),yl);

for i=1:3
    text(xsA,ysA,abc(i+3),'fontsize',fsA,'Unit','normalized','fontname',fn,'parent',h(i));
%     set(h(i),'ylim',yl(i,:));
end

i =1;
text(.5,yst,'Binary','fontsize',fst,'Unit','normalized','fontname',fnt,'parent',h(i),'HorizontalAlignment','Center','fontweight','bold');

end

function [o,x]=timeseries_lin
simcat = 'basic';
n = 10;
w = .01;

fname = sprintf('lin_n%d.mat',n);

p = [ones(1,2) repmat([-1 1],1,2) -ones(1,4)];
p = [p -p p -p];

pipedir = getdefaults('pipedir');
fdir = fullfile(pipedir,simcat); makedir(fdir);
fname = fullfile(fdir,fname);
if ~exist(fname,'file')
    nb = length(p);
    N  =  nb*n;
    x  = nan(1,N);
    o  = nan(1,N);
    for i=1:nb
        ii = (i-1)*n + (1:n);
        x(ii)  = p(i);
        o(ii)  = x(ii) + sqrt(w)*randn(1,n);
    end
    
    save(fname,'o','x');
end
data = load(fname);
o = data.o';
x = data.x';
end

function [o, x] = timeseries_bin
simcat = 'basic';
n = 10;
b = .8;
fname = sprintf('bin_n%d.mat',n);

p = [b*ones(1,2) repmat([1-b b],1,2) (1-b)*ones(1,4)];
p = [1-p p 1-p p];

pipedir = getdefaults('pipedir');
fdir = fullfile(pipedir,simcat); makedir(fdir);
fname = fullfile(fdir,fname);
if ~exist(fname,'file')        
    nb = length(p);  
    N  =  nb*n;    
    x  = nan(1,N);    
    o  = zeros(1,N);

    t0 = 0;
    for i=1:nb
        ii = t0 + (1:n);
        x(ii) = p(i);
        ni = randperm(n);
        ni = ni(1: round(p(i)*n));
        o(ii(ni))  = 1;
        t0 = t0 + n;
    end
    
    save(fname,'o','x');
end
data = load(fname);
o = data.o';
x = data.x';

end
