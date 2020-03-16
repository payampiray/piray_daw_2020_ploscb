function fig7
simcat = 'basic';
pipedir = getdefaults('pipedir');
fname = fullfile(pipedir,simcat,'bin.mat');
data = load(fname);
o = data.o(1:80)';
x = data.x(1:80)';

tx = [.1 .5 .1;.5 .5 .1;.1 .2 .1;.1 .5 .2];

T   = length(o);
val = nan(T,size(tx,1));
vol = nan(T,size(tx,1));
for j=1:size(tx,1)
    [m1, ~, v1] = vkf_bin(o,tx(j,1),tx(j,2),tx(j,3));
    val(:,j) = m1;
    vol(:,j) = v1;
end
fig_plot(x,vol,val,tx);

end

function fig_plot(x,v,m,tx)

fpos0 = [0.2    0.0800    .85*1.0000    .65*0.8133];

fn  = getdefaults('fn');
fnt = getdefaults('fnt');
fst = getdefaults('fst');
fsA = getdefaults('fsA');
xsA = getdefaults('xsA');
ysA = getdefaults('ysA');
abc = getdefaults('abc');
%---------

figure;
set(gcf,'units','normalized');
set(gcf,'position',fpos0);

%------------
pnames = {'\lambda','v_0','\omega'};

yl = [0 .59;-3 3];

nr = 2;
nc = size(m,2);

for j = 1:size(m,2)
    sub_plts = [j j+size(m,2)];
    
    h(j,:) = sim_plot(nr,nc,sub_plts,x,v(:,j),m(:,j),yl); %#ok<AGROW>
    
    st = sprintf('Scenario %d',j);
    title(h(j),st,'fontsize',fst,'fontname',fnt );
    
    st = sprintf('$%s=%0.1f$, $%s=%0.1f$, $%s=%0.1f$',pnames{1},tx(j,1),pnames{2},tx(j,2),pnames{3},tx(j,3));
    text(.5,1.12, st,'fontsize',16,'Unit','normalized','fontname',fnt,'parent',h(j,2),...
         'HorizontalAlignment','center','VerticalAlignment','middle','fontweight','bold','Interpreter','latex');
    
    for i=1:2
    end     
    text(xsA,ysA,abc(j),'fontsize',fsA,'Unit','normalized','fontname',fn,'parent',h(j,1));
    
end

end

function h = sim_plot(nr,nc,sub_plts,x,v,m,yl)
fn = getdefaults('fn');
fsy = getdefaults('fsy');

cols = getdefaults('colmap');
lcol = [.8 .8 .8];
%---------

if isrow(x)
    dx = [0 diff(x)~=0];
else
    dx = [0; diff(x)~=0];
end
ix = (find(dx));
tt = 1:length(x);

%---------

h(1) = subplot(nr,nc,sub_plts(1)); 
plot(v,'color',cols(1,:),'linewidth',2); hold on;
if ~any(isnan(yl(1,:)))
    ylim(yl(1,:)); 
else
    ym = get(gca,'ylim');
    ym(2) = 1.01*ym(2);
    set(gca,'ylim',ym);
end
ym = get(gca,'ylim');
for i=1:length(ix)
    plot([tt(ix(i)); tt(ix(i))],ym','color',lcol,'linewidth',2);
end
plot(v,'color',cols(1,:),'linewidth',2); hold on;
set(gca,'fontname',fn);
ylabel('Volatility','fontsize',fsy);
set(gca,'ticklength', [0 0]);

h(2) = subplot(nr,nc,sub_plts(2)); 
plot(m,'color',cols(1,:),'linewidth',2); hold on;
if ~any(isnan(yl(2,:))), ylim(yl(2,:)); end
ym = get(gca,'ylim');
for i=1:length(ix)
    plot([tt(ix(i)); tt(ix(i))],ym','color',lcol,'linewidth',2);
end
plot(m,'color',cols(1,:),'linewidth',2); hold on;
set(gca,'fontname',fn);
ylabel('State predictions','fontsize',fsy); hold on;
set(gca,'ticklength', [0 0]);
xlabel('Trial','fontsize',fsy);
end