function h = sim_B_plot(nr,nc,sub_plts,x,y,v,m,yl)
fs = getdefaults('fs');
fn = getdefaults('fn');
fsy = getdefaults('fsy');
alf = getdefaults('alf');
fsA = getdefaults('fsA');
xsA = -.05 + getdefaults('xsA');
ysA = getdefaults('ysA');
ysA = 1.05;
abc = getdefaults('abc');

cols = getdefaults('colmap');
lcol = [.8 .8 .8];
%---------

if isrow(x)
    dx = [0 diff(x)~=0];
else
    dx = [0; diff(x)~=0];
end
ix = (find(dx));
ttvol = ix(1):(ix(end)+ix(2)-ix(1)-1);

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
% for t= ttvol(1):1:ttvol(end)
%     plot([t; t],ym','color',[.8 .8 .8],'linewidth',2); hold on;    
% end
for i=1:length(ix)
    plot([tt(ix(i)); tt(ix(i))],ym','color',lcol,'linewidth',2);
end
plot(v,'color',cols(1,:),'linewidth',2); hold on;
set(gca,'fontname',fn);
% legend({'Volatility estimate','True change-points'},'location','northeast','fontsize',fsl);
% ylabel('$v_t$','Interpreter','latex','fontsize',fsy);
ylabel('Volatility','fontsize',fsy);
set(gca,'ticklength', [0 0]);

h(2) = subplot(nr,nc,sub_plts(2)); 
plot(m,'color',cols(1,:),'linewidth',2); hold on;
plot(x,'color',cols(2,:),'linewidth',1); hold on;
if all(~isnan(y))
    plot(y,'.','color','k');
end
set(gca,'fontname',fn);
ylabel('Predictions','fontsize',fsy); hold on;
set(gca,'ticklength', [0 0]);
xlabel('Trial','fontsize',fsy);
if ~any(isnan(yl(2,:))), ylim(yl(2,:)); end
end