function hb = errorbar1xN(mx,ex,labels,basevalue,colmap)
if nargin<4, basevalue = 0; end
if nargin<5, colmap = [.9 .9 .9]; end

% N = length(mx);
barwidth = 0.25;

if ~isrow(mx), mx = mx'; end
if ~isrow(mx), error('mx should be a vector'); end


if size(ex,2)~=size(mx,2), ex=ex'; end
if size(ex,2)~=size(mx,2), error('ex is not matched with mx'); end

if size(ex,1)==1
    e    = [mx-ex; mx+ex];    
elseif size(ex,1)>1
%     e = [mx-ex(1,:); mx+ex(2,:)];
    e = [ex(1,:); ex(2,:)];
end

% figure;
hb = bar(mx,barwidth,'FaceColor',colmap(1,:),'EdgeColor','k','linewidth',1,'basevalue',basevalue);
hold on;
a =1:length(mx); %get(gca,'xtick');
for i=1:length(mx)    
    plot([a(i);a(i)],e(:,i),'-','color','k','linewidth',2);
end
if ~isempty(labels)
set(gca,'xticklabel',labels,'xcolor','k','fontsize',14);
end
set(gca,'xlim',[a(1)-.5 a(end)+.5]);
set(gca,'ticklength', [0 0]);

% % yrng = get(gca,'ytick');
% % ystp = diff(yrng); ystp = ystp(1);
% % 
% % sm = max(e(2,:))+0.5*ystp;
% % % if ~isempty(asterisk)    
% % %     for i=1:length(mx)    
% % %         if asterisk(i)
% % %             plot(a(i),sm,'*','color','k');
% % %         end
% % %     end    
% % % end
% % 
% % ymax = yrng(end)+.95*ystp;
% % ymin = yrng(2)-0.95*ystp;
% % % set(gca,'ylim',[ymin;ymax]);
% % % set(gca,'XTickLabelRotation',45);
% % set(gca,'box','off');
% % set(gca,'ticklength', [0 0]);
% % % % pos = get(gcf,'position');
% % % % pos(1:2) = [100 100];
% % % % pos(3:4) = rs*pos(3:4);
% % % % % set(gcf,'position',pos);    

end
