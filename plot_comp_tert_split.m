function plot_comp_tert_split(cr,vpmr,pomr,p,plim,ylab,split_var_name)

tt=tiledlayout(2,2);
title(tt,split_var_name,'1st vs 3rd tertile comparison')
nexttile;
selx=vpmr(p(vpmr,1)>=plim);
if ~isempty(selx)
    plot([1 3],cr(selx,[1 3],2,1),'Color',[.5 .5 .5]);hold on
end
selx=vpmr(p(vpmr,1)<plim);
if ~isempty(selx)
    plot([1 3],cr(selx,[1 3],2,1),'Color','k');hold off
end
xlim([.5 3.5]);box off;
title('vpm puff')
ylabel(ylab)
set(gca,'XTick',[1 3])
xlabel('tertile')

nexttile;
selx=pomr(p(pomr,1)>=plim);
if ~isempty(selx)
    plot([1 3],cr(selx,[1 3],2,1),'Color',[.5 .5 .5]);hold on
end
selx=pomr(p(pomr,1)<plim);
if ~isempty(selx)
    plot([1 3],cr(selx,[1 3],2,1),'Color','k');hold off
end
xlim([.5 3.5]);box off;
title('pom puff')
ylabel(ylab)
set(gca,'XTick',[1 3])
xlabel('tertile')

nexttile;
selx=vpmr(p(vpmr,2)>=plim);
if ~isempty(selx)
    plot([1 3],cr(selx,[1 3],2,2),'Color',[.5 .5 .5]);hold on
end
selx=vpmr(p(vpmr,2)<plim);
if ~isempty(selx)
    plot([1 3],cr(selx,[1 3],2,2),'Color','k');hold off
end
xlim([.5 3.5]);box off;
title('vpm touch')
ylabel(ylab)
set(gca,'XTick',[1 3])
xlabel('tertile')
nexttile;
selx=pomr(p(pomr,2)>=plim);
if ~isempty(selx)
    plot([1 3],cr(selx,[1 3],2,2),'Color',[.5 .5 .5]);hold on
end
selx=pomr(p(pomr,2)<plim);
if ~isempty(selx)
    plot([1 3],cr(selx,[1 3],2,2),'Color','k');hold off
end
xlim([.5 3.5]);box off;
title('pom touch')
ylabel(ylab)
xlabel('tertile')
set(gca,'XTick',[1 3])