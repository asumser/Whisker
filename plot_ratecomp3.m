function plot_ratecomp3(Rates,indiv_sig,RateLabels,CellTypes,DeflectionTypeLabels,plotMeanStd,YlabelName)
%%
%figure
Rates_in=Rates;

maxRate=max(Rates(CellTypes,:),[],'all','omitnan');
if all(Rates(CellTypes,:)>0,'all')
minRate=-.1*maxRate;
else
minRate=min(Rates(CellTypes,:),[],'all','omitnan');
end


off=0;
for DT=1:numel(DeflectionTypeLabels)
    selRates=Rates(CellTypes,:,DT);
    selRatesSig=indiv_sig(CellTypes,DT);
    
    [p]=signrank(selRates(~any(isnan(selRates),2),1),selRates(~any(isnan(selRates),2),2));
    if sum(selRatesSig)>0
    plot(off+[1:2],selRates(selRatesSig,:)','k'); hold on
    end
    if any(~selRatesSig)
    plot(off+[1:2],selRates(~selRatesSig,:)','Color',[.75 .75 .75]); hold on
    end
    switch plotMeanStd
        case 'meanstd'
        errorbar(off+[.85 2.15],mean(selRates,1,'omitnan'),std(selRates,0,1,'omitnan'),'r')
        case 'mean'
            line(off+[.85 2.15],mean(selRates,1,'omitnan'),'Color','r')
    end
    
    errorbar(off+1.5,maxRate*1.05,.65,'horizontal','r')
    text(off+1.5,maxRate*1.05,['p=' num2str(p)],'HorizontalAlignment','center','VerticalAlignment','bottom')
    
    text(off+1.5,maxRate*.9,DeflectionTypeLabels{DT},'FontSize',14,'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','middle')
    
    off=off+2;
    
    %
end
%
off=0;
if numel(DeflectionTypeLabels)>1
    for c=1:2
        [p]=ranksum(Rates_in(CellTypes,c,1),Rates_in(CellTypes,c,2));
        errorbar(off+1.85,maxRate*1.15+(off)*maxRate*0.05,1,'horizontal','r')
        text(off+1.85,maxRate*1.15+(off)*maxRate*0.05,['p=' num2str(p)],'HorizontalAlignment','center','VerticalAlignment','bottom')
        off=off+1.3;
    end
end

ylim([minRate-5 maxRate*1.3])
xlim([.75 2*numel(DeflectionTypeLabels)+.25])
set(gca,'Xtick',1:2*numel(DeflectionTypeLabels))
%
set(gca,'Xticklabel',RateLabels(repmat(1:2,1,2)))
box off
ylabel(YlabelName)
%%

