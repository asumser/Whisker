function [fig]=plot_single_cell_example(trigs_e,trig_names,psth_e,spikes,trace_e,sem_trace_e,trace_name,ppms,bins,Name)
RastN=100;
pbins=bins(1:end-1)+diff(bins);
for s=1:numel(trigs_e)
    [SR{s}]=get_raster_single(trigs_e{s}, spikes,ppms,-bins(1),bins(end));
    isr_resF{s}=smoothdata(psth_e{s},2,'gaussian',milliseconds(3),'SamplePoints',milliseconds(pbins));
end
m=max(cat(1,isr_resF{:}),[],'all');

tt=tiledlayout(5,2);
title(tt,Name)
for s=1:2
    nexttile(1+(s-1),[2 1]);
    scatter(SR{s}(SR{s}(:,2)<=RastN,1),SR{s}(SR{s}(:,2)<=RastN,2),10,'k.');box off;
    xlim(bins([1 end]));ylim([1 RastN]);
    ylabel('trials')
    title(trig_names{s})

    nexttile(5+(s-1));
    yyaxis left
   if ~any(isnan(sem_trace_e(1,:,s)))
        semplot = cat(2,trace_e(1,:,s) + sem_trace_e(1,:,s),fliplr(trace_e(1,:,s) - sem_trace_e(1,:,s)));
        fill(cat(2,pbins,fliplr(pbins)), semplot,[.7 .7 .7]);
        hold on;
    end
     plot(pbins,trace_e(1,:,s)); hold off;
    ylabel(trace_name{1})
    yyaxis right
    if ~any(isnan(sem_trace_e(2,:,s)))
        semplot = cat(2,trace_e(2,:,s) + sem_trace_e(2,:,s),fliplr(trace_e(2,:,s) - sem_trace_e(2,:,s)));
        fill(cat(2,pbins,fliplr(pbins)), semplot,[.7 .7 .7]);
        hold on;
    end
     plot(pbins,trace_e(2,:,s)); hold off;
    ylabel(trace_name{2})
    xlim(pbins([1 end]));

    nexttile(7+(s-1),[2 1]);
    ax(s)=bar(pbins,isr_resF{s},1,'FaceColor',ones(1,3)/2);box off;
    xlim(pbins([1 end]));
    ylabel('Rate (Hz)')
    xlabel('Time (ms)')
    ylim([0 1.2*m])
end