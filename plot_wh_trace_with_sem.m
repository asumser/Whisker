function plot_wh_trace_with_sem(trace_e_in,sem_trace_e_in,sel,selname,trace_name,relTime,Name)

tt=tiledlayout(1,2);
title(tt,Name)
for s=1:2
    nexttile(s);

    trace_e=cat(3,trace_e_in{sel{s}});
    trace_e=trace_e-mean(trace_e(:,relTime<0,:),2);
     trace_e=mean(trace_e,3,'omitnan');
    sem_trace=cat(3,sem_trace_e_in{sel{s}});
    sem_trace=mean(sem_trace,3,'omitnan');
     semplot = cat(2,trace_e + sem_trace,fliplr(trace_e - sem_trace));
    fill(cat(2,relTime,fliplr(relTime)), semplot,[.7 .7 .7],'EdgeColor','none');
        hold on;
   
     plot(relTime,trace_e); hold off;box off
    ylabel(trace_name)
    xlim(relTime([1 end]));
    title(selname{s})

end