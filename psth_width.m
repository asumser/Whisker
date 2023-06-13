parameter='Puff';

meth='gaussian';
contPlot={MeanWp,MeanP,MeanAp};
exclude={'Pole';'Light';'Exclude';'Grooming'};include=[];
Name='Puff';
[P_PSTH_instrate]=plot_Pop_PSTH_wrapper3(parameter, DiscreteData,puff_triglength_limit,ppms,psth_time_before_trig,psth_time_after_trig,exclude,psth_safety_margin,psth_binsize,[],Name,{vpmr;pomr},{'VPM';'POm'},include,'latency',-inf,contPlot,{'meanbefore','scaleminmax','mono'},'baseline corr. whisker angle');
pbins=-psth_time_before_trig:psth_binsize:psth_time_after_trig;
pbins=pbins(1:end-1)+psth_binsize/2;
%touch psths
parameter='Touch';Name='Touch';
exclude={'Puff';'Light';'Exclude';'Grooming'};
contPlot={MeanWt,MeanT,MeanAt};
[T_PSTH_instrate]=plot_Pop_PSTH_wrapper3(parameter, DiscreteData,touch_triglength_limit,ppms,psth_time_before_trig,psth_time_after_trig,exclude,psth_safety_margin,psth_binsize,[],parameter,{vpmr;pomr},{'VPM';'POm'},include,'latency',-inf,contPlot,{'meanbefore','scaleminmax','mono'},'baseline corr. whisker angle');

%%
smoothwidth=milliseconds(5);%ms
psth_peaks=nan(size(T_PSTH_instrate,1),3,2);
MP=smoothdata(P_PSTH_instrate,2,meth,smoothwidth,'SamplePoints',milliseconds(pbins));
MPs=MP(:,pbins>0 & pbins<75)-mean(MP(:,pbins<=0),2);
MPnoise=std(MPs(:,pbins<=0),0,2);
MT=smoothdata(T_PSTH_instrate,2,meth,smoothwidth,'SamplePoints',milliseconds(pbins));
MTs=MT(:,pbins>0 & pbins<75)-mean(MT(:,pbins<=0),2);
MTnoise=std(MTs(:,pbins<=0),0,2);
for n=1:size(MPs,1)
    if any([vpmra;pomra]==n)
        [psth_peaks(n,1,1),psth_peaks(n,2,1),psth_peaks(n,3,1),~] = findpeaks(MPs(n,:),pbins(pbins>0 & pbins<75),'MinPeakHeight',.99*max(MPs(n,:),[],2),'NPeaks',1,'WidthReference','halfheight');
    [h,l,w,~] = findpeaks(MTs(n,:),pbins(pbins>0 & pbins<75),'MinPeakHeight',.99*max(MTs(n,:),[],2),'NPeaks',1,'WidthReference','halfheight');
   if ~isempty(h)
    psth_peaks(n,1,2)=h;
   psth_peaks(n,2,2)=l;
   psth_peaks(n,3,2)=w;
   end
    % height loc width
    end
end
%%
mean_psth_param=cat(1,mean(psth_peaks(vpmr,:,:),1,'omitnan'),mean(psth_peaks(pomr,:,:),1,'omitnan'));
std_psth_param=cat(1,std(psth_peaks(vpmr,:,:),0,1,'omitnan'),std(psth_peaks(pomr,:,:),0,1,'omitnan'))./sqrt([numel(vpmr);numel(pomr)]);

fprintf('PUFF: vpm: peak rate=%.1f +- %.1f Hz, peak latency= %.1f +- %.1f ms, peak width= %.1f +- %.1f ms\n PUFF: POm: peak rate=%.1f +- %.1f Hz, peak latency= %.1f +- %.1f ms, peak width= %.1f +- %.1f ms\n',...
    mean_psth_param(1,1,1),std_psth_param(1,1,1),mean_psth_param(1,2,1),std_psth_param(1,2,1),mean_psth_param(1,3,1),std_psth_param(1,3,1),...
    mean_psth_param(2,1,1),std_psth_param(2,1,1),mean_psth_param(2,2,1),std_psth_param(2,2,1),mean_psth_param(2,3,1),std_psth_param(2,3,1));
for c=1:3
    for d=1:2
   [pnp(c,d)]=ranksum(psth_peaks(vpmr,c,d),psth_peaks(pomr,c,d));
    end
end
nucnames=categorical({'VPM';'POm'},{'VPM';'POm'},{'VPM';'POm'});
figure('Position',[ 910   176   1200   588])
tt=tiledlayout(1,3);
title(tt,'Puff peak params')
nexttile;
bar(nucnames,mean_psth_param(:,1,1));hold on;errorbar(nucnames,mean_psth_param(:,1,1),std_psth_param(:,1,1),'.k','CapSize',0); box off;
ylabel('peak height (Hz)')
title(sprintf('p=%.4g',pnp(1,1)))
nexttile;
bar(nucnames,mean_psth_param(:,2,1));hold on;errorbar(nucnames,mean_psth_param(:,2,1),std_psth_param(:,2,1),'.k','CapSize',0); box off;
ylabel('peak latency (ms)')
title(sprintf('p=%.4g',pnp(2,1)))
nexttile;
bar(nucnames,mean_psth_param(:,3,1));hold on;errorbar(nucnames,mean_psth_param(:,3,1),std_psth_param(:,3,1),'.k','CapSize',0); box off;
ylabel('peak width (ms)')
title(sprintf('p=%.4g',pnp(3,1)))
if write_figs;exportgraphics (gcf,[figdir Fig_no 'peak_params_puff.pdf'],'BackgroundColor','none','ContentType','vector');end


fprintf('TOUCH: vpm: peak rate=%.1f +- %.1f Hz, peak latency= %.1f +- %.1f ms, peak width= %.1f +- %.1f ms\n TOUCH: POm: peak rate=%.1f +- %.1f Hz, peak latency= %.1f +- %.1f ms, peak width= %.1f +- %.1f ms\n',...
    mean_psth_param(1,1,2),std_psth_param(1,1,2),mean_psth_param(1,2,2),std_psth_param(1,2,2),mean_psth_param(1,3,2),std_psth_param(1,3,2),...
    mean_psth_param(2,1,2),std_psth_param(2,1,2),mean_psth_param(2,2,2),std_psth_param(2,2,2),mean_psth_param(2,3,2),std_psth_param(2,3,2));

figure('Position',[ 910   176   1200   588])
tt=tiledlayout(1,3);
title(tt,'Touch peak params')
nexttile;
bar(nucnames,mean_psth_param(:,1,2));hold on;errorbar(nucnames,mean_psth_param(:,1,2),std_psth_param(:,1,2),'.k','CapSize',0); box off;
ylabel('peak height (Hz)')
title(sprintf('p=%.4g',pnp(1,2)))
nexttile;
bar(nucnames,mean_psth_param(:,2,2));hold on;errorbar(nucnames,mean_psth_param(:,2,2),std_psth_param(:,2,2),'.k','CapSize',0); box off;
ylabel('peak latency (ms)')
title(sprintf('p=%.4g',pnp(2,2)))
nexttile;
bar(nucnames,mean_psth_param(:,3,2));hold on;errorbar(nucnames,mean_psth_param(:,3,2),std_psth_param(:,3,2),'.k','CapSize',0); box off;
ylabel('peak width (ms)')
title(sprintf('p=%.4g',pnp(3,2)))
if write_figs;exportgraphics (gcf,[figdir Fig_no 'peak_params_touch.pdf'],'BackgroundColor','none','ContentType','vector');end

% %%
% psth_peaksR=permute(psth_peaks,[1 3 2]);
% psth_peaksR(:,:,1)=10*log10(psth_peaksR(:,:,1)./cat(2,MPnoise,MTnoise));
% write_figs=false;
% figure('Position',[ 910   176   1200   588])
% tiledlayout(1,3);
% nexttile;
% plot_ratecomp2(psth_peaksR(:,:,1),true(size(psth_peaksR,[1 2])),{'Puff','Touch'},{vpmr;pomr},{'VPM';'POm'},'mean','peak height');title('strength')
% nexttile;
% plot_ratecomp2(psth_peaksR(:,:,2),true(size(psth_peaksR,[1 2])),{'Puff','Touch'},{vpmr;pomr},{'VPM';'POm'},'mean','peak height');title('latency')
% nexttile;
% plot_ratecomp2(psth_peaksR(:,:,3),true(size(psth_peaksR,[1 2])),{'Puff','Touch'},{vpmr;pomr},{'VPM';'POm'},'mean','peak height');title('width')
% 
% 
% 
% if write_figs;exportgraphics (gcf,[figdir Fig_no 'ReliabilityTouchVsPuff.pdf'],'BackgroundColor','none','ContentType','vector');end
% %%

