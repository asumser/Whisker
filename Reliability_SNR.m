%%
parameter='Puff';
exclude={'Pole';'Light';'Exclude';'Grooming'};
triglength_limit=[28 32];%ms
raster_binsize=1;
tempDD=struct2cell(DiscreteData);
spikes=squeeze(tempDD(1,1,:));clear temp*
[RasterHP,P_trig]=get_raster_pop_no_self_safety(parameter,spikes, DiscreteData,triglength_limit,ppms,rate_time_before_trig,rate_time_after_trig,exclude,psth_safety_margin,raster_binsize,isi,[]);
[RasterPP,~]=get_raster_pop_no_self_safety(parameter,P_trig, DiscreteData,triglength_limit,ppms,rate_time_before_trig,rate_time_after_trig,exclude,psth_safety_margin,raster_binsize,isi,[]);

parameter='Touch';
exclude={'Puff';'Light';'Exclude';'Grooming'};
triglength_limit=[-inf inf];%ms
[RasterHT,T_trig]=get_raster_pop_no_self_safety(parameter,spikes, DiscreteData,triglength_limit,ppms,rate_time_before_trig,rate_time_after_trig,exclude,psth_safety_margin,raster_binsize,isi,[]);
%%
win=20;
SNR=nan(numel(RasterHP),2);
Reliability=nan(numel(RasterHP),2);
for n=1:numel(RasterHP)
    tempR=RasterHP{n};
    relT=-psth_time_before_trig:psth_time_after_trig;
    temp_trials=zeros(size(P_trig{n},1),numel(relT));
    for s=1:size(tempR,1)
        temp_trials(tempR(s,2),relT>=tempR(s,1)-win/2 & relT<tempR(s,1)+win/2)=temp_trials(tempR(s,2),relT>=tempR(s,1)-win/2 & relT<tempR(s,1)+win/2)+1;
    end
    temp_trials=temp_trials(:,relT>0);
    SNR(n,1)=var(mean(temp_trials,1),0,2)./mean(var(temp_trials,0,2),1);
    % NT=size(temp_trials,1);
    % Sparse=(1-(sum(NeuronTestActMean/NT,2).^2./sum(NeuronTestActMean.^2/NT,2))).*(1/(1-1/NT));
    %Reliability=nan(size(NeuronTestActMean,3),1);
    inc=triu(true(size(temp_trials,1)),1);
    c=corr(temp_trials');
    Reliability(n,1)=mean(c(inc),'all','omitnan');
    if size(T_trig{n},1)>0
    tempR=RasterHT{n};
    relT=-psth_time_before_trig:psth_time_after_trig;
    temp_trials=zeros(size(T_trig{n},1),numel(relT));
    for s=1:size(tempR,1)
        temp_trials(tempR(s,2),relT>=tempR(s,1)-win/2 & relT<tempR(s,1)+win/2)=temp_trials(tempR(s,2),relT>=tempR(s,1)-win/2 & relT<tempR(s,1)+win/2)+1;
    end
    temp_trials=temp_trials(:,relT>0);
    SNR(n,2)=var(mean(temp_trials,1),0,2)./mean(var(temp_trials,0,2),1);
    % NT=size(temp_trials,1);
    % Sparse=(1-(sum(NeuronTestActMean/NT,2).^2./sum(NeuronTestActMean.^2/NT,2))).*(1/(1-1/NT));
    %Reliability=nan(size(NeuronTestActMean,3),1);
    inc=triu(true(size(temp_trials,1)),1);
    c=corr(temp_trials');
    Reliability(n,2)=mean(c(inc),'all','omitnan');
    end

end
figure, scatter(SNR(pomr,1),SNR(pomr,2),[],'bo');hold on;scatter(SNR(vpmr,1),SNR(vpmr,2),[],'ro');line([0 .75],[0 .75],'Color','k');axis equal
%%
figure('Position',[674   573   886   367])
tiledlayout(1,2);
nexttile;
plot(SNR(pomr,1),SNR(pomr,2),'bo');hold on;
plot(SNR(vpmr,1),SNR(vpmr,2),'ro');
line([0 .75],[0 .75],'Color','k');axis equal;box off
xlabel('Puff');
ylabel('Touch')
title('SNR')
legend('POm','VPM','Location','northwest');legend boxoff
xlim([0 1]);
ylim([0 1]);
nexttile;
plot(Reliability(pomr,1),Reliability(pomr,2),'bo');hold on;
plot(Reliability(vpmr,1),Reliability(vpmr,2),'ro');
line([0 .75],[0 .75],'Color','k');axis equal;box off
title('Reliability')
xlabel('Puff');
xlim([0 1]);
ylim([0 1]);
ylabel('Touch')
legend('POm','VPM','Location','northwest');legend boxoff
%%
figure('Position',[ 910   176   500   588])
plot_ratecomp2(Reliability,true(size(Reliability)),{'Puff','Touch'},{vpmr;pomr},{'VPM';'POm'},'mean','Reliability');title('Reliability')
if write_figs;exportgraphics (gcf,[figdir Fig_no 'ReliabilityTouchVsPuff.pdf'],'BackgroundColor','none','ContentType','vector');end

figure('Position',[ 910   176   500   588])
plot_ratecomp2(SNR,true(size(SNR)),{'Puff','Touch'},{vpmr;pomr},{'VPM';'POm'},'mean','SNR');title('SNR')
if write_figs;exportgraphics (gcf,[figdir Fig_no 'SNRTouchVsPuff.pdf'],'BackgroundColor','none','ContentType','vector');end