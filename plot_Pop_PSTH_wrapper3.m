function [isr_res]=plot_Pop_PSTH_wrapper3(parameter, DiscreteData,triglength_limit,ppms,time_before_trig,time_after_trig,exclude,safety_margin,binsize,trigstartend,Name,sel,selnames,force_inclusive,sortby,scalem,varargin)



if ischar(parameter)
    [trig,trigL]=select_trigs(parameter,DiscreteData,triglength_limit,ppms,trigstartend);
    trig1=trig;
    for R=1:size(DiscreteData,2);[trig{R}]=exclude_reg_puff(trig{R});end
    %  trig=cleanup_self_safety(time_before_trig*ppms,time_after_trig*ppms,trig); %only isolated
    trig=cleanup_sections_exclusive(DiscreteData,exclude,safety_margin*ppms,trig);
    if ~isempty(force_inclusive)
        [trig]=cleanup_sections_inclusive(DiscreteData,force_inclusive{1},force_inclusive{2}*ppms,trig);
    end
    for R=1:size(DiscreteData,2)
        td=trig{R}(:)-trig1{R}(:)';
        [a,b]=find(td==0);
        trigL{R}=trigL{R}(b);
        
        %trig{R}=trig{R}-time_before_trig*ppms;
        
    end
else
    trig=parameter(:,1);
    %     for R=1:size(DiscreteData,2)
    %          trig{R}=trig{R}-time_before_trig*ppms;
    %     end
    if size(parameter,2)>1
        trigL=parameter(:,2);
    else
        trigL=parameter(:,1);
    end
end
tempDD=struct2cell(DiscreteData);
spikes=squeeze(tempDD(1,1,:));
clear temp*
bins=(-time_before_trig:binsize:time_after_trig);
maxT=time_after_trig/1000;
[isr_res,p,lat]=DD_pop_PSTH(trig,spikes,time_before_trig,time_after_trig,binsize,ppms,find(bins<=0,1,'last'),'none',maxT,trigL);
isr_res=cat(2,isr_res{:})';
% n_zeta=4;
% zeta_p=ones(numel(spikes),1);
% zeta_latency=ones(numel(spikes),4);
% for R=1:size(DiscreteData,2)
%     if ~isempty(trig{R})
%         [z1,z2,z3,z4]  = getZeta(spikes{R}(:)/(ppms*1000),[trig{R}(:)/(ppms*1000) trig{R}(:)/(ppms*1000)+trigL{R}(:)/(ppms*1000)],maxT+time_before_trig/1000,100,0,n_zeta,[-inf maxT+time_before_trig/1000],false);
%         if ~isempty(z4)
%             zeta_p(R)=z1;
%             zeta_latency(R,:)=z2;
%             sZETA(R)=z3;
%             sRate(R)=z4;
%             sR{R,1}=sRate(R).vecT;
%             sR{R,2}=sRate(R).vecRate;
%         end
%     end
% end
% t_res=linspace(0,maxT+time_before_trig/1000,300);
% isr_res=nan(numel(zeta_p),numel(t_res));
% for R=1:size(sRate,2)
%     [~,ia,~]=unique(sRate(R).vecT);
%     if ~isempty(ia)
%         isr_res(R,:)=interp1(sRate(R).vecT(ia),sRate(R).vecRate(ia),t_res,'linear');
%     end
% end

smoothwidth=milliseconds(5);%ms
meth='gaussian';
pbins=bins(1:end-1)+diff(bins);
isr_resF=smoothdata(isr_res,2,meth,smoothwidth,'SamplePoints',milliseconds(pbins));

%peakTime=cat(1,sRate(:).dblPeakTime);
tiledlayout(5,numel(sel))
mscale=prctile(prctile(isr_res(cat(1,sel{:}),:),99,2),100);
m=0;
%t_res_corr=t_res*1000-time_before_trig;
for s=1:numel(sel)
    H=isr_resF(sel{s},:)';
    cellnotnan=sum(isnan(H))==0;
    
    switch sortby
        case 'none'
            sortval=cellnotnan;
        case 'latency'
            [~,peakTimeSel]=max(movmean(H,5,1),[],1);
            [~,sortval]=sort(peakTimeSel);
            cellnotnan=cellnotnan(sortval);
            sortval=sortval(cellnotnan);
    end
    %H2=log2(H./mean(H(pbins<0,:),'omitnan'));
    H2=max(log2(max(H-mean(H(pbins<0,:),'omitnan'),0)),-inf);
    nexttile(s,[2 1]);
    % imagesc(pbins,1:sum(cellnotnan),(H(:,cellnotnan)-mean(H(baseline,cellnotnan)))');
    %imagesc(pbins,1:sum(cellnotnan),(H(:,cellnotnan)./std(H(:,cellnotnan),0,1))');
    %imagesc(t_res*1000-time_before_trig,1:sum(cellnotnan),zscore(H(:,sortval),0,1)',[0 5]);
    imagesc(pbins,1:sum(cellnotnan),H2(:,sortval)',[0 10]);
    %set(gca,'ColorScale','log')
    colormap(gray);
    cb(s)=colorbar;cb(s).Label.String='dRate';
    cb(s).TickLabels=cellstr(num2str(cb(s).Ticks(:).^2));
    %  set(gca,'CLim',[0 0.8]*max(abs(get(gca,'CLim'))))
    % set(gca,'CLim',[0 5])
    title(selnames{s});box off;set(gca,'YTick',[]);
    %H=H./sum(H,1,'omitnan');
    %H=H./std(H,0,1,'omitnan');
    m=max(m,max(mean(H,2,'omitnan')+std(H,0,2,'omitnan')/sqrt(sum(~isnan(H(1,:))))));
    h(s)=nexttile(s+3*numel(sel),[2 1]);
    
    bar(pbins,mean(H,2,'omitnan'),1,'facecolor',[.5 .5 .5],'edgecolor','none'); hold on
    plot(pbins,mean(H,2,'omitnan')+std(H,0,2,'omitnan')/sqrt(sum(~isnan(H(1,:)))),'k--');
    plot(pbins,mean(H,2,'omitnan')-std(H,0,2,'omitnan')/sqrt(sum(~isnan(H(1,:)))),'k--');
    set(gca,'xminortick','on','box','off')
    xlim(pbins([1 end]))
    ylabel('rate (Hz)')
    xlabel(['time (s)'])
    if size(varargin,1)~=0
        w(s)=nexttile(s+2*numel(sel));
        a=1;
        ccont=cat(1,varargin{1}{a}{sel{s}});
        if strcmp(varargin{2}{a},'meanbefore')
            ccont=ccont-mean(ccont(:,pbins<0),2);
        end
        ccont_mean=mean(ccont,1,'omitnan');
        
        plot(pbins,ccont_mean,'r');box off
        xlim(pbins([1 end]))
        scaley{s}=ylim;
        ylabel(varargin{3})
        if numel(varargin{1})>1
            yyaxis right
            for a=2:numel(varargin{1})
                ccont=cat(1,varargin{1}{a}{sel{s}});
                if strcmp(varargin{2}{a},'meanbefore')
                    ccont=ccont-mean(ccont(:,pbins<0),2);
                end
                ccont_mean=mean(ccont,1,'omitnan');
                if strcmp(varargin{2}{a},'scaleminmax')
                    h(s)=nexttile(s+3*numel(sel),[2 1]);
                    ccont_mean=ccont_mean-min(ccont_mean);
                    ccont_mean=ccont_mean./max(ccont_mean);
                    ccont_mean=normalize(ccont_mean,'range',[0 m]);
                elseif contains(varargin{2}{a},'scale')
                    w(s)=nexttile(s+2*numel(sel));
                    scale=str2double(varargin{2}{a}(6:end));
                    ccont_mean=ccont_mean*scale;
                elseif contains(varargin{2}{a},'mono')
                    w(s)=nexttile(s+2*numel(sel));
                    ccont_mean(:)=mean(ccont_mean);
                end
                plot(pbins,ccont_mean,'b-');box off
                xlim(pbins([1 end]))
            end
        end
        %         burst_var=cat(2,varargin{1}{sel{s}});
        %         if numel(varargin)>1 && varargin{2}
        %             burst_var=burst_var./sum(burst_var,1,'omitnan');
        %         end
        %         yyaxis right
        %         m=max(m,max(mean(burst_var,2,'omitnan')));
        %
        %         b(s)=stairs(pbins-diff(pbins(1:2))/2,mean(burst_var,2,'omitnan'),'r');
        %         ylabel(varargin{3})
        
    end
    title([selnames{s} ' Pop ' Name ' Responses (N=' num2str(sum(~isnan(H(1,:)))) ')'])
    
end
scaley2=cat(1,scaley{:});
scaley3=[min(scaley2(:,1)) max(scaley2(:,2))];
if m<scalem
    m=scalem;
end
%     title([selnames{s} ' Pop ' Name ' Responses (N=' num2str(sum(~isnan(H(1,:)))) ')'])
for s=1:numel(sel)
    if size(varargin,1)~=0
        %
        set(h(s),'YLim',[0 m]);
        yyaxis (w(s),'left')
        set(w(s),'YLim',scaley3);
        yyaxis (w(s),'right')
        temp=get(w(s),'YLim');
        temp(1)=0;
        set(w(s),'YLim',temp);
    else
        set(h(s),'YLim',[0 m]);
    end
end
