split_time_before_trig=50;%ms
split_time_after_trig=50;%ms
split_safety_margin=150;%time_after_trig;
isi=10;
Wdname='D:\Whisker\Wdec\Wdec_by_field.mat';
split_binsize=1;%ms
smoothwidth=milliseconds(10);%ms
meth='gaussian';
bins=-split_time_before_trig:split_binsize:split_time_after_trig;
pbins=bins(1:end-1);
pbins=pbins+mean(diff(pbins)/2);
binsR=[-50 0 50];
split_var_nameC={'Interval';'Curvature';'Acceleration';'Velocity'};
% splits={[0 25 50 75 100];[0 100/3 200/3 100];[0 50 100]};
% splitname={'Quartile';'Tertile';'median'};
splits={[0 100/3 200/3 100]};
splitname={'Tertile'};
%Puff
parameter='Puff';
exclude={'Pole';'Light';'Exclude';'Grooming'};
triglength_limit=[28 32];%ms
tempDD=struct2cell(DiscreteData);
spikes=squeeze(tempDD(1,1,:));clear temp*
[RasterHP,P_trig]=get_raster_pop_no_self_safety(parameter,spikes, DiscreteData,triglength_limit,ppms,split_time_before_trig,split_time_after_trig,exclude,split_safety_margin,split_binsize,isi,[]);
[RasterPP,~]=get_raster_pop_no_self_safety(parameter,P_trig, DiscreteData,triglength_limit,ppms,split_time_before_trig,split_time_after_trig,exclude,split_safety_margin,split_binsize,isi,[]);

parameter='Touch';
exclude={'Puff';'Light';'Exclude';'Grooming'};
triglength_limit=[-inf inf];%ms
[RasterHT,T_trig]=get_raster_pop_no_self_safety(parameter,spikes, DiscreteData,triglength_limit,ppms,split_time_before_trig,split_time_after_trig,exclude,split_safety_margin,split_binsize,isi,[]);


selc={vpmr;pomr};
selN={'VPM';'POm'};
selS={'Puff';'Touch'};
%%

for SN=1:numel(split_var_nameC)
    split_var_name=split_var_nameC{SN};
    % smoothker1=[1:1000 999:-1:0];
    % smoothkerLU=linspace(-smoothwidth/2,smoothwidth/2,numel(smoothker1));
    % smoothker2=interp1(smoothkerLU,smoothker1,bins,'linear',0);
    % smoothker2=smoothker2./sum(smoothker2);

    %


    switch split_var_name
        case 'Interval'
            sp_tr_T=cell(size(T_trig));
            for n=1:numel(RasterHT)
                if ~isempty(T_trig{n})

                    sp_tr_T{n}=cat(1,T_trig{n}(1),diff(T_trig{n}))/20000;
                end
            end
            sp_tr_P=cellfun(@(x) cat(1,x(1),diff(x))/20000,P_trig,'UniformOutput',0);

        case {'Angle','Setpoint'}
            a=load(Wdname, split_var_name);
            a=a.(split_var_name);
            a=cellfun(@(x) x-median(x),a,'UniformOutput',0);
            [~,sp_tr_P]=get_trig_cont_pop(P_trig,pbins(pbins<25 & pbins>-25),a,20,DiscreteData);
            sp_tr_P=cellfun(@(x) mean(x,2,'omitnan')',sp_tr_P,'UniformOutput',0);
            [~,sp_tr_T]=get_trig_cont_pop(T_trig,pbins(pbins<0 & pbins>-50),a,20,DiscreteData);
            sp_tr_T=cellfun(@(x) mean(x,2,'omitnan')',sp_tr_T,'UniformOutput',0);
        case 'Phase'
            a=load(Wdname, split_var_name);
            load(Wdname, 'Amp_extrema')
            a=a.(split_var_name);
            for x=1:numel(a)
                a{x}(Amp{x}<1.5)=nan;%cellfun(@(x) x-median(x),a,'UniformOutput',0);

            end
            [~,sp_tr_P]=get_trig_cont_pop(P_trig,pbins(pbins<10 & pbins>-10),a,20,DiscreteData);
            sp_tr_P=cellfun(@(x) mean(x,2,'omitnan')',sp_tr_P,'UniformOutput',0);
            [~,sp_tr_T]=get_trig_cont_pop(T_trig,pbins(pbins<0 & pbins>-50),a,20,DiscreteData);
            sp_tr_T=cellfun(@(x) mean(x,2,'omitnan')',sp_tr_T,'UniformOutput',0);
        case 'Curvature'
            a=load(Wdname, split_var_name);
            a=a.(split_var_name);
            a=cellfun(@(x) x-median(x),a,'UniformOutput',0);
            [~,sp_tr_P]=get_trig_cont_pop(P_trig,pbins(pbins<25 & pbins>0),a,20,DiscreteData);
            sp_tr_P=cellfun(@(x) mean(x,2,'omitnan')',sp_tr_P,'UniformOutput',0);
            [~,sp_tr_T]=get_trig_cont_pop(T_trig,pbins(pbins<25 & pbins>0),a,20,DiscreteData);
            sp_tr_T=cellfun(@(x) mean(x,2,'omitnan')',sp_tr_T,'UniformOutput',0);
        case {'Acceleration','Velocity','Amp_extrema'}
            a=load(Wdname, split_var_name);
            a=a.(split_var_name);
            %a=cellfun(@(x) x-median(x),a,'UniformOutput',0);
            [~,sp_tr_P]=get_trig_cont_pop(P_trig,pbins(pbins<25 & pbins>-25),a,20,DiscreteData);
            sp_tr_P=cellfun(@(x) mean(x,2,'omitnan')',sp_tr_P,'UniformOutput',0);
            [~,sp_tr_T]=get_trig_cont_pop(T_trig,pbins(pbins<25 & pbins>-25),a,20,DiscreteData);
            sp_tr_T=cellfun(@(x) mean(x,2,'omitnan')',sp_tr_T,'UniformOutput',0);

    end

    sp_tr_T=cellfun(@(x) reshape(x,[],1),sp_tr_T,'UniformOutput',0);
    sp_tr_P=cellfun(@(x) reshape(x,[],1),sp_tr_P,'UniformOutput',0);
    sp_tr_Tall=cat(1,sp_tr_T{cat(1,selc{:})});
    sp_tr_Pall=cat(1,sp_tr_P{cat(1,selc{:})});

    %Ilim=prctile(cat(2,sp_tr_Pall,sp_tr_Tall),[0 100/3 200/3 100]);
    %Ilim=prctile(cat(2,sp_tr_Pall,sp_tr_Tall),[0 25 50 75 100]);
    %Ilim=prctile(cat(2,sp_tr_Pall,sp_tr_Tall),[0 100/2 100]);

    %
    %Ilim=[0 0.25 1.1 inf];
    %

    %Ilim=[0 inf];
    %
    for S=1:numel(splits)
        temp1=cat(1,sp_tr_Pall,sp_tr_Tall);
        temp1(isnan(temp1))=[];
        Ilim=prctile(temp1,splits{S});





        labels=cell(size(Ilim,2)-1,1);
        for l=1:numel(Ilim)-2
            labels{l}=sprintf('v<= %.3f',Ilim(l+1));
        end
        labels{end}=sprintf('v> %.3f',Ilim(end-1));




        ST=cell(numel(T_trig),numel(Ilim)-1,2);
        STH=ST;
        STR=cell(numel(T_trig),numel(binsR)-1,2);
        STHc=ST;
        for n=1:numel(T_trig)
            if ~isempty(RasterHP{n}) && ~isempty(sp_tr_P{n}) && ~any(n==[80 81])
                ix=discretize(sp_tr_P{n},Ilim)';
                sb=ix(RasterHP{n}(:,2));
                for x=1:numel(Ilim)-1

                    ST{n,x,1}=RasterHP{n}(sb==x,1);
                    STH{n,x,1}=histcounts(ST{n,x,1},bins)./sum(ix==x);
                    STR{n,x,1}=histcounts(ST{n,x,1},binsR)./sum(ix==x)./(diff(binsR)/1000);
                    if sum(ix==x) <=1
                        STH{n,x,2}(:)=nan;
                        STR{n,x,2}(:)=nan;
                    end
                    %      STHc{n,x,1}=conv(STH{n,x,1},smoothker2,'same');
                end
            else
                STH(n,:,1)={nan(1,numel(pbins))};
                STR(n,:,2)={nan(1,numel(binsR)-1)};
                %   STHc(n,:,2)={nan(1,numel(pbins))};
            end
            if ~isempty(RasterHT{n}) && ~isempty(sp_tr_T{n}) && ~any(n==[80 81])
                ix=discretize(sp_tr_T{n},Ilim);
                sb=ix(RasterHT{n}(:,2));
                for x=1:numel(Ilim)-1

                    ST{n,x,2}=RasterHT{n}(sb==x,1);
                    STH{n,x,2}=histcounts(ST{n,x,2},bins)./sum(ix==x);
                    STR{n,x,2}=histcounts(ST{n,x,2},binsR)./sum(ix==x)./(diff(binsR)/1000);
                    if sum(ix==x) <=1
                        STH{n,x,2}(:)=nan;
                        STR{n,x,2}(:)=nan;
                    end
                    %  STHc{n,x,2}=conv(STH{n,x,2},smoothker2,'same');
                end
            else
                STH(n,:,2)={nan(1,numel(pbins))};
                STR(n,:,2)={nan(1,numel(binsR)-1)};
                %   STHc(n,:,2)={nan(1,numel(pbins))};
            end
        end
        %
        stn=cellfun(@numel,ST);
        normSN=ones(numel(T_trig),1);
        STR(80,:,2)={nan(1,numel(binsR)-1)};
        STR(81,:,2)={nan(1,numel(binsR)-1)};
        STR(80,:,1)={nan(1,numel(binsR)-1)};
        STR(81,:,1)={nan(1,numel(binsR)-1)};
        STRs=STR;


        %normSN=1./sum(stn,[2 3]);
        C=colororder;
        % C=repmat((0:.2:.8)',1,3);

        fig1= figure('Position',[50         453        1059         425]);
        tx=0;
        ML=0;
        fillpbins=[pbins fliplr(pbins)]';
        tt=tiledlayout(2,2);

        for s=1:2

            for nx=1:2
                n=selc{nx};
                tx=tx+1;
                tind(tx)=nexttile(tx);
                temp= cell2mat(STH(n,:,s))./normSN(n);
                temp=reshape(temp,size(temp,1),numel(bins)-1,size(STH,2));
                %         SEM=permute(std(temp,0,1,'omitnan')./sqrt(sum(~any(isnan(temp),[2 3]))),[2 3 1]);
                %         SEM=smoothdata(SEM,1,meth,smoothwidth,'SamplePoints',milliseconds(pbins));
                M=permute(mean(temp,1,'omitnan'),[2 3 1]);
                M=smoothdata(M,1,meth,smoothwidth,'SamplePoints',milliseconds(pbins));
                for b=1:size(M,2)
                    %             SEMp=cat(1,M(:,b)+SEM(:,b),flipud(M(:,b)-SEM(:,b)));
                    %             fill(fillpbins, SEMp, C(b,:),'FaceAlpha',.5,'EdgeColor','none');
                    hold on;
                    %plot(pbins',conv(M(:,b),smoothker2,'same'),'Color',C(b,:));box off;
                    plot(pbins',M(:,b),'Color',C(b,:));box off;

                end

                ML=max(ML,max(M,[],'all'));
                legend(labels);legend('boxoff');
                xlabel('time around deflection (ms)')
                ylabel({'Population mean';'spike prob. per bin'})
                title(sprintf('%s (N=%u) %s',selN{nx},sum(any(~all(isnan(temp),2),3)),selS{s}));
            end
        end
        title(tt,sprintf('PSTHs, split by %s of %s',splitname{S},split_var_name))
        %
        for tx=1:4
            set(tind(tx),'YLim',[0 1.1*ML])
        end
        exportgraphics(fig1,[figdir 'SplitQ_PT_' split_var_name '_' splitname{S} '_' datestr(now,'yymmdd') '.pdf'],'BackgroundColor','none','ContentType','vector')

        %%
        fig2= figure('Position',[ 50   453   397   425]);
        tt=tiledlayout(3,1);
        nexttile;
        switch split_var_name
            case  'Interval'
            vbins=[linspace(Ilim(1),1.25,50) Ilim(end)];
            pvbins=vbins(1:end-1)+median(diff(vbins))/2;
            otherwise
            vbins=linspace(Ilim(1),Ilim(end),51);
            pvbins=vbins(1:end-1)+median(diff(vbins))/2;

        end

        H1=histcounts(sp_tr_Pall,vbins,'Normalization','probability');
        H2=histcounts(sp_tr_Tall,vbins,'Normalization','probability');
        bar(pvbins,H1,1,'FaceColor',C(1,:),'EdgeColor','none');hold on
        bar(pvbins,H2,1,'FaceColor',C(2,:),'EdgeColor','none','FaceAlpha',.7);hold on
        line(Ilim([2 2]),[0 1.25*max([H1 H2])],'Color','k');
        line(Ilim([end-1 end-1]),1.25*[0 max([H1 H2])],'Color','k');
        box off
        legend(selS,'Location','eastoutside');legend boxoff
        title(split_var_name)

        STR=cellfun(@(x) x(2)-x(1),STRs);
        for nx=1:2
            nexttile;
            n=selc{nx};
            tempR=STR(n,:,:);
            tempR(isinf(tempR))=nan;
            mean_tempR=squeeze(mean(tempR,1,'omitnan'));
            sem_tempR=squeeze(std(tempR,0,1,'omitnan')./sqrt(sum(~isnan(tempR))));
            b=bar(1:size(mean_tempR,1),mean_tempR);hold on;
            x = nan(size(mean_tempR,2), size(mean_tempR,1));
            for i = 1:size(mean_tempR,2)
                x(i,:) = b(i).XEndPoints;
            end

            errorbar(x',mean_tempR,sem_tempR,'k','linestyle','none');
            hold off
            ylabel('rate change (Hz)')
            xlabel('quantile')
            legend(selS,'Location','eastoutside');legend boxoff
            box off
            title(selN{nx})


        end
        %%
        exportgraphics(fig2,[figdir 'SplitQ_PTm_' split_var_name '_' splitname{S} '_' datestr(now,'yymmdd') '.pdf'],'BackgroundColor','none','ContentType','vector')

    end
end
%%
% figure('Position',[273   110   880   762])
% for n=1:numel(RasterHP)
%     if size(RasterHP{n},1)>10
%
%         %[IPIs,IPIsortIdx]=sort(IPI{n});
%
%         scatter(RasterPP{n}(:,1),RasterPP{n}(:,2),'ro');hold on
%         ax(1)=scatter(RasterHP{n}(:,1),RasterHP{n}(:,2),'k.');box off;hold off;
%         %tIPI=IPI{n}(IPIsortIdx(RasterHP{n}(:,2)))/20;
%
%         %plot(-IPIs,1:numel(P_trig{n}));hold off
%         yt=get(gca,'YTick');
%         ipit=interp1(1:numel(IPIs),IPIs,yt(1:end-1),'linear','extrap');
%         set(gca,'YTick',yt(1:end-1),'YTickLabel',ipit);
%         xlim([-time_before_trig time_after_trig]);
%         ylabel('Inter Puff interval (ms)')
%         drawnow;
%         ginput(1);
%     end
% end
% %%
%
%
% for n=1:numel(RasterHT)
%     if size(RasterHT{n},1)>10
%         figure('Position',[273   110   880   762])
%         sp_tr_T=cat(2,T_trig{n}(1),diff(T_trig{n}));
%         [ITIs,ITIsortIdx]=sort(sp_tr_T);
%         ITIs=ITIs/20;
%         ax(1)=scatter(RasterHT{n}(:,1),IPIsortIdx(RasterHT{n}(:,2)),'k.');box off;
%         yt=get(gca,'YTick');
%         ipit=interp1(1:numel(ITIs),ITIs,yt(1:end-1),'linear','extrap');
%         set(gca,'YTick',yt(1:end-1),'YTickLabel',ipit);
%         xlim([-time_before_trig time_after_trig]);
%         ylabel('Inter Touch interval (ms)')
%     end
% end
