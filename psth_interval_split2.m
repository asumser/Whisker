split_time_before_trig=50;%ms
split_time_after_trig=50;%ms
split_safety_margin=150;%time_after_trig;
isi=10;
Wdname='R:\Whisker\Wdec\Wdec_by_field.mat';
split_binsize=1;%ms
smoothwidth=milliseconds(7.5);%ms
meth='gaussian';
bins=-split_time_before_trig:split_binsize:split_time_after_trig;
pbins=bins(1:end-1);
pbins=pbins+mean(diff(pbins)/2);
binsR=[-30 5 35];
splits={[0 100/3 200/3 100],[0:20:100]};
splitname={'Tertile','Pentile'};


% splits={[0 25 50 75 100];[0 100/3 200/3 100];[0 50 100]};
% splitname={'Quartile';'Tertile';'median'};


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
split_var_nameC={'Acceleration';'Curvature';'Interval';'Velocity';'Setp';'Amp';'Angle';'Phase'};
mkabs=false;
normbefore=false;
fign_add='';
%%
split_var_nameC={'Acceleration';'Curvature';'Velocity';'Setp';'Angle';'Phase'};
mkabs=true;
normbefore=false;
fign_add='Abs_';
%%
split_var_nameC={'Acceleration';'Curvature';'Interval';'Velocity';'Setp';'Phase'};
mkabs=true;
normbefore=true;
fign_add='Abs_Norm';
%%
split_var_nameC={'Acceleration';'Curvature';'Interval';'Velocity';'Setp';'Amp';'Angle';'Phase'};
mkabs=false;
normbefore=true;
fign_add='Norm_';
%%
close all
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

        case {'Angle','Setp','Phase','Curvature','Acceleration','Velocity','Amp'}
            switch split_var_name
                case {'Angle','Setp','Amp'}
                    eval(sprintf('a=%s;',split_var_name));
                    switch split_var_name
                        case {'Angle','Setp'}
                            a=cellfun(@(x) x-median(x,'omitnan'),a,'UniformOutput',0);
                    end
                case {'Acceleration','Velocity'}
                    a=load(Wdname, split_var_name);
                    a=a.(split_var_name);
                case 'Phase'
                    a=Phase;
                    for x=1:numel(a)
                        % a{x}(smoothdata(Amp{x},2,'movmedian',.5,'SamplePoints',milliseconds(1:numel(Amp{x})))<3)=nan;%cellfun(@(x) x-median(x),a,'UniformOutput',0);
                        a{x}(Amp{x}<3)=nan;%cellfun(@(x) x-median(x),a,'UniformOutput',0);
                    end
                case 'Curvature'
                    a=load(Wdname, split_var_name);
                    a=a.(split_var_name);
                    a=cellfun(@(x) x-median(x,'omitnan'),a,'UniformOutput',0);
            end
            pbinsx=pbins(pbins<25 & pbins>-25);
            [~,sp_tr_P]=get_trig_cont_pop(P_trig,pbinsx,a,20,DiscreteData);
            [~,sp_tr_T]=get_trig_cont_pop(T_trig,pbinsx,a,20,DiscreteData);
            xT=cellfun(@isempty,sp_tr_T);
            xP=cellfun(@isempty,sp_tr_P);
            sp_tr_P(xP)={nan(size(pbinsx))};
            sp_tr_T(xT)={nan(size(pbinsx))};
            if mkabs
                sp_tr_P=cellfun(@abs,sp_tr_P,'UniformOutput',0);
                sp_tr_T=cellfun(@abs,sp_tr_T,'UniformOutput',0);
            end
            if normbefore
                sp_tr_P=cellfun(@(x) mean(x(:,pbinsx>0),2,'omitnan')-mean(x(:,pbinsx<0),2,'omitnan'),sp_tr_P,'UniformOutput',0);
                sp_tr_T=cellfun(@(x) mean(x(:,pbinsx>0),2,'omitnan')-mean(x(:,pbinsx<0),2,'omitnan'),sp_tr_T,'UniformOutput',0);
            else
                sp_tr_P=cellfun(@(x) mean(x(:,pbinsx>0),2,'omitnan'),sp_tr_P,'UniformOutput',0);
                sp_tr_T=cellfun(@(x) mean(x(:,pbinsx>0),2,'omitnan'),sp_tr_T,'UniformOutput',0);
            end
            sp_tr_P(xP)={[]};
            sp_tr_T(xT)={[]};
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
    %%
    for S=1:numel(splits)
        temp1=cat(1,sp_tr_Pall,sp_tr_Tall);
        temp1(isnan(temp1))=[];
        Ilim=prctile(temp1,splits{S});
        ix1=discretize(temp1,Ilim);
        mean_value_disc=nan(max(ix1),1);
        for j=1:max(ix1)
            mean_value_disc(j)=mean(temp1(ix1==j));
        end




        labels=cell(size(Ilim,2)-1,1);
        for l=1:numel(Ilim)-2
            labels{l}=sprintf('v<= %.3f',Ilim(l+1));
        end
        labels{end}=sprintf('v> %.3f',Ilim(end-1));




        ST=cell(numel(T_trig),numel(Ilim)-1,2);
        ST2=ST;
        STH=ST;
        STR=cell(numel(T_trig),numel(binsR)-1,2);
        STHc=ST;
        %puff trial split
        Climx=nan(numel(T_trig),numel(Ilim));
        for n=1:numel(T_trig)
            Clim=prctile(cat(1,sp_tr_P{n},sp_tr_T{n}),splits{S});
            Climx(n,:)=Clim;
            if ~isempty(RasterHP{n}) && ~isempty(sp_tr_P{n}) && ~any(n==[80 81])
                ix=discretize(sp_tr_P{n},Ilim)';
                sb=ix(RasterHP{n}(:,2));
                for x=1:numel(Ilim)-1

                    ST{n,x,1}=RasterHP{n}(sb==x,1);
                    ST2{n,x,1}=RasterHP{n}(sb==x,:);
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
            %touch trial split
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

        %%
        %normSN=1./sum(stn,[2 3]);
        %C=colororder;
        C=repmat(linspace(0,.8,numel(Ilim)-1)',1,3);

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
        exportgraphics(fig1,[figdir fign_add 'SplitQ_PT_' split_var_name '_' splitname{S} '_' datestr(now,'yymmdd') '.pdf'],'BackgroundColor','none','ContentType','vector')

        %%
        fig2= figure('Position',[56   196   641   883]);
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
        %         line(Ilim([2 2]),[0 1.25*max([H1 H2])],'Color','k');
        %         line(Ilim([end-1 end-1]),1.25*[0 max([H1 H2])],'Color','k');
        line(repmat(Ilim(:),1,2),[0 1.25*max([H1 H2])],'Color','k');
        box off
        legend(selS,'Location','eastoutside');legend boxoff
        title(split_var_name)

        STR=cellfun(@(x) x(2)-x(1),STRs);
        %STR=cellfun(@(x) log2(x(2)./x(1)),STRs);
        %  STR=cellfun(@(x) 2.*(x(2)-x(1))./x(1),STRs);
        for nx=1:2
            nt(nx)=nexttile;
            n=selc{nx};
            tempR=STR(n,:,:);
            tempR(isinf(tempR))=nan;
            tempp=nan(size(tempR,3),size(tempR,2));
            for pt=1:size(tempR,3)
                for tx=1:size(tempR,2)
                    if any(~isnan(tempR(:,tx,pt)))
                    tempp(pt,tx)=signrank(tempR(:,tx,pt),[],'tail','right');
                    end
                end
            end
            tempp=tempp(:);
            tempind=cell(size(tempp));
            for np=1:numel(tempp)
                if tempp(np)<.001
                    tempind{np}='***';
                elseif tempp(np)<.01
                    tempind{np}='**';
                elseif    tempp(np)<.05
                    tempind{np}='*';
                else
                    tempind{np}='n.s.';
                end
            end

            temppc=nan(1,size(tempR,2));

            for tx=1:size(tempR,2)
                if any(all(~isnan(tempR(:,tx,[1 2])),3))
                temppc(1,tx)=ranksum(tempR(:,tx,1),tempR(:,tx,2));
                end
            end

            tempindc=cell(size(temppc));
            for np=1:numel(temppc)
                if temppc(np)<.001
                    tempindc{np}='***';
                elseif temppc(np)<.01
                    tempindc{np}='**';
                elseif    temppc(np)<.05
                    tempindc{np}='*';
                else
                    tempindc{np}='n.s.';
                end
            end

            mean_tempR=squeeze(mean(tempR,1,'omitnan'));
            sem_tempR=squeeze(std(tempR,0,1,'omitnan')./sqrt(sum(~isnan(tempR))));
            b=bar(1:size(mean_tempR,1),mean_tempR);hold on;
            x = nan(size(mean_tempR,2), size(mean_tempR,1));
            for i = 1:size(mean_tempR,2)
                x(i,:) = b(i).XEndPoints;
            end

            errorbar(x',mean_tempR,sem_tempR,'k','linestyle','none');
            text(sort(x(:)),reshape(mean_tempR'+sem_tempR'+1,[],1),tempind,'HorizontalAlignment','center')

            line(x,repmat(1.15*max(mean_tempR+sem_tempR,[],'all'),size(x,1),size(x,2)),'color','k')
            text(mean(x,1),repmat(1.15*max(mean_tempR+sem_tempR,[],'all'),1,size(x,2)),tempindc,'HorizontalAlignment','center','VerticalAlignment','bottom')
            yl=ylim;
            yl(2)=1.35*max(mean_tempR+sem_tempR,[],'all');
            ylim(yl);
            hold off
            ylabel('rate change (Hz)')
            xlabel('quantile')
            legend(selS,'Location','eastoutside');legend boxoff
            box off
            cR=nan(2,2);
            for pt=1:2
                tempR2=tempR(:,:,pt);
                tempR2=tempR2(~(sum(isnan(tempR2),2)>=size(tempR2,2)-1),:);

                temp_mvd=repmat(mean_value_disc',size(tempR2,1),1);
                %           temp_mvd=repmat(1:3,size(tempR2,1),1);
                tempR2=tempR2(:);temp_mvd=temp_mvd(:);
                re=~isnan(tempR2);
                tempR2=tempR2(re);
                temp_mvd=temp_mvd(re);
                [c,p]=corrcoef(temp_mvd,tempR2);
                cR(pt,1)=c(2,1);
                cR(pt,2)=p(2,1);
            end
            title(selN{nx},sprintf('c(P)=%.2f, p(P)=%.6f, c(T)=%.2f, p(T)=%.6f',cR(1,1),cR(1,2),cR(2,1),cR(2,2)))
        end
        %%
        exportgraphics(fig2,[figdir fign_add 'SplitQ_PTm_' split_var_name '_' splitname{S} '_' datestr(now,'yymmdd') '.pdf'],'BackgroundColor','none','ContentType','vector')

    end
end
%%
%% single cells
%normSN=1./sum(stn,[2 3]);
%C=colororder;
C=repmat(linspace(0,.8,numel(Ilim)-1)',1,3);

fig1= figure('Position',[50         453        1059         425]);
tx=0;
ML=0;
fillpbins=[pbins fliplr(pbins)]';
nx=2;
tt=tiledlayout(numel(selc{nx}),1);
s=1;


for n=1:numel(selc{nx})
    tx=tx+1;
    tind(tx)=nexttile(tx);
    temp= cell2mat(STH(selc{nx}(n),:,s))./normSN(selc{nx}(n));
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

end

title(tt,sprintf('PSTHs, split by %s of %s',splitname{S},split_var_name))
%
for tx=1:numel(selc{nx})
    set(tind(tx),'YLim',[0 1.1*ML])
end
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
