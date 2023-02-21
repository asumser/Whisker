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
binsR=[-35 0 35];
splits={[0 100/3 200/3 100]};
splitname={'Tertile'};


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
Diffs=cell(numel(split_var_nameC),5);
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

        ST=cell(numel(T_trig),numel(Ilim)-1,2);
        pdrs=nan(numel(T_trig),2);
        pdks=pdrs;
        pdksn=pdrs;
        pdrsn=pdrs;
        meanratesplit=nan(numel(T_trig),numel(Ilim)-1,2,2);
        %puff trial split
        Climx=nan(numel(T_trig),numel(Ilim),2);
        for n=1:numel(T_trig)
            Clim=prctile(sp_tr_P{n},splits{S});
            Climx(n,:,1)=Clim;
            Clim(1)=-Inf;
            Clim(end)=Inf;
             Clim=prctile(sp_tr_P{n},splits{S});
            if ~isempty(RasterHP{n}) && ~isempty(sp_tr_P{n}) && ~any(n==[80 81])
                ix=quantileranks(sp_tr_P{n},numel(Clim)-1)';
                for x=1:numel(Ilim)-1
                    ixx=find(ix==x);
                    tempR=zeros(numel(ixx),2);
                    for xx=1:numel(ixx)
                        tempR(xx,:)= histcounts(RasterHP{n}(RasterHP{n}(:,2)==ixx(xx),1),binsR)./(diff(binsR)/1000);
                    end
                    ST{n,x,1}=tempR;
                    meanratesplit(n,x,:,1)=mean(tempR);
                end
                pdrsn(n,1)=ranksum(ST{n,1,1}(:,2)-mean(ST{n,1,1}(:,1)),ST{n,end,1}(:,2)-mean(ST{n,end,1}(:,1)));
                [~,pdksn(n,1)]=kstest2(ST{n,1,1}(:,2)-mean(ST{n,1,1}(:,1)),ST{n,end,1}(:,2)-mean(ST{n,end,1}(:,1)));
                pdrs(n,1)=ranksum(ST{n,1,1}(:,2),ST{n,end,1}(:,2));
                [~,pdks(n,1)]=kstest2(ST{n,1,1}(:,2),ST{n,end,1}(:,2));

            end
            %touch trial split
            if ~isempty(RasterHT{n}) && numel(sp_tr_T{n})>numel(splits{S}) && ~any(n==[80 81])
                Clim=prctile(sp_tr_T{n},splits{S});
                Climx(n,:,2)=Clim;
                ix=discretize(sp_tr_T{n},Clim);
                for x=1:numel(Ilim)-1
                    ixx=find(ix==x);
                    tempR=zeros(numel(ixx),2);
                    for xx=1:numel(ixx)
                        tempR(xx,:)= histcounts(RasterHT{n}(RasterHT{n}(:,2)==ixx(xx),1),binsR)./(diff(binsR)/1000);
                    end
                    ST{n,x,2}=tempR;
                    meanratesplit(n,x,:,2)=mean(tempR);
                end
                pdrsn(n,2)=ranksum(ST{n,1,2}(:,2)-mean(ST{n,1,2}(:,1)),ST{n,end,2}(:,2)-mean(ST{n,end,2}(:,1)));
                [~,pdksn(n,2)]=kstest2(ST{n,1,2}(:,2)-mean(ST{n,1,2}(:,1)),ST{n,end,2}(:,2)-mean(ST{n,end,2}(:,1)));
                pdrs(n,2)=ranksum(ST{n,1,2}(:,2),ST{n,end,2}(:,2));
                [~,pdks(n,2)]=kstest2(ST{n,1,2}(:,2),ST{n,end,2}(:,2));




            end
        end
        %%
        p=pdrs;
        plim=.5;
        ylab='rate';
        %         cr=meanratesplit-meanratesplit(:,:,1,:);
        cr=meanratesplit;
        comp13tert=figure('Position',[355   144   549   679]);
        plot_comp_tert_split(cr,vpmr,pomr,p,plim,ylab,split_var_name)
        exportgraphics(comp13tert,[figdir fign_add 'Indiv_' split_var_name '_' splitname{S} '_Comparison_' datestr(now,'yymmdd') '.pdf'],'BackgroundColor','none','ContentType','vector')

        p=pdksn;
        plim=.5;
        ylab='rate diff';
        cr=meanratesplit-meanratesplit(:,:,1,:);

        comp13tertnorm=figure('Position',[929   132   631   691]);
        plot_comp_tert_split(cr,vpmr,pomr,p,plim,ylab,split_var_name)
        exportgraphics(comp13tertnorm,[figdir fign_add 'Indiv_Diff_' split_var_name '_' splitname{S} '_Comparison_' datestr(now,'yymmdd') '.pdf'],'BackgroundColor','none','ContentType','vector')

        %%
Diffs{SN,1}=meanratesplit;
Diffs{SN,2}=pdrs;
Diffs{SN,3}=pdks;
Diffs{SN,4}=pdrsn;
Diffs{SN,5}=pdksn;
    end
end
%%
X=3;
pX=cat(3,Diffs{:,X});
pXs=pX<.05;
pXsn=squeeze(cat(2,mean(pXs(vpmr,:,:),1),mean(pXs(pomr,:,:),1)))';
bC=categorical(split_var_nameC,split_var_nameC,split_var_nameC);
figure, bar(bC,pXsn*100)
ylabel('% significantly different 1st & 3rd tertile responses')
legend('vpm puff','vpm touch','pom puff','pom touch')