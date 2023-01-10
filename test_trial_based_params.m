splits2={linspace(0,100,4);linspace(0,100,11);linspace(0,100,26);linspace(0,100,101)};
ty='Spearman';
 mkabs=false;
    mkindiv=false;
for S=2


   
    indiv_split_response=cell(2,numel(splits2{S})-1,numel(sp_tr_T));
    ks_p_TP_split=nan(numel(sp_tr_T),numel(splits2{S})-1);
    ks_p_TP13_split=nan(numel(sp_tr_T),2);
    corr_TP_split=nan(numel(sp_tr_T),2,3);
    temp1=cat(1,sp_tr_Pall,sp_tr_Tall);
    if mkabs;temp1=abs(temp1);end
    temp1(isnan(temp1))=[];
    Ilim=prctile(temp1,splits2{S});
    for n=1:numel(sp_tr_T)
        tempparam=cat(1,sp_tr_P{n},sp_tr_T{n});
        if mkabs;tempparam=abs(tempparam);end
        if ~isempty(tempparam)
            if mkindiv
                Clim=prctile(tempparam,splits2{S});
            else
                Clim=Ilim;
                Clim=[linspace(prctile(temp1,5),prctile(temp1,95),12)];
            end

            ixp=discretize(tempparam,Clim)';
            ix=discretize(sp_tr_P{n},Clim)';
            mean_param_value_discr=nan(3,numel(Clim)-1);
            for c=1:numel(splits2{S})-1
                cix=find(ix==c);
                mean_param_value_discr(3,c)=mean(tempparam(ixp==c));
                mean_param_value_discr(1,c)=mean(sp_tr_P{n}(ix==c));
                temp_num=nan(numel(cix),2);
                for nix=1:numel(cix)
                    temp=RasterHP{n}(RasterHP{n}(:,2)==cix(nix),1);
                    temp_num(nix,:)=[sum(temp<=0) sum(temp>0)];

                end
                if ~isempty(temp_num)
                indiv_split_response{1,c,n}=temp_num(:,2)-mean(temp_num(:,1));
                else
indiv_split_response{1,c,n}=[];
                end
            end
            ix=discretize(sp_tr_T{n},Clim)';
            for c=1:numel(splits2{S})-1
                cix=find(ix==c);
                mean_param_value_discr(2,c)=mean(sp_tr_T{n}(ix==c));
                temp_num=nan(numel(cix),2);
                for nix=1:numel(cix)
                    temp=RasterHT{n}(RasterHT{n}(:,2)==cix(nix),1);
                    temp_num(nix,:)=[sum(temp<=0) sum(temp>0)];

                end
                indiv_split_response{2,c,n}=temp_num(:,2)-mean(temp_num(:,1));
            end
            for pt=1:2

                tempRx=indiv_split_response(pt,:,n);
                tempPx=[];
                for c=1:numel(splits2{S})-1
                    
                    tempPx=[tempPx;repmat(mean_param_value_discr(pt,c),numel(tempRx{c}),1)];
                end
                tempRx=cat(1,tempRx{:});
                if ~isempty(tempRx)
                    [corr_TP_split(n,1,pt),corr_TP_split(n,2,pt)]=corr(tempPx,tempRx,'type',ty);
                end
                if pt==1 && n==selc{1}(6)
figure('Position',[270   849   560   420])
scatter(tempPx+rand(size(tempPx))/5000,tempRx)
title(sprintf('corr=%.3f',corr_TP_split(n,1,pt)))
                end
            end

            tempRx=indiv_split_response(:,:,n);
            tempPx=[];
            for c=1:numel(splits2{S})-1
                tempRx{1,c}=cat(1,tempRx{:,c});
                tempRx{2,c}=[];
                tempPx=[tempPx;repmat(mean_param_value_discr(3,c),numel(tempRx{1,c}),1)];
            end
            tempRx=cat(1,tempRx{:});
            if ~isempty(tempRx)
                [corr_TP_split(n,1,3),corr_TP_split(n,2,3)]=corr(tempPx,tempRx,'type','Spearman');
            end


%             for c=1:numel(splits2{S})-1
%                 if ~isempty(indiv_split_response{1,c,n}) && ~isempty(indiv_split_response{2,c,n})
%                     [~,ks_p_TP_split(n,c)]=kstest2(indiv_split_response{1,c,n},indiv_split_response{2,c,n});
%                 end
%             end
%             if ~isempty(indiv_split_response{1,1,n}) && ~isempty(indiv_split_response{1,end,n})
%                 [~,ks_p_TP13_split(n,1)]=kstest2(indiv_split_response{1,1,n},indiv_split_response{1,end,n});
%             end
%             if ~isempty(indiv_split_response{2,1,n}) && ~isempty(indiv_split_response{2,end,n})
%                 [~,ks_p_TP13_split(n,2)]=kstest2(indiv_split_response{2,1,n},indiv_split_response{2,end,n});
%             end
        end
    end
    %
    corr_TP_split2=cat(2,corr_TP_split(:,:,1),corr_TP_split(:,:,2),corr_TP_split(:,:,3));
    corrV=corr_TP_split2(selc{1},:);
    corrP=corr_TP_split2(selc{2},:);
    %
    figure
    scatter(1,corrV(corrV(:,2)<.05,1),'r*');hold on;
    scatter(1,corrV(corrV(:,2)>=.05,1),'ro');hold on;
    scatter(2,corrV(corrV(:,4)<.05,3),'r*');hold on;
    scatter(2,corrV(corrV(:,4)>=.05,3),'ro');hold on;
    scatter(3,corrV(corrV(:,6)<.05,5),'r*');hold on;
    scatter(3,corrV(corrV(:,6)>=.05,5),'ro');hold on;
    scatter(4,corrP(corrP(:,2)<.05,1),'b*');hold on;
    scatter(4,corrP(corrP(:,2)>=.05,1),'bo');hold on;
    scatter(5,corrP(corrP(:,4)<.05,3),'b*');hold on;
    scatter(5,corrP(corrP(:,4)>=.05,3),'bo');hold on;
    scatter(6,corrP(corrP(:,6)<.05,5),'b*');hold on;
    scatter(6,corrP(corrP(:,6)>=.05,5),'bo');hold on;
    xlim([.5 6.5])
    % figure
    line([.5 6.5],[0 0],'Color','k');hold on
    plot(1:3,corrV(:,[1 3 5]),'r');hold on
    plot(4:6,corrP(:,[1 3 5]),'b');hold on
    title(num2str(S));

    %
    figure
    b=-1:.2:1;
    tt=tiledlayout(2,2);
    title(tt,split_var_name)
    nexttile;
    histogram(corrV(:,1),b,'FaceColor',[.7 .7 .7]);hold on;
    histogram(corrV(corrV(:,2)<.05,1),b,'FaceColor',[.2 .2 .2],'FaceAlpha',.5);hold on;
    title('correlation vpm puff')
    nexttile;
    histogram(corrV(:,3),b,'FaceColor',[.7 .7 .7]);hold on;
    histogram(corrV(corrV(:,4)<.05,3),b,'FaceColor',[.2 .2 .2],'FaceAlpha',.5);hold on;
    title('correlation vpm touch')
    nexttile;
    histogram(corrP(:,1),b,'FaceColor',[.7 .7 .7]);hold on;
    histogram(corrP(corrP(:,2)<.05,1),b,'FaceColor',[.2 .2 .2],'FaceAlpha',.5);hold on;
    title('correlation pom puff')
    nexttile;
    histogram(corrP(:,3),b,'FaceColor',[.7 .7 .7]);hold on;
    histogram(corrP(corrP(:,4)<.05,3),b,'FaceColor',[.2 .2 .2],'FaceAlpha',.5);hold on;
    title('correlation pom touch')

     %%
    figp=figure('Position',[1000         781         754         557]);
    pcorrV=[sum(corrV(:,[1 3])<0 & corrV(:,[2 4])<0.05);sum(corrV(:,[1 3])>=0 & corrV(:,[2 4])<0.05);sum(corrV(:,[2 4])>=0.05)];
     pcorrP=[sum(corrP(:,[1 3])<0 & corrP(:,[2 4])<0.05);sum(corrP(:,[1 3])>=0 & corrP(:,[2 4])<0.05);sum(corrP(:,[2 4])>=0.05)];
    
    tt=tiledlayout(2,2);
    title(tt,split_var_name)
    nexttile;
    pie(pcorrV(:,1))
    title('correlation vpm puff')
    nexttile;
    pie(pcorrV(:,2))
    title('correlation vpm touch')
    nexttile;
    pie(pcorrP(:,1))
    title('correlation pom puff')
    nexttile;
    pie(pcorrP(:,2))
    title('correlation pom touch')
    lgd = legend({'sig. anticorrelated','sig. correlated','nonsig.'});
    lgd.Layout.Tile = 'east';
    exportgraphics(figp,[figdir 'test_trial_corr_' split_var_name '_' datestr(now,'yymmdd') '.pdf'],'BackgroundColor','none','ContentType','vector')

end
%%
cp=cellfun(@mean,indiv_split_response);
cpp=cp(:,:,selc{2})./.05;
cppm=mean(cpp,3,'omitnan');
%%
ks13c=ks_p_TP13_split(selc{1},:)<.05
%%
figure
% histogram(Pt{1,2},'Normalization','cdf');hold on;
% histogram(Pt{2,2},'Normalization','cdf');hold on;
ecdf(Pt{1,3});hold on;
ecdf(Pt{2,3});hold on;
[H,p]=kstest2(Pt{1,3},Pt{2,3})

%%
sb=ix(RasterHP{n}(:,2));
RasterHP{n}(:,3)=sb';

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