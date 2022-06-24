%%
include_wh_params={'Angle';'Velocity';'Acceleration';'Curvature'};
W_dir='F:\Whisker\Wdec\';
kin_window=[-20 -5;0 15];
resp_window=[5 40];
min_trial=10;
[Rkin,Wkin,TTkin]=get_kinematic_comparison(trigsP,trigsT,DiscreteData,W_dir,kin_window,resp_window,min_trial,ppms,include_wh_params,false);
%
WA=cat(1,Wkin{:});
Wp=prctile(WA,[2.5;97.5],1);
binN=7;
Wb=nan(binN,size(WA,2));
for d=1:size(WA,2)
    Wb(:,d)=linspace(Wp(1,d),Wp(2,d),binN)';
end
Wb_plot=Wb(1:end-1,:)+diff(Wb(1:2,:),1,1);
WkinB=cell(size(Wkin));
for r=1:numel(Wkin)
    if ~isempty(Wkin{r})
       WkinB{r}=nan(size(Wkin{r}));
        for d=1:size(WA,2)
       WkinB{r}(:,d)=discretize(Wkin{r}(:,d),Wb(:,d));    
       end
    end
end
%

RB=cell(size(Wkin));
for r=1:numel(Wkin)
    if ~isempty(Wkin{r})
        RB{r}=nan(binN-1,2,size(WA,2));
        for d=1:size(WA,2)
            temp=cat(2,WkinB{r}(:,d),TTkin{r},Rkin{r});
            temp(isnan(temp(:,1)) | temp(:,1)==0,:)=[];
            temp2=accumarray(temp(:,1:end-1),temp(:,end),[binN-1 2],@mean,nan);
            N2=accumarray(temp(:,1:end-1),1,[binN-1 2],@sum,0);
            temp2(N2<3)=nan;
             RB{r}(:,:,d)=temp2;
        end
    end
end
%%
close all
RBp=RB(pomr);
for d=1:size(WA,2)
    figure('Position',[100        -220        2210        1028]);
    t=tiledlayout('flow');
    title(t,['POm ' include_wh_params{d}])
    for r=1:numel(RBp)
        if ~isempty(RBp{r})
            nexttile;
            plot(Wb_plot(:,d),RBp{r}(:,:,d));box off;
            xlim(Wb_plot([1 end],d))
            title(sprintf('cell %u',pomr(r)))
            ylabel('#sp resp.')
            xlabel('kin.value')
        end
    end
    legend('puff','touch'); legend('boxoff');
   %export_fig([figdir 'POm_binned_kinematic_' include_wh_params{d} '.pdf'],'-native')
end
%
RBp=RB(vpmr);
for d=1:size(WA,2)
    figure('Position',[100        -220        2210        1028]);
    t=tiledlayout('flow');
    title(t,['VPM ' include_wh_params{d}])
    for r=1:numel(RBp)
        if ~isempty(RBp{r})
            nexttile;
            plot(Wb_plot(:,d),RBp{r}(:,:,d));box off;
            title(sprintf('cell %u',vpmr(r)))
            ylabel('#sp resp.')
            xlabel('kin.value')
            xlim(Wb_plot([1 end],d))
        end
    end
        legend('puff','touch'); legend('boxoff')
     %   export_fig([figdir 'POm_binned_kinematic_' include_wh_params{d} '.pdf'],'-native')
end
%%
close all
RBp=mean(cat(4,RB{pomr}),4,'omitnan');
figure('Position',[110         193        1707         683]);
t=tiledlayout('flow');
title(t,'POm')
for d=1:size(WA,2)
    nexttile;
    
    
    
    plot(Wb_plot(:,d),RBp(:,:,d));box off;
    xlim(Wb_plot([1 end],d))
    %title(sprintf('cell %u',pomr(r)))
    ylabel('#sp resp.')
    xlabel('kin.value')
    title([include_wh_params{d}])
    legend('puff','touch'); legend('boxoff');
   
end
 export_fig([figdir 'POm_mean_binned_kinematic.pdf'],'-native')
%
%close all
RBp=mean(cat(4,RB{vpmr}),4,'omitnan');
figure('Position',[110         193        1707         683]);
t=tiledlayout('flow');
title(t,'VPM')
for d=1:size(WA,2)
    nexttile;
    
    
    
    plot(Wb_plot(:,d),RBp(:,:,d));box off;
    xlim(Wb_plot([1 end],d))
    %title(sprintf('cell %u',pomr(r)))
    ylabel('#sp resp.')
    xlabel('kin.value')
    title([include_wh_params{d}])
    legend('puff','touch'); legend('boxoff');
   
end
export_fig([figdir 'VPM_mean_binned_kinematic.pdf'],'-native')
%%
close all
RBp=RB(pomr);
for d=1:size(WA,2)
    figure('Position',[100        -220        2210        1028]);
    t=tiledlayout('flow');
    title(t,['POm ' include_wh_params{d}])
    for r=1:numel(RBp)
        if ~isempty(RBp{r})
            nexttile;
            Rtemp=Rkin{pomr(r)}+rand(size(Rkin{pomr(r)},1),1)/3;
                scatter(Wkin{pomr(r)}(TTkin{pomr(r)}==1,d),Rtemp(TTkin{pomr(r)}==1,1),15,'b');hold on
                scatter(Wkin{pomr(r)}(TTkin{pomr(r)}==2,d),Rtemp(TTkin{pomr(r)}==2,1),15,'r');hold off
                ylabel('Response spikes + noise')
                xlabel('kinematic val.')
            
            xlim(Wb_plot([1 end],d))
            title(sprintf('cell %u',pomr(r)))
           
        end
    end
    legend('puff','touch'); legend('boxoff');
  % export_fig([figdir 'POm_scatter_kinematic_' include_wh_params{d} '.pdf'],'-native')
end
%
RBp=RB(vpmr);
for d=1:size(WA,2)
    figure('Position',[100        -220        2210        1028]);
    t=tiledlayout('flow');
    title(t,['VPM ' include_wh_params{d}])
    for r=1:numel(RBp)
        if ~isempty(RBp{r})
            nexttile;
            Rtemp=Rkin{vpmr(r)}+rand(size(Rkin{vpmr(r)},1),1)/3;
                scatter(Wkin{vpmr(r)}(TTkin{vpmr(r)}==1,d),Rtemp(TTkin{vpmr(r)}==1,1),15,'b');hold on
                scatter(Wkin{vpmr(r)}(TTkin{vpmr(r)}==2,d),Rtemp(TTkin{vpmr(r)}==2,1),15,'r');hold off
                ylabel('Response spikes + noise')
                xlabel('kinematic val.')
            
            xlim(Wb_plot([1 end],d))
            title(sprintf('cell %u',vpmr(r)))
           
        end
    end
    legend('puff','touch'); legend('boxoff');
 %  export_fig([figdir 'VPM_scatter_kinematic_' include_wh_params{d} '.pdf'],'-native')
end