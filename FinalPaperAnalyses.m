%% prepare
DataPath='D:\active_passive_touch_data\';
ppms=20;

load([DataPath 'ephys_database.mat'])
load([DataPath 'DiscreteData.mat'])
load([DataPath 'WhiskerBehavData.mat'])
minFillQ=4;
%
figdir=['C:\Data\Whisker\Paperfigs_matlab\' datestr(now,'YYmmDD') '_50_excllowtouch\'];
set(0,'defaultfigurecolor',[1 1 1]);
set(0, 'DefaultAxesBox', 'off')
set(0,'defaultAxesFontName', 'Arial')
set(0,'defaultTextFontName', 'Arial')
write_figs=true;
if write_figs
    mkdir(figdir)
end
% settings
psth_time_before_trig=25;%ms
psth_time_after_trig=75;%ms
psth_binsize=1;%ms
reassign_response_types=true;
psth_safety_margin=psth_time_after_trig;%ms
rate_time_before_trig=50;%ms
rate_time_after_trig=rate_time_before_trig;%ms
rate_binsize=rate_time_before_trig;%ms
sig=.05;
isi=10;
puff_triglength_limit=[28 32];%ms
touch_triglength_limit=[0 inf];%ms
phasebinN=32;
amp_thresh=3;

rate_bins=-rate_time_before_trig:rate_binsize:rate_time_after_trig;
psth_bins=-psth_time_before_trig:psth_binsize:psth_time_after_trig;
psth_relativeTime=psth_bins(1:end-1)+diff(psth_bins);


% calculate supporting data
clear MeanP MeanPL MeanT
%
%Puff
parameter='Puff';exclude={'Pole';'Light';'Exclude';'Grooming'};include=[];
if ~exist('MeanP','var') 
    [~,~,~,~,trigsP,zeta_psth_P]=get_PSTH_pop(parameter, DiscreteData,puff_triglength_limit,ppms,psth_time_before_trig,psth_time_after_trig,exclude,psth_safety_margin,psth_binsize,isi,'none',[]);
    [MeanPA,trialPA]=get_trig_cont_pop(trigsP,[-50:0],Amp,20,DiscreteData);
    pA=cellfun(@(x) mean(x,2),trialPA,'UniformOutput',0);
    trigsPW=cellfun(@(x,y) x(y>3),trigsP,pA,'UniformOutput',0);
    trigsPnW=cellfun(@(x,y) x(y<=3),trigsP,pA,'UniformOutput',0);
    [PuffV]=mk_vector_from_start_length(DiscreteData,'Puff');
    [MeanWp,TrialWp]=get_trig_cont_pop(trigsP,psth_relativeTime,Angle,20,DiscreteData);
    SEM_Wp=cellfun(@(x) std(x,0,1)./sqrt(size(x,1)),TrialWp,'UniformOutput',false);
    [MeanP]=get_trig_cont_pop(trigsP,psth_relativeTime,PuffV,1,DiscreteData);
    [MeanAp,TrialAp]=get_trig_cont_pop(trigsP,psth_relativeTime,Amp,20,DiscreteData);
    [MeanWpw,TrialWpw]=get_trig_cont_pop(trigsPW,psth_relativeTime,Angle,20,DiscreteData);
    SEM_Wpw=cellfun(@(x) std(x,0,1)./sqrt(size(x,1)),TrialWpw,'UniformOutput',false);
    [MeanWpnw,TrialWpnw]=get_trig_cont_pop(trigsPnW,psth_relativeTime,Angle,20,DiscreteData);
    SEM_Wpnw=cellfun(@(x) std(x,0,1)./sqrt(size(x,1)),TrialWpnw,'UniformOutput',false);
    [MeanApw]=get_trig_cont_pop(trigsPW,psth_relativeTime,Amp,20,DiscreteData);
    [MeanApnw]=get_trig_cont_pop(trigsPnW,psth_relativeTime,Amp,20,DiscreteData);
    [MeanPW]=get_trig_cont_pop(trigsPW,psth_relativeTime,PuffV,1,DiscreteData);
    [MeanPnW]=get_trig_cont_pop(trigsPnW,psth_relativeTime,PuffV,1,DiscreteData);
    
end
%Puff Light
parameter='Puff';exclude={'Pole';'Exclude';'Grooming'};include={{'Light'};[-5 35]};
if ~exist('MeanPL','var') 
    [~,~,~,~,trigsPL]=get_PSTH_pop(parameter, DiscreteData,puff_triglength_limit,ppms,psth_time_before_trig,psth_time_after_trig,exclude,psth_safety_margin,psth_binsize,isi,'none',[],include{1},include{2});
    [MeanWpL,TrialWpL]=get_trig_cont_pop(trigsPL,psth_relativeTime,Angle,20,DiscreteData);
    SEM_WpL=cellfun(@(x) std(x,0,1)./sqrt(size(x,1)),TrialWpL,'UniformOutput',false);
    [PuffV]=mk_vector_from_start_length(DiscreteData,'Puff');
    [MeanPL]=get_trig_cont_pop(trigsPL,psth_relativeTime,PuffV,1,DiscreteData);
    [MeanL,trialL]=get_trig_cont_pop(trigsPL,psth_relativeTime,'Light',1,DiscreteData);
    [MeanAL,TrialAL]=get_trig_cont_pop(trigsPL,psth_relativeTime,Amp,20,DiscreteData);
    hasLight=~cellfun(@isempty,MeanPL);
end
% Touch
parameter='Touch';exclude={'Puff';'Light';'Exclude';'Grooming'};include=[];
if ~exist('MeanT','var') 
    [~,~,~,~,trigsT,zeta_psth_T,trigsT_l]=get_PSTH_pop(parameter, DiscreteData,touch_triglength_limit,ppms,psth_time_before_trig,psth_time_after_trig,exclude,psth_safety_margin,psth_binsize,isi,'none',[]);
    [MeanWt,TrialWt]=get_trig_cont_pop(trigsT,psth_relativeTime,Angle,20,DiscreteData);
    SEM_Wt=cellfun(@(x) std(x,0,1)./sqrt(size(x,1)),TrialWt,'UniformOutput',false);
    [MeanT]=get_trig_cont_pop(trigsT,psth_relativeTime,'Touch',1,DiscreteData);
    [MeanAt,TrialAt]=get_trig_cont_pop(trigsT,psth_relativeTime,Amp,20,DiscreteData);
end

%puff rates
parameter='Puff';
exclude={'Pole';'Light';'Exclude';'Grooming'};include=[];
[H_Puff,~,p_Puff,~,p_trigs,~,p_trigL]=get_PSTH_pop(parameter, DiscreteData,puff_triglength_limit,ppms,rate_time_before_trig,rate_time_after_trig,exclude,psth_safety_margin,rate_binsize,isi,'left',[]);
p_Puff=p_Puff{1}(:,2);
sig_Puff=p_Puff<sig;
H_P=cat(2,H_Puff{:})';

%touch rates
parameter='Touch';Name='Touch';
exclude={'Puff';'Light';'Exclude';'Grooming'};
[~,~,~,~,ttrigs,zeta_T,t_trigL]=get_PSTH_pop(parameter, DiscreteData,touch_triglength_limit,ppms,rate_time_before_trig,rate_time_after_trig,exclude,psth_safety_margin,rate_binsize,isi,'left',[]);
[H_Touch,~,p_Touch]=get_PSTH_pop(parameter, DiscreteData,touch_triglength_limit,ppms,rate_time_before_trig,rate_time_after_trig,exclude,psth_safety_margin,rate_binsize,isi,'left',[]);
p_Touch=p_Touch{1}(:,2);
sig_Touch=p_Touch<sig;
H_T=cat(2,H_Touch{:})';
ntouches=cellfun(@numel,trigsT);
%
has_Touch=ntouches>=10;%find(~isnan(p_Touch));%

if reassign_response_types
    typemat=false(numel(p_Touch),2+2+2+2);
 vpm=find(RecDB.DefinedNucleus=='VPM' & RecDB.FillQuality<=minFillQ & RecDB.RecordingType=='juxta');
 pom=find(RecDB.DefinedNucleus=='POm' & RecDB.FillQuality<=minFillQ & RecDB.RecordingType=='juxta');
typemat(vpm,1)=true;
typemat(pom,2)=true;
typemat(has_Touch,3)=true;
typemat(hasLight,4)=true;
typemat(sig_Puff,5)=true;
typemat(sig_Touch,6)=true;
typemat(zeta_output.p_Puff<sig,7)=true;
typemat(zeta_output.p_Touch<sig,8)=true;

vpmr=find(all(typemat(:,[1 3 5]),2));
pomr=find(all(typemat(:,[2 3 5]),2));


vpmnpr=find(all(cat(2,typemat(:,[1 3]),~typemat(:,5)),2));
pomnpr=find(all(cat(2,typemat(:,[2 3]),~typemat(:,5)),2));

vpmra=find(all(typemat(:,[1 5]),2));
pomra=find(all(typemat(:,[2 5]),2));

vpmz=find(all(typemat(:,[1 3 7]),2));
pomz=find(all(typemat(:,[2 3 7]),2));

vpmza=find(all(typemat(:,[1 7]),2));
pomza=find(all(typemat(:,[2 7]),2));

vpmrz=find(all(typemat(:,[1 3 5 7]),2));
pomrz=find(all(typemat(:,[2 3 5 7]),2));

vpmrza=find(all(typemat(:,[1 5 7]),2));
pomrza=find(all(typemat(:,[2 5 7]),2));

vpmL=find(all(typemat(:,[1 4]),2));
pomL=find(all(typemat(:,[2 4]),2));

vpmLa=find(all(typemat(:,[1 4 5]),2));
pomLa=find(all(typemat(:,[2 4 5]),2));

vpmLr=find(all(typemat(:,[1 3 4 5]),2));
pomLr=find(all(typemat(:,[2 3 4 5]),2));
end
fprintf('test interval: %u ms \n vpmr:%u, vpmra:%u , vpmLa:%u \n pomr:%u, pomra:%u , pomLa:%u \n',rate_time_after_trig,numel(vpmr),numel(vpmra),numel(vpmL),numel(pomr),numel(pomra),numel(pomL))
VN={'condition 1','condition 2','VPM cell IDs','POm cell IDs', 'VPM N','VPM N Animals','POm N','POm N Animals','Total N Animals',...
    'VPM mean condition 1 base','VPM sem condition 1 base','VPM mean condition 1 resp','VPM sem condition 1 resp','VPM p-value condition 1 response vs baseline',...
    'VPM mean condition 2 base','VPM sem condition 2 base','VPM mean condition 2 resp','VPM sem condition 2 resp','VPM p-value condition 2 response vs baseline',...
    'VPM p-value baseline across conditions','VPM p-value response across conditions',...
    'POm mean condition 1 base','POm sem condition 1 base','POm mean condition 1 resp','POm sem condition 1 resp','POm p-value condition 1 response vs baseline',...
    'POm mean condition 2 base','POm sem condition 2 base','POm mean condition 2 resp','POm sem condition 2 resp','POm p-value condition 2 response vs baseline',...
    'POm p-value baseline across conditions','POm p-value response across conditions',...
    'VPM mean modulation condition 1','VPM sem modulation condition 1','VPM p-value modulation cond1 vs zero','VPM mean modulation condition 2','VPM sem modulation condition 2','VPM p-value modulation cond2 vs zero','VPM p-value modulation across conditions',...
    'POm mean modulation condition 1','POm sem modulation condition 1','POm p-value modulation cond1 vs zero','POm mean modulation condition 2','POm sem modulation condition 2','POm p-value modulation cond2 vs zero','POm p-value modulation across conditions',...
    'p-value modulation condition 1 across nuclei', 'p-value modulation condition 2 across nuclei',...
    'VPM mean fsl condition 1','VPM sem fsl condition 1','VPM mean fsl condition 2','VPM sem fsl condition 2','VPM p-value fsl across conditions',...
    'POm mean fsl condition 1','POm sem fsl condition 1','POm mean fsl condition 2','POm sem fsl condition 2','POm p-value fsl across conditions',...
    'p-value fsl condition 1 across nuclei', 'p-value fsl condition 2 across nuclei'}';


SumStats=struct();%cell2table(cat(2,{''},{''},{nan(4,1)},{nan(4,1)},num2cell(nan(1,numel(VN)-4))),'VariableNames',VN);
%%
A={'N neurons'
'N animals'
'Baseline mean'
'Baseline SEM'
'p-Value Baseline across'
'Response mean'
'Response SEM'
'p-Value Response across'
'p-Value Response  = Baseline'
'Mean modulation'
'Mean modulation SEM'
'p-Value modulation = zero'
'p-Value modulation across'
'First spike latency mean'
'First spike latency SEM'
'p-Value First spike latency across'
'N neurons'
'N animals'
'Baseline mean'
'Baseline SEM'
'p-Value Baseline across'
'Response mean'
'Response SEM'
'p-Value Response across'
'p-Value Response  = Baseline'
'Mean modulation'
'Mean modulation SEM'
'p-Value modulation = zero'
'p-Value modulation across'
'First spike latency mean'
'First spike latency SEM'
'p-Value First spike latency across'
'N animals total'
'p-Value Baseline'
'p-Value Response'
'p-Value modulation'
'p-Value First spike latency'};
A=cat(2,cat(1,repmat({'VPM'},16,1),repmat({'POm'},16,1),repmat({'VPM vs. POm'},5,1)),A,num2cell(nan(37,10)));

SumStatTable=cell2table(A,'VariableNames',{'Nucleus','ValName','Puff','Touch','Q-Puff','W-Puff','Quiescence spont','Whisking spont','nL-Puff','L-Puff','nL-Whisking','L-Whisking'});


%% Figure 2
Fig_no='Fig2_';
SumStatLine=1;

SumStats(SumStatLine).condition_1='Puff';
SumStats(SumStatLine).condition_2='Touch';
SumStats(SumStatLine).VPM_cell_IDs=vpmr';
SumStats(SumStatLine).POm_cell_IDs=pomr';
[N_animals]=get_N_animals(SumStats(SumStatLine).VPM_cell_IDs,SumStats(SumStatLine).POm_cell_IDs,RecDB);
SumStats(SumStatLine).VPM_N=numel(vpmr);
SumStats(SumStatLine).VPM_N_animals=N_animals(3);
SumStats(SumStatLine).POm_N=numel(pomr);
SumStats(SumStatLine).POm_N_animals=N_animals(4);
SumStats(SumStatLine).Total_N_animals=N_animals(1);
C=3;
SumStatTable{1,C}=SumStats(SumStatLine).VPM_N;
SumStatTable{2,C}=SumStats(SumStatLine).VPM_N_animals;
SumStatTable{17,C}=SumStats(SumStatLine).POm_N;
SumStatTable{18,C}=SumStats(SumStatLine).POm_N_animals;
SumStatTable{33,C}=SumStats(SumStatLine).Total_N_animals;

%puff psths
parameter='Puff';
contPlot={MeanWp,MeanP,MeanAp};
exclude={'Pole';'Light';'Exclude';'Grooming'};include=[];
Name='Puff';
figure('Position',[38 260 1413 525],'Color','White');
[P_PSTH_instrate]=plot_Pop_PSTH_wrapper3(parameter, DiscreteData,puff_triglength_limit,ppms,psth_time_before_trig,psth_time_after_trig,exclude,psth_safety_margin,psth_binsize,[],Name,{vpmr;pomr},{'VPM';'POm'},include,'latency',-inf,contPlot,{'meanbefore','scaleminmax','mono'},'baseline corr. whisker angle');
if write_figs;exportgraphics(gcf,[figdir Fig_no 'Pop_' parameter '_instRateResp_' datestr(now,'yymmdd') '.pdf'],'BackgroundColor','none');end
fig2a2=figure('Position',[-1731         836        1322         171],'Color','White');
Name=strcat(parameter,' trig whisking');
plot_wh_trace_with_sem(MeanWp,SEM_Wp,{vpmr;pomr},{'VPM';'POm'},'wh angle',psth_relativeTime,Name)
if write_figs;exportgraphics(fig2a2,[figdir Fig_no 'Pop_' parameter '_wh_with_sem_' datestr(now,'yymmdd') '.pdf'],'BackgroundColor','none');end


%touch psths
parameter='Touch';Name='Touch';
exclude={'Puff';'Light';'Exclude';'Grooming'};
contPlot={MeanWt,MeanT,MeanAt};
figure('Position',[38 260 1413 525],'Color','White');
[T_PSTH_instrate]=plot_Pop_PSTH_wrapper3(parameter, DiscreteData,touch_triglength_limit,ppms,psth_time_before_trig,psth_time_after_trig,exclude,psth_safety_margin,psth_binsize,[],parameter,{vpmr;pomr},{'VPM';'POm'},include,'latency',-inf,contPlot,{'meanbefore','scaleminmax','mono'},'baseline corr. whisker angle');
if write_figs;exportgraphics(gcf,[figdir Fig_no 'Pop_' parameter '_instRateResp_' datestr(now,'yymmdd') '.pdf'],'BackgroundColor','none','ContentType','vector');end
fig2a2=figure('Position',[-1731         836        1322         171],'Color','White');
Name=strcat(parameter,' trig whisking');
plot_wh_trace_with_sem(MeanWt,SEM_Wt,{vpmr;pomr},{'VPM';'POm'},'wh angle',psth_relativeTime,Name)
if write_figs;exportgraphics(fig2a2,[figdir Fig_no 'Pop_' parameter '_wh_with_sem_' datestr(now,'yymmdd') '.pdf'],'BackgroundColor','none');end



% summary puff / touch
figure('Position',[910   176   500   588])
[pI,pC,mS]=plot_ratecomp3(cat(3,H_P,H_T),cat(2,sig_Puff,sig_Touch),{'0';'1';'0';'1'},vpmr,{'Puff','Touch'},'mean','Rate [Hz]');
title('Puff/Touch Response VPM', 'resp cells w touch')
if write_figs;exportgraphics(gcf,[figdir Fig_no 'RatesPT_VPM_RespondingCellsWithTouch.pdf'], 'BackgroundColor','none','ContentType','vector');end

SumStats(SumStatLine).VPM_mean_cond1_base=mS(1,1,1);
SumStats(SumStatLine).VPM_sem_cond1_base=mS(1,1,2);
SumStats(SumStatLine).VPM_mean_cond1_resp=mS(1,2,1);
SumStats(SumStatLine).VPM_sem_cond1_resp=mS(1,2,2);
SumStats(SumStatLine).VPM_p_cond1=pI(1);
SumStats(SumStatLine).VPM_mean_cond2_base=mS(2,1,1);
SumStats(SumStatLine).VPM_sem_cond2_base=mS(2,1,2);
SumStats(SumStatLine).VPM_mean_cond2_resp=mS(2,2,1);
SumStats(SumStatLine).VPM_sem_cond2_resp=mS(2,2,2);
SumStats(SumStatLine).VPM_p_cond2=pI(2);
SumStats(SumStatLine).VPM_p_base_cond12=pC(1);
SumStats(SumStatLine).VPM_p_resp_cond12=pC(2);




%
figure('Position',[910   176   500   588])
[pI,pC,mS]=plot_ratecomp3(cat(3,H_P,H_T),cat(2,sig_Puff,sig_Touch),{'0';'1';'0';'1'},pomr,{'Puff','Touch'},'mean','Rate [Hz]');
title('Puff/Touch Response POm', 'resp cells w touch')
if write_figs;exportgraphics(gcf,[figdir Fig_no 'RatesPT_POm_RespondingCellsWithTouch.pdf'], 'BackgroundColor','none','ContentType','vector');end

SumStats(SumStatLine).POm_mean_cond1_base=mS(1,1,1);
SumStats(SumStatLine).POm_sem_cond1_base=mS(1,1,2);
SumStats(SumStatLine).POm_mean_cond1_resp=mS(1,2,1);
SumStats(SumStatLine).POm_sem_cond1_resp=mS(1,2,2);
SumStats(SumStatLine).POm_p_cond1=pI(1);
SumStats(SumStatLine).POm_mean_cond2_base=mS(2,1,1);
SumStats(SumStatLine).POm_sem_cond2_base=mS(2,1,2);
SumStats(SumStatLine).POm_mean_cond2_resp=mS(2,2,1);
SumStats(SumStatLine).POm_sem_cond2_resp=mS(2,2,2);
SumStats(SumStatLine).POm_p_cond2=pI(2);
SumStats(SumStatLine).POm_p_base_cond12=pC(1);
SumStats(SumStatLine).POm_p_resp_cond12=pC(2);



pRc=nan(2,2);
for c=1:2
pRc(1,c)=ranksum(H_P(vpmr,c,1),H_P(pomr,c,1));
end
for c=1:2
pRc(2,c)=ranksum(H_T(vpmr,c,1),H_T(pomr,c,1));
end
%
% touch puff comp
p_PuffPTr=sig_test_comp_trig({p_trigs;ttrigs},DiscreteData,rate_bins,2,'both',ppms);
p_PuffPTb=sig_test_comp_trig({p_trigs;ttrigs},DiscreteData,rate_bins,1,'both',ppms);
sig_PTr=p_PuffPTr<sig;
sig_PTb=p_PuffPTb<sig;
tempD=log2(cat(2,H_P(:,2)./H_P(:,1),H_T(:,2)./H_T(:,1)));
figure('Position',[ 910   176   500   588])
plot_ratecomp2(tempD,sig_PTr,{'Puff','Touch'},{vpmr;pomr},{'VPM';'POm'},'mean','log2 Rate fold change');title('Touch/Puff fold change Response, responsive cells');
set(gca,'YTick',log2([0.5 1:10 12:2:20 24:4:64]))
set(gca,'YTickLabel',2.^get(gca,'YTick'))
if write_figs;exportgraphics (gcf,[figdir Fig_no 'foldchangeRatesTouchVsPuffRespondingCells.pdf'], 'BackgroundColor','none','ContentType','vector');end

tempD=cat(2,mod_depth(H_P),mod_depth(H_T));
figure('Position',[ 910   176   500   588])
[pI,pC,mS,p0]=plot_ratecomp2(tempD,sig_PTr,{'Puff','Touch'},{vpmr;pomr},{'VPM';'POm'},'mean','mod depth');title('Touch/Puff modulation, responsive cells')

 if write_figs;exportgraphics (gcf,[figdir Fig_no 'modulationTouchVsPuffRespondingCellsLog.pdf'], 'BackgroundColor','none','ContentType','vector');end


tempD=cat(2,H_P(:,2)./H_P(:,1),H_T(:,2)./H_T(:,1));
figure('Position',[ 910   176   500   588])
plot_ratecomp2(tempD,sig_PTr,{'Puff','Touch'},{vpmr;pomr},{'VPM';'POm'},'mean','Rate fold change');title('Touch/Puff fold change Response, responsive cells')
 set(gca,'YScale','log')
 set(gca,'YTick',[0.5 1:10 12:2:20 25 30 35])
 if write_figs;exportgraphics (gcf,[figdir Fig_no 'foldchangeRatesTouchVsPuffRespondingCellsLog.pdf'], 'BackgroundColor','none','ContentType','vector');end

SumStats(SumStatLine).VPM_mean_modulation_cond1=mS(1,1,1);
SumStats(SumStatLine).VPM_sem_modulation_cond1=mS(1,1,2);
SumStats(SumStatLine).VPM_p_modulation_cond10=p0(1,1);
SumStats(SumStatLine).VPM_mean_modulation_cond2=mS(1,2,1);
SumStats(SumStatLine).VPM_sem_modulation_cond2=mS(1,2,2);
SumStats(SumStatLine).VPM_p_modulation_cond20=p0(1,2);
SumStats(SumStatLine).VPM_p_modulation_cond12=pI(1);
SumStats(SumStatLine).POm_mean_modulation_cond1=mS(2,1,1);
SumStats(SumStatLine).POm_sem_modulation_cond1=mS(2,1,2);
SumStats(SumStatLine).POm_p_modulation_cond10=p0(2,1);
SumStats(SumStatLine).POm_mean_modulation_cond2=mS(2,2,1);
SumStats(SumStatLine).POm_sem_modulation_cond2=mS(2,2,2);
SumStats(SumStatLine).POm_p_modulation_cond20=p0(2,2);
SumStats(SumStatLine).POm_p_modulation_cond12=pI(2);
SumStats(SumStatLine).VPMPOm_p_modulation_cond1=pC(1);
SumStats(SumStatLine).VPMPOm_p_modulation_cond2=pC(2);
SumStats(SumStatLine).VPMPOm_p_cond1_base=pRc(1,1);
SumStats(SumStatLine).VPMPOm_p_cond2_base=pRc(2,1);
SumStats(SumStatLine).VPMPOm_p_cond1_resp=pRc(1,2);
SumStats(SumStatLine).VPMPOm_p_cond2_resp=pRc(2,2);


% summary no puff / touch
figure('Position',[265         633        1220         587])
tt=tiledlayout(1,3);
title(tt,'non-puff-responsive neurons')
nexttile(1);
[pI,pC,mS]=plot_ratecomp3(cat(3,H_P,H_T),cat(2,sig_Puff,sig_Touch),{'0';'1';'0';'1'},vpmnpr,{'Puff','Touch'},'mean','Rate [Hz]');
title('Puff/Touch Response VPM')
nexttile(2);
[pI,pC,mS]=plot_ratecomp3(cat(3,H_P,H_T),cat(2,sig_Puff,sig_Touch),{'0';'1';'0';'1'},pomnpr,{'Puff','Touch'},'mean','Rate [Hz]');
title('Puff/Touch Response POm')
nexttile(3);
tempD=cat(2,mod_depth(H_P),mod_depth(H_T));
[pI,pC,mS,p0]=plot_ratecomp2(tempD,true(size(sig_PTr)),{'Puff','Touch'},{vpmnpr;pomnpr},{'VPM';'POm'},'mean','mod depth');
title('Touch/Puff modulation')
 if write_figs;exportgraphics (gcf,[figdir 'SuppNonPuffRespTouch.pdf'], 'BackgroundColor','none','ContentType','vector');end


%%
% Supp Figure 2
Fig_no='EDFig2_';
lat_binsize=1;
lat_time_before_trig=0;
lat_time_after_trig=75;
lat_safety_margin=lat_time_after_trig;
parameter='Puff';
exclude={'Pole';'Light';'Exclude';'Grooming'};
triglength_limit=[28 32];%ms
tempDD=struct2cell(DiscreteData);
spikes=squeeze(tempDD(1,1,:));clear temp*
[RasterHPL,P_trig]=get_raster_pop_no_self_safety(parameter,spikes, DiscreteData,triglength_limit,ppms,lat_time_before_trig,lat_time_after_trig,exclude,lat_safety_margin,lat_binsize,isi,[]);

parameter='Touch';
exclude={'Puff';'Light';'Exclude';'Grooming'};
triglength_limit=[-inf inf];%ms
[RasterHTL,T_trig]=get_raster_pop_no_self_safety(parameter,spikes, DiscreteData,triglength_limit,ppms,lat_time_before_trig,lat_time_after_trig,exclude,lat_safety_margin,lat_binsize,isi,[]);
first_spike_times=cell(numel(RasterHPL),2);
fsl_p=ones(numel(RasterHPL),1);
for n=1:numel(RasterHPL)
    [a,ia,ib]=unique(RasterHPL{n}(:,2));
    first_spike_times{n,1}=nan(size(P_trig{n},1),1);
    first_spike_times{n,1}(a)=RasterHPL{n}(ia,1);
    [a,ia,ib]=unique(RasterHTL{n}(:,2));
    first_spike_times{n,2}=nan(size(T_trig{n},1),1);
    first_spike_times{n,2}(a)=RasterHTL{n}(ia,1);
if all(cellfun(@(x) any(~isnan(x)),first_spike_times(n,:)))
     [fsl_p(n)]=ranksum(first_spike_times{n,1},first_spike_times{n,2});
end 
end
%
FSLavg=cellfun(@(x) mean(x,'omitnan'),first_spike_times);
FSLstd=cellfun(@(x) std(x,0,'omitnan'),first_spike_times);
FSLmed=cellfun(@(x) median(x,'omitnan'),first_spike_times);
FSLiqr=cellfun(@(x) iqr(x),first_spike_times);
FSLany=cellfun(@(x) mean(~isnan(x)),first_spike_times);
FSLsig=true(size(FSLavg));%repmat(fsl_p<.05,1,2);%
figure('Position',[ 910   176   500   588])
[pI,pC,mS]=plot_ratecomp2(FSLavg,FSLsig,{'Puff','Touch'},{vpmr;pomr},{'VPM';'POm'},'mean','first spike latency (ms)');title('mean first spike latency')
if write_figs;exportgraphics (gcf,[figdir Fig_no 'firstspikelatency.pdf'],'BackgroundColor','none','ContentType','vector');end

SumStats(SumStatLine).VPM_mean_fsl_cond1=mS(1,1,1);
SumStats(SumStatLine).VPM_sem_fsl_cond1=mS(1,1,2);
SumStats(SumStatLine).VPM_mean_fsl_cond2=mS(1,2,1);
SumStats(SumStatLine).VPM_sem_fsl_cond2=mS(1,2,2);
SumStats(SumStatLine).VPM_p_fsl_cond12=pI(1);
SumStats(SumStatLine).POm_mean_fsl_cond1=mS(2,1,1);
SumStats(SumStatLine).POm_sem_fsl_cond1=mS(2,1,2);
SumStats(SumStatLine).POm_mean_fsl_cond2=mS(2,2,1);
SumStats(SumStatLine).POm_sem_fsl_cond2=mS(2,2,2);
SumStats(SumStatLine).POm_p_fsl_cond12=pI(2);
SumStats(SumStatLine).VPMPOm_p_fsl_cond1=pC(1);
SumStats(SumStatLine).VPMPOm_p_fsl_cond2=pC(2);


XCond=3;
SumStatTable{3,XCond}=SumStats(SumStatLine).VPM_mean_cond1_base;
SumStatTable{4,XCond}=SumStats(SumStatLine).VPM_sem_cond1_base;
SumStatTable{6,XCond}=SumStats(SumStatLine).VPM_mean_cond1_resp;
SumStatTable{7,XCond}=SumStats(SumStatLine).VPM_sem_cond1_resp;
SumStatTable{9,XCond}=SumStats(SumStatLine).VPM_p_cond1;
SumStatTable{3,XCond+1}=SumStats(SumStatLine).VPM_mean_cond2_base;
SumStatTable{4,XCond+1}=SumStats(SumStatLine).VPM_sem_cond2_base;
SumStatTable{6,XCond+1}=SumStats(SumStatLine).VPM_mean_cond2_resp;
SumStatTable{7,XCond+1}=SumStats(SumStatLine).VPM_sem_cond2_resp;
SumStatTable{9,XCond+1}=SumStats(SumStatLine).VPM_p_cond2;
SumStatTable{5,XCond}=SumStats(SumStatLine).VPM_p_base_cond12;
SumStatTable{8,XCond}=SumStats(SumStatLine).VPM_p_resp_cond12;
SumStatTable{16+3,XCond}=SumStats(SumStatLine).POm_mean_cond1_base;
SumStatTable{16+4,XCond}=SumStats(SumStatLine).POm_sem_cond1_base;
SumStatTable{16+6,XCond}=SumStats(SumStatLine).POm_mean_cond1_resp;
SumStatTable{16+7,XCond}=SumStats(SumStatLine).POm_sem_cond1_resp;
SumStatTable{16+9,XCond}=SumStats(SumStatLine).POm_p_cond1;
SumStatTable{16+3,XCond+1}=SumStats(SumStatLine).POm_mean_cond2_base;
SumStatTable{16+4,XCond+1}=SumStats(SumStatLine).POm_sem_cond2_base;
SumStatTable{16+6,XCond+1}=SumStats(SumStatLine).POm_mean_cond2_resp;
SumStatTable{16+7,XCond+1}=SumStats(SumStatLine).POm_sem_cond2_resp;
SumStatTable{16+9,XCond+1}=SumStats(SumStatLine).POm_p_cond2;
SumStatTable{16+5,XCond}=SumStats(SumStatLine).POm_p_base_cond12;
SumStatTable{16+8,XCond}=SumStats(SumStatLine).POm_p_resp_cond12;
SumStatTable{10,XCond}=SumStats(SumStatLine).VPM_mean_modulation_cond1;
SumStatTable{11,XCond}=SumStats(SumStatLine).VPM_sem_modulation_cond1;
SumStatTable{12,XCond}=SumStats(SumStatLine).VPM_p_modulation_cond10;
SumStatTable{10,XCond+1}=SumStats(SumStatLine).VPM_mean_modulation_cond2;
SumStatTable{11,XCond+1}=SumStats(SumStatLine).VPM_sem_modulation_cond2;
SumStatTable{12,XCond+1}=SumStats(SumStatLine).VPM_p_modulation_cond20;
SumStatTable{13,XCond}=SumStats(SumStatLine).VPM_p_modulation_cond12;
SumStatTable{16+10,XCond}=SumStats(SumStatLine).POm_mean_modulation_cond1;
SumStatTable{16+11,XCond}=SumStats(SumStatLine).POm_sem_modulation_cond1;
SumStatTable{16+12,XCond}=SumStats(SumStatLine).POm_p_modulation_cond10;
SumStatTable{16+10,XCond+1}=SumStats(SumStatLine).POm_mean_modulation_cond2;
SumStatTable{16+11,XCond+1}=SumStats(SumStatLine).POm_sem_modulation_cond2;
SumStatTable{16+12,XCond+1}=SumStats(SumStatLine).POm_p_modulation_cond20;
SumStatTable{16+13,XCond}=SumStats(SumStatLine).POm_p_modulation_cond12;
SumStatTable{36,XCond}=SumStats(SumStatLine).VPMPOm_p_modulation_cond1;
SumStatTable{36,XCond+1}=SumStats(SumStatLine).VPMPOm_p_modulation_cond2;
SumStatTable{34,XCond}=SumStats(SumStatLine).VPMPOm_p_cond1_base;
SumStatTable{34,XCond+1}=SumStats(SumStatLine).VPMPOm_p_cond2_base;
SumStatTable{35,XCond}=SumStats(SumStatLine).VPMPOm_p_cond1_resp;
SumStatTable{35,XCond+1}=SumStats(SumStatLine).VPMPOm_p_cond2_resp;
SumStatTable{14,XCond}=SumStats(SumStatLine).VPM_mean_fsl_cond1;
SumStatTable{15,XCond}=SumStats(SumStatLine).VPM_sem_fsl_cond1;
SumStatTable{14,XCond+1}=SumStats(SumStatLine).VPM_mean_fsl_cond2;
SumStatTable{15,XCond+1}=SumStats(SumStatLine).VPM_sem_fsl_cond2;
SumStatTable{16,XCond}=SumStats(SumStatLine).VPM_p_fsl_cond12;
SumStatTable{16+14,XCond}=SumStats(SumStatLine).POm_mean_fsl_cond1;
SumStatTable{16+15,XCond}=SumStats(SumStatLine).POm_sem_fsl_cond1;
SumStatTable{16+14,XCond+1}=SumStats(SumStatLine).POm_mean_fsl_cond2;
SumStatTable{16+15,XCond+1}=SumStats(SumStatLine).POm_sem_fsl_cond2;
SumStatTable{16+16,XCond}=SumStats(SumStatLine).POm_p_fsl_cond12;
SumStatTable{37,XCond}=SumStats(SumStatLine).VPMPOm_p_fsl_cond1;
SumStatTable{37,XCond+1}=SumStats(SumStatLine).VPMPOm_p_fsl_cond2;
% 
 figure('Position',[ 910   176   500   588])
[pIs,pCs,mSs]=plot_ratecomp2(FSLstd,FSLsig,{'Puff','Touch'},{vpmr;pomr},{'VPM';'POm'},'mean','std first spike latency (ms)');title('std first spike latency')
 if write_figs;exportgraphics (gcf,[figdir Fig_no 'stdfirstspikelatency.pdf'],'BackgroundColor','none','ContentType','vector');end
% 
% figure('Position',[ 910   176   500   588])
% plot_ratecomp2(FSLmed,FSLsig,{'Puff','Touch'},{vpmr;pomr},{'VPM';'POm'},'mean','first spike latency (ms)');title('median first spike latency')
% if write_figs;exportgraphics (gcf,[figdir Fig_no 'firstspikelatencymedian.pdf'],'BackgroundColor','none','ContentType','vector');end
% 
% figure('Position',[ 910   176   500   588])
% plot_ratecomp2(FSLiqr,FSLsig,{'Puff','Touch'},{vpmr;pomr},{'VPM';'POm'},'mean','std first spike latency (ms)');title('iqr first spike latency')
% if write_figs;exportgraphics (gcf,[figdir Fig_no 'iqrfirstspikelatency.pdf'],'BackgroundColor','none','ContentType','vector');end
% 
% figure('Position',[ 910   176   500   588])
% plot_ratecomp2(FSLany,FSLsig,{'Puff','Touch'},{vpmr;pomr},{'VPM';'POm'},'mean','response probability');title('spike in window')
% if write_figs;exportgraphics (gcf,[figdir Fig_no 'responseprobability.pdf'],'BackgroundColor','none','ContentType','vector');end
%

%% Figure 1
Fig_no='Fig1_';
cell_example=vpmr(1);
trigs_e=[trigsP(cell_example);trigsT(cell_example)];
trig_names={'Puff';'Touch'};
psth_e={P_PSTH_instrate(cell_example,:);T_PSTH_instrate(cell_example,:)};
trace_e=cat(3,cat(1,MeanP{cell_example},MeanWp{cell_example}),cat(1,MeanT{cell_example},MeanWt{cell_example}));
sem_trace_e=cat(3,cat(1,nan(size(MeanP{cell_example})),std(TrialWp{cell_example},0,1)./sqrt(size(TrialWp{cell_example},1))),cat(1,nan(size(MeanT{cell_example})),std(TrialWt{cell_example},0,1)./sqrt(size(TrialWt{cell_example},1))));
trace_name={'trig','whisker angle'};
fig1v=figure('Position',[10         245         574         412],'Color','White');Name='VPM example';
plot_single_cell_example(trigs_e,trig_names,psth_e,DiscreteData(cell_example).Spikes,trace_e,sem_trace_e,trace_name,ppms,psth_bins,Name);
exportgraphics(fig1v,[figdir Fig_no 'example_VPM_PuffTouch.pdf'],'BackgroundColor','none','ContentType','vector')

%
cell_example=pomr(10);
trigs_e=[trigsP(cell_example);trigsT(cell_example)];
trig_names={'Puff';'Touch'};
psth_e={P_PSTH_instrate(cell_example,:);T_PSTH_instrate(cell_example,:)};
trace_e=cat(3,cat(1,MeanP{cell_example},MeanWp{cell_example}),cat(1,MeanT{cell_example},MeanWt{cell_example}));
sem_trace_e=cat(3,cat(1,nan(size(MeanP{cell_example})),std(TrialWp{cell_example},0,1)./sqrt(size(TrialWp{cell_example},1))),cat(1,nan(size(MeanT{cell_example})),std(TrialWt{cell_example},0,1)./sqrt(size(TrialWt{cell_example},1))));
% sem_trace_e=cat(3,cat(1,nan(size(MeanP{cell_example})),std(TrialWp{cell_example},0,1)),cat(1,nan(size(MeanT{cell_example})),std(TrialWt{cell_example},0,1)));
trace_name={'trig','whisker angle'};
fig1p=figure('Position',[10         245         574         412],'Color','White');Name='POm example';
plot_single_cell_example(trigs_e,trig_names,psth_e,DiscreteData(cell_example).Spikes,trace_e,sem_trace_e,trace_name,ppms,psth_bins,Name);
exportgraphics(fig1p,[figdir Fig_no 'example_POm_PuffTouch.pdf'],'BackgroundColor','none','ContentType','vector')


fig1c=figure;
coords=RecDB{:,'CellAtlasCoords'};
scatter3(coords(pom,1),-coords(pom,2),-coords(pom,3),5,'red','filled');hold on;
scatter3(coords(pomr,1),-coords(pomr,2),-coords(pomr,3),25,'red','filled');hold on;
scatter3(coords(vpm,1),-coords(vpm,2),-coords(vpm,3),5,'blue','filled');hold on;
scatter3(coords(vpmr,1),-coords(vpmr,2),-coords(vpmr,3),25,'blue','filled');hold on;
sp=100*round(min(coords([vpm;pom],1:3).*[1 -1 -1])/100)+[-100 100 100];
quiver3(sp(1),sp(2),sp(3),200,0,0,'off','color','k');
quiver3(sp(1),sp(2),sp(3),0,0,200,'off','color','k');
quiver3(sp(1),sp(2),sp(3),0,200,0,'off','color','k');
axis equal
view(gca,23.5702,16.4257)
exportgraphics(fig1c,[figdir Fig_no 'cell_locations_all.pdf'],'BackgroundColor','none','ContentType','vector')
fig1vv=figure;
coords=RecDB{:,'CellAtlasCoords'};
%scatter3(coords(pom,1),-coords(pom,2),-coords(pom,3),5,'red','filled');hold on;
scatter3(coords(pomr,1),-coords(pomr,2),-coords(pomr,3),25,'red','filled');hold on;
%scatter3(coords(vpm,1),-coords(vpm,2),-coords(vpm,3),5,'blue','filled');hold on;
scatter3(coords(vpmr,1),-coords(vpmr,2),-coords(vpmr,3),25,'blue','filled');hold on;
sp=100*round(min(coords([vpm;pom],1:3).*[1 -1 -1])/100)+[-100 100 100];
quiver3(sp(1),sp(2),sp(3),200,0,0,'off','color','k');
quiver3(sp(1),sp(2),sp(3),0,0,200,'off','color','k');
quiver3(sp(1),sp(2),sp(3),0,200,0,'off','color','k');
axis equal
view(gca,23.5702,16.4257)
exportgraphics(fig1vv,[figdir Fig_no 'cell_locations_r.pdf'],'BackgroundColor','none','ContentType','vector')

% 
% raw_data_example=vpmr(5);
% tE=(1:DiscreteData(raw_data_example).LengthInd)/(ppms*1000);
% tW=downsample(tE,20);
% load(['F:\Whisker\AnalysisMatFiles\' RecDB.Properties.RowNames{raw_data_example} '\' RecDB.Properties.RowNames{raw_data_example} 'analysis.mat'],'filteredResponse')
% %
% figure
% tiledlayout(2,1)
% a(1)=nexttile(1);
% plot(tE,filteredResponse.data,'k');box off;hold on
% plot(DiscreteData(raw_data_example).PuffStart/(ppms*1000),max(filteredResponse.data),'b.');
% plot(DiscreteData(raw_data_example).TouchStart/(ppms*1000),max(filteredResponse.data),'r.');
% a(2)=nexttile(2);
% plot(tW,Angle{raw_data_example},'Color',[.8 .8 .8]);box off;hold on
% plot(DiscreteData(raw_data_example).PuffStart/(ppms*1000),max(Angle{raw_data_example}),'b.');
% plot(DiscreteData(raw_data_example).TouchStart/(ppms*1000),max(Angle{raw_data_example}),'r.');
% 
% linkaxes(a,'x')

%% Figure 3
Fig_no='Fig3_';
SumStatLine=2;
SumStats(SumStatLine).condition_1='Quiescience';
SumStats(SumStatLine).condition_2='Whisking';
SumStats(SumStatLine).VPM_cell_IDs=vpmra';
SumStats(SumStatLine).POm_cell_IDs=pomra';
[N_animals]=get_N_animals(SumStats(SumStatLine).VPM_cell_IDs,SumStats(SumStatLine).POm_cell_IDs,RecDB);
SumStats(SumStatLine).VPM_N=numel(vpmra);
SumStats(SumStatLine).VPM_N_animals=N_animals(3);
SumStats(SumStatLine).POm_N=numel(pomra);
SumStats(SumStatLine).POm_N_animals=N_animals(4);
SumStats(SumStatLine).Total_N_animals=N_animals(1);
C=8;
SumStatTable{1,C}=SumStats(SumStatLine).VPM_N;
SumStatTable{2,C}=SumStats(SumStatLine).VPM_N_animals;
SumStatTable{17,C}=SumStats(SumStatLine).POm_N;
SumStatTable{18,C}=SumStats(SumStatLine).POm_N_animals;
SumStatTable{33,C}=SumStats(SumStatLine).Total_N_animals;
shuffN=1e4;
outsafety=250;%ms
insafety=25;
isi=10;
include={'Whisking'};
exclude={'Puff';'Touch';'Light';'Exclude';'Grooming'};
[RatesW,RatesWp]=get_DD_Rates2(DiscreteData,exclude,include, outsafety,insafety,ppms,shuffN);
RatesWp(RatesWp>.5)=RatesWp(RatesWp>.5)-1;

outsafety=250;%ms
insafety=25;
isi=10;
include={'NoWhisking'};
exclude={'Puff';'Touch';'Light';'Exclude';'Grooming'};
[RatesQ,RatesQp]=get_DD_Rates2(DiscreteData,exclude,include, outsafety,insafety,ppms,shuffN);
RatesQp(RatesQp>.5)=RatesQp(RatesQp>.5)-1;

figure('Position',[ 910   176   500   588])
[pI,pC,mS]=plot_ratecomp2(cat(2,RatesQ,RatesW),abs(RatesWp)<.05,{'Q','W'},{vpmra;pomra},{'VPM';'POm'},'mean','Rate (Hz)');title('Whisking Rate, responsive cells')
if write_figs;exportgraphics (gcf,[figdir Fig_no 'WhiskingRatesGenrealRespondingCellsWTouch.pdf'], 'BackgroundColor','none','ContentType','vector');end
%

%%

SumStats(SumStatLine).VPM_mean_cond1_resp=mS(1,1,1);
SumStats(SumStatLine).VPM_sem_cond1_resp=mS(1,1,2);
SumStats(SumStatLine).VPM_mean_cond2_resp=mS(1,2,1);
SumStats(SumStatLine).VPM_sem_cond2_resp=mS(1,2,2);
SumStats(SumStatLine).VPM_p_resp_cond12=pI(1);
SumStats(SumStatLine).POm_mean_cond1_resp=mS(2,1,1);
SumStats(SumStatLine).POm_sem_cond1_resp=mS(2,1,2);
SumStats(SumStatLine).POm_mean_cond2_resp=mS(2,2,1);
SumStats(SumStatLine).POm_sem_cond2_resp=mS(2,2,2);
SumStats(SumStatLine).POm_p_resp_cond12=pI(2);
SumStats(SumStatLine).VPMPOm_p_cond1_resp=pC(1);
SumStats(SumStatLine).VPMPOm_p_cond2_resp=pC(2);


MD_spontWhisk=mod_depth(RatesQ,RatesW);
figure('Position',[ 910   176   500   588])
[pI,pC,mS,p0]=plot_ratecomp2(repmat(MD_spontWhisk,1,2),RatesWp<.05,{'spont Whisk','spont Whisk'},{vpmra;pomra},{'VPM';'POm'},'mean','mod depth');title('Whisking modulation, responsive cells')

 if write_figs;exportgraphics (gcf,[figdir Fig_no 'modulationWhiskingRespondingCells.pdf'], 'BackgroundColor','none','ContentType','vector');end


SumStats(SumStatLine).VPM_mean_modulation_cond2=mS(1,2,1);
SumStats(SumStatLine).VPM_sem_modulation_cond2=mS(1,2,2);
SumStats(SumStatLine).VPM_p_modulation_cond20=p0(1,2);
SumStats(SumStatLine).POm_mean_modulation_cond2=mS(2,2,1);
SumStats(SumStatLine).POm_sem_modulation_cond2=mS(2,2,2);
SumStats(SumStatLine).POm_p_modulation_cond20=p0(2,2);
SumStats(SumStatLine).VPMPOm_p_modulation_cond2=pC(2);

XCond=8;
SumStatTable{3,XCond}=SumStats(SumStatLine).VPM_mean_cond1_resp;
SumStatTable{4,XCond}=SumStats(SumStatLine).VPM_sem_cond1_resp;
SumStatTable{6,XCond}=SumStats(SumStatLine).VPM_mean_cond2_resp;
SumStatTable{7,XCond}=SumStats(SumStatLine).VPM_sem_cond2_resp;
 SumStatTable{9,XCond}=SumStats(SumStatLine).VPM_p_resp_cond12;
% SumStatTable{8,XCond}=SumStats(SumStatLine).VPM_p_resp_cond12;
SumStatTable{16+3,XCond}=SumStats(SumStatLine).POm_mean_cond1_resp;
SumStatTable{16+4,XCond}=SumStats(SumStatLine).POm_sem_cond1_resp;
SumStatTable{16+6,XCond}=SumStats(SumStatLine).POm_mean_cond2_resp;
SumStatTable{16+7,XCond}=SumStats(SumStatLine).POm_sem_cond2_resp;
SumStatTable{16+9,XCond}=SumStats(SumStatLine).POm_p_resp_cond12;
%SumStatTable{16+8,XCond}=SumStats(SumStatLine).POm_p_resp_cond12;
SumStatTable{10,XCond}=SumStats(SumStatLine).VPM_mean_modulation_cond2;
SumStatTable{11,XCond}=SumStats(SumStatLine).VPM_sem_modulation_cond2;
SumStatTable{12,XCond}=SumStats(SumStatLine).VPM_p_modulation_cond20;
% SumStatTable{13,XCond}=SumStats(SumStatLine).VPM_p_modulation_cond12;
SumStatTable{16+10,XCond}=SumStats(SumStatLine).POm_mean_modulation_cond2;
SumStatTable{16+11,XCond}=SumStats(SumStatLine).POm_sem_modulation_cond2;
SumStatTable{16+12,XCond}=SumStats(SumStatLine).POm_p_modulation_cond20;
% SumStatTable{16+13,XCond}=SumStats(SumStatLine).POm_p_modulation_cond12;
SumStatTable{36,XCond}=SumStats(SumStatLine).VPMPOm_p_modulation_cond2;
SumStatTable{34,XCond}=SumStats(SumStatLine).VPMPOm_p_cond1_resp;
SumStatTable{35,XCond}=SumStats(SumStatLine).VPMPOm_p_cond2_resp;



exclude={'Pole';'Light';'Exclude';'Grooming';'Puff'};
include=[];
safety_margin=[100 200];
min_t=500;%minimal time in bin

phasebinN=16;
ampbins=[0:5:30 inf];
abssetpbins=[0:5:30 inf];
absangbins=[0:5:40 inf];
setpbins=[-90 -25:5:30 inf];
angbins=[-90 -30:5:45 inf];
wbins={phasebinN;ampbins;abssetpbins;absangbins;setpbins;angbins};
[Num_bin_occ,rates_to_plot,plot_x,name_vars,PhaseMod,whisker_corrs,name_whisker_corrs,neuron_correlations_with_params]=get_whisking_in_air_modulation(wbins,min_t,DiscreteData,exclude,include,safety_margin,[],ppms,amp_thresh,Amp,Setp,Angle,Phase);


cell_select={vpmra;pomra};
cell_select_names={'VPM';'POm'};
cell_colors={'r';'b'};
figure('Position',[148         198        1593         688])
tt= tiledlayout(1,2);
title(tt,'mean whisking in air correlations amp')
v=1;
co=nan(size(rates_to_plot{v},1),2);
for n=1:size(rates_to_plot{v},1)
    tempresp=rates_to_plot{v}(n,:)';
    tempnan=~isnan(tempresp);
    tempnan(end)=false;
    if sum(tempnan)>2
        [co(n,1),co(n,2)]=corr(plot_x{v}(tempnan)',tempresp(tempnan));
    end
end
cosig=co(:,2)<sig;
for c=1:numel(cell_select)
    nexttile(c);

    plot(plot_x{v},rates_to_plot{v}(cell_select{c}(~cosig(cell_select{c})),:),'Color',[.75 .75 .75]);hold on
    plot(plot_x{v},rates_to_plot{v}(cell_select{c}(cosig(cell_select{c})),:),'Color',[.25 .25 .25]);hold on
    mA=mean(rates_to_plot{v}(cell_select{c},:),1,'omitnan');
    [cox, coxp]=corr(plot_x{v}',mA');
    plot(plot_x{v},mA,cell_colors{c});hold off;box off
    ylabel('mean Rate (Hz)');
    xlabel(name_vars{v});
    title(cell_select_names{c},sprintf('%s c=%.2f, p=%.3g',name_vars{v},cox,coxp ))
    xlim(plot_x{v}([1 end]));
end

if write_figs;exportgraphics (gcf,[figdir Fig_no 'response_relative_to_whisking_in_air_amp.pdf'], 'BackgroundColor','none','ContentType','vector');end


figure('Position',[213         667        1509         563],'Color','white')
tt=tiledlayout(4,1);
cell_select=vpmra(6:7);
title(tt,'phase examples')
plotx=linspace(-pi,pi,100);
for r=1:numel(cell_select)
nexttile
bar(plot_x{6},rates_to_plot{6}(cell_select(r),:),1,'FaceColor',[.7 .7 .7],'EdgeColor','k');hold on;
plot(plotx,PhaseMod.fitresult{cell_select(r)}(plotx),'Color','b','LineWidth',3)
legend off
xlim([-pi pi])
box off;
title('VPM')
end
cell_select=pomra([2 8]);
for r=1:numel(cell_select)
nexttile
bar(plot_x{6},rates_to_plot{6}(cell_select(r),:),1,'FaceColor',[.7 .7 .7],'EdgeColor','k');hold on;
plot(plotx,PhaseMod.fitresult{cell_select(r)}(plotx),'Color','m','LineWidth',3)
legend off
xlim([-pi pi])
box off;
title('POm')
end
if write_figs;exportgraphics (gcf,[figdir Fig_no 'phase_lock_example_vpmpom.pdf'], 'BackgroundColor','none','ContentType','vector');end

cell_select={vpmra;pomra};
figure('Position',[680   558   560   420],'Color','white')
tiledlayout(2,3);
nexttile(1,[2 2]);
polarscatter(PhaseMod.pref_phase(cell_select{1}(PhaseMod.p(cell_select{1})>.05)),PhaseMod.SNR(cell_select{1}(PhaseMod.p(cell_select{1})>.05)),75,'m.');hold on
polarscatter(PhaseMod.pref_phase(cell_select{1}(PhaseMod.p(cell_select{1})<=.05)),PhaseMod.SNR(cell_select{1}(PhaseMod.p(cell_select{1})<=.05)),75,'r.');hold on
polarscatter(PhaseMod.pref_phase(cell_select{2}(PhaseMod.p(cell_select{2})>.05)),PhaseMod.SNR(cell_select{2}(PhaseMod.p(cell_select{2})>.05)),75,'c.');hold on
polarscatter(PhaseMod.pref_phase(cell_select{2}(PhaseMod.p(cell_select{2})<=.05)),PhaseMod.SNR(cell_select{2}(PhaseMod.p(cell_select{2})<=.05)),75,'b.');hold off

title('phase SNR','preferred phase');
nexttile(3);
histogram(PhaseMod.SNR(cell_select{1}),0:.25:max(PhaseMod.SNR(cat(1,cell_select{1},cell_select{2})))+.2,'FaceColor','r');hold on;box off
title('phase SNR');
nexttile(6);
histogram(PhaseMod.SNR(cell_select{2}),0:.25:max(PhaseMod.SNR(cat(1,cell_select{1},cell_select{2})))+.2,'FaceColor','b');hold on;box off
if write_figs;exportgraphics (gcf,[figdir Fig_no 'phase_mod_whisking_in_air.pdf'], 'BackgroundColor','none','ContentType','vector');end


%%


SumStatLine=3;
SumStats(SumStatLine).condition_1='Puff nonWhisking';
SumStats(SumStatLine).condition_2='Puff Whisking';
SumStats(SumStatLine).VPM_cell_IDs=vpmr';
SumStats(SumStatLine).POm_cell_IDs=pomr';
[N_animals]=get_N_animals(SumStats(SumStatLine).VPM_cell_IDs,SumStats(SumStatLine).POm_cell_IDs,RecDB);
SumStats(SumStatLine).VPM_N=numel(vpmr);
SumStats(SumStatLine).VPM_N_animals=N_animals(3);
SumStats(SumStatLine).POm_N=numel(pomr);
SumStats(SumStatLine).POm_N_animals=N_animals(4);
SumStats(SumStatLine).Total_N_animals=N_animals(1);

C=5;
SumStatTable{1,C}=SumStats(SumStatLine).VPM_N;
SumStatTable{2,C}=SumStats(SumStatLine).VPM_N_animals;
SumStatTable{17,C}=SumStats(SumStatLine).POm_N;
SumStatTable{18,C}=SumStats(SumStatLine).POm_N_animals;
SumStatTable{33,C}=SumStats(SumStatLine).Total_N_animals;

parameter='Puff';
exclude={'Pole';'Light';'Exclude';'Grooming'};
include=[];
puff_triglength_limit=[28 32];%ms
% non Whisking
contPlot={MeanWpnw,MeanPnW,MeanApnw};
trigLen=cellfun(@(x) ones(1,numel(x))*30*20,trigsPnW,'UniformOutput',0);
parameter=[trigsPnW trigLen];
figure('Position',[38 260 1413 525],'Color','White');Name='Puff noWh';
plot_Pop_PSTH_wrapper3(parameter, DiscreteData,puff_triglength_limit,ppms,psth_time_before_trig,psth_time_after_trig,exclude,psth_safety_margin,psth_binsize,[],Name,{vpmr;pomr},{'VPM';'POm'},include,'latency',-inf,contPlot,{'meanbefore','scaleminmax','mono'},'baseline corr. whisker angle');
if write_figs;exportgraphics(gcf,[figdir Fig_no 'Pop_Puff_nW_instRateResp_' datestr(now,'yymmdd') '.pdf'],'BackgroundColor','none','ContentType','vector');end

[H_PuffnW,~,p_PuffnW]=get_PSTH_pop(parameter, DiscreteData,puff_triglength_limit,ppms,rate_time_before_trig,rate_time_after_trig,exclude,psth_safety_margin,rate_binsize,isi,'left',[]);
H_PnW=cat(2,H_PuffnW{:})';
p_PuffnW=p_PuffnW{1}(:,2);
sig_PuffnW=p_PuffnW<sig;

% Whisking


contPlot={MeanWpw,MeanPW,MeanApw};
trigLen=cellfun(@(x) ones(1,numel(x))*30*20,trigsPW,'UniformOutput',0);
parameter=[trigsPW trigLen];
figure('Position',[38 260 1413 525],'Color','White');Name='Puff Wh';
plot_Pop_PSTH_wrapper3(parameter, DiscreteData,puff_triglength_limit,ppms,psth_time_before_trig,psth_time_after_trig,exclude,psth_safety_margin,psth_binsize,[],Name,{vpmr;pomr},{'VPM';'POm'},include,'latency',-inf,contPlot,{'meanbefore','scaleminmax','mono'},'baseline corr. whisker angle');
if write_figs;exportgraphics(gcf,[figdir Fig_no 'Pop_Puff_W_instRateResp_' datestr(now,'yymmdd') '.pdf'],'BackgroundColor','none');end

[H_PuffW,~,p_PuffW]=get_PSTH_pop(parameter, DiscreteData,puff_triglength_limit,ppms,rate_time_before_trig,rate_time_after_trig,exclude,psth_safety_margin,rate_binsize,isi,'left',[]);
H_PW=cat(2,H_PuffW{:})';
p_PuffW=p_PuffW{1}(:,2);
sig_PuffW=p_PuffW<sig;




figure('Position',[910   176   500   588])
[pI,pC,mS]=plot_ratecomp3(cat(3,H_PnW,H_PW),cat(2,sig_PuffnW,sig_PuffW),{'0';'1';'0';'1'},vpmr,{'QPuff','WPuff'},'mean','Rate [Hz]');
title('Q/W Puff Response VPM', 'resp cells w touch')
if write_figs;exportgraphics(gcf,[figdir Fig_no 'RatesQWPuff_VPM_RespondingCellsWithTouch.pdf'], 'BackgroundColor','none','ContentType','vector');end

SumStats(SumStatLine).VPM_mean_cond1_base=mS(1,1,1);
SumStats(SumStatLine).VPM_sem_cond1_base=mS(1,1,2);
SumStats(SumStatLine).VPM_mean_cond1_resp=mS(1,2,1);
SumStats(SumStatLine).VPM_sem_cond1_resp=mS(1,2,2);
SumStats(SumStatLine).VPM_p_cond1=pI(1);
SumStats(SumStatLine).VPM_mean_cond2_base=mS(2,1,1);
SumStats(SumStatLine).VPM_sem_cond2_base=mS(2,1,2);
SumStats(SumStatLine).VPM_mean_cond2_resp=mS(2,2,1);
SumStats(SumStatLine).VPM_sem_cond2_resp=mS(2,2,2);
SumStats(SumStatLine).VPM_p_cond2=pI(2);
SumStats(SumStatLine).VPM_p_base_cond12=pC(1);
SumStats(SumStatLine).VPM_p_resp_cond12=pC(2);




figure('Position',[910   176   500   588])
[pI,pC,mS]=plot_ratecomp3(cat(3,H_PnW,H_PW),cat(2,sig_PuffnW,sig_PuffW),{'0';'1';'0';'1'},pomr,{'QPuff','WPuff'},'mean','Rate [Hz]');
title('Q/W Puff Response POm', 'resp cells w touch')
if write_figs;exportgraphics(gcf,[figdir Fig_no 'RatesQWPuff_POm_RespondingCellsWithTouch.pdf'], 'BackgroundColor','none','ContentType','vector');end
SumStats(SumStatLine).POm_mean_cond1_base=mS(1,1,1);
SumStats(SumStatLine).POm_sem_cond1_base=mS(1,1,2);
SumStats(SumStatLine).POm_mean_cond1_resp=mS(1,2,1);
SumStats(SumStatLine).POm_sem_cond1_resp=mS(1,2,2);
SumStats(SumStatLine).POm_p_cond1=pI(1);
SumStats(SumStatLine).POm_mean_cond2_base=mS(2,1,1);
SumStats(SumStatLine).POm_sem_cond2_base=mS(2,1,2);
SumStats(SumStatLine).POm_mean_cond2_resp=mS(2,2,1);
SumStats(SumStatLine).POm_sem_cond2_resp=mS(2,2,2);
SumStats(SumStatLine).POm_p_cond2=pI(2);
SumStats(SumStatLine).POm_p_base_cond12=pC(1);
SumStats(SumStatLine).POm_p_resp_cond12=pC(2);


pRc=nan(2,2);
for c=1:2
pRc(1,c)=ranksum(H_PnW(vpmr,c,1),H_PnW(pomr,c,1));
end
for c=1:2
pRc(2,c)=ranksum(H_PW(vpmr,c,1),H_PW(pomr,c,1));
end
SumStats(SumStatLine).VPMPOm_p_cond1_base=pRc(1,1);
SumStats(SumStatLine).VPMPOm_p_cond2_base=pRc(2,1);
SumStats(SumStatLine).VPMPOm_p_cond1_resp=pRc(1,2);
SumStats(SumStatLine).VPMPOm_p_cond2_resp=pRc(2,2);

p_PuffQWr=sig_test_comp_trig({trigsPnW;trigsPW},DiscreteData,rate_bins,2,'both',ppms);
sig_PuffQWr=p_PuffQWr<sig;
p_PuffQWb=sig_test_comp_trig({trigsPnW;trigsPW},DiscreteData,rate_bins,1,'both',ppms);
sig_PuffQWb=p_PuffQWb<sig;
% p_PuffQWq=sig_test_comp_trig({trigsPnW;trigsPW},DiscreteData,rate_bins,'quot','both',ppms);
% sig_PuffQWq=p_PuffQWq<sig;

% 
% tempD=cat(3,cat(2,H_PnW(:,2)./H_PnW(:,1),H_PW(:,2)./H_PW(:,1)),log2(cat(2,H_PnW(:,2)./H_PnW(:,1),H_PW(:,2)./H_PW(:,1))));

tempD=log2(cat(2,H_PnW(:,2)./H_PnW(:,1),H_PW(:,2)./H_PW(:,1)));
figure('Position',[ 910   176   500   588])
plot_ratecomp2(tempD,sig_PuffQWr,{'QPuff','WPuff'},{vpmr;pomr},{'VPM';'POm'},'mean','log2 Rate fold change');title('Q/W Puff fold change Response, responsive cells')
if write_figs;exportgraphics (gcf,[figdir Fig_no 'foldchangeRatesQWRespondingCellsWT.pdf'], 'BackgroundColor','none','ContentType','vector');end

tempD=cat(2,mod_depth(H_PnW),mod_depth(H_PW));
figure('Position',[ 910   176   500   588])
[pI,pC,mS,p0]=plot_ratecomp2(tempD,sig_PuffQWr,{'QPuff','WPuff'},{vpmr;pomr},{'VPM';'POm'},'mean','mod depth');title('Q/W Puff modulation, responsive cells')

 if write_figs;exportgraphics (gcf,[figdir Fig_no 'modulationQWRespondingCellsLog.pdf'], 'BackgroundColor','none','ContentType','vector');end

tempD=cat(2,H_PnW(:,2)./H_PnW(:,1),H_PW(:,2)./H_PW(:,1));
figure('Position',[ 910   176   500   588])
plot_ratecomp2(tempD,sig_PuffQWr,{'QPuff','WPuff'},{vpmr;pomr},{'VPM';'POm'},'mean','log2 Rate fold change');title('Q/W Puff fold change Response, responsive cells')
 set(gca,'YScale','log')
 set(gca,'YTick',[0.5 1:10 12:2:20 30 40 50 60])
 ylim([0.5 80])
 if write_figs;exportgraphics (gcf,[figdir Fig_no 'foldchangeRatesQWRespondingCellsLog.pdf'], 'BackgroundColor','none','ContentType','vector');end


SumStats(SumStatLine).VPM_mean_modulation_cond1=mS(1,1,1);
SumStats(SumStatLine).VPM_sem_modulation_cond1=mS(1,1,2);
SumStats(SumStatLine).VPM_p_modulation_cond10=p0(1,1);
SumStats(SumStatLine).VPM_mean_modulation_cond2=mS(1,2,1);
SumStats(SumStatLine).VPM_sem_modulation_cond2=mS(1,2,2);
SumStats(SumStatLine).VPM_p_modulation_cond20=p0(1,2);
SumStats(SumStatLine).VPM_p_modulation_cond12=pI(1);

SumStats(SumStatLine).POm_mean_modulation_cond1=mS(2,1,1);
SumStats(SumStatLine).POm_sem_modulation_cond1=mS(2,1,2);
SumStats(SumStatLine).POm_p_modulation_cond10=p0(2,1);
SumStats(SumStatLine).POm_mean_modulation_cond2=mS(2,2,1);
SumStats(SumStatLine).POm_sem_modulation_cond2=mS(2,2,2);
SumStats(SumStatLine).POm_p_modulation_cond20=p0(2,2);
SumStats(SumStatLine).POm_p_modulation_cond12=pI(2);

SumStats(SumStatLine).VPMPOm_p_modulation_cond1=pC(1);
SumStats(SumStatLine).VPMPOm_p_modulation_cond2=pC(2);

%% Supp Figure 3
Fig_no='EDFig3_';
lat_binsize=1;
lat_time_before_trig=0;
lat_time_after_trig=75;
lat_safety_margin=lat_time_after_trig;


puff_triglength_limit=[28 32];%ms
% non Whisking
trigLen=cellfun(@(x) ones(1,numel(x))*30*20,trigsPnW,'UniformOutput',0);
parameter=[trigsPnW trigLen];
tempDD=struct2cell(DiscreteData);
spikes=squeeze(tempDD(1,1,:));clear temp*
[RasterHPL_nW,P_trig]=get_raster_pop_no_self_safety(parameter,spikes, DiscreteData,triglength_limit,ppms,lat_time_before_trig,lat_time_after_trig,exclude,lat_safety_margin,lat_binsize,isi,[]);
%whisking
trigLen=cellfun(@(x) ones(1,numel(x))*30*20,trigsPW,'UniformOutput',0);
parameter=[trigsPW trigLen];
[RasterHPL_W,T_trig]=get_raster_pop_no_self_safety(parameter,spikes, DiscreteData,triglength_limit,ppms,lat_time_before_trig,lat_time_after_trig,exclude,lat_safety_margin,lat_binsize,isi,[]);
first_spike_times=cell(numel(RasterHPL_nW),2);
fsl_p=ones(numel(RasterHPL_nW),1);
for n=1:numel(RasterHPL_nW)
    [a,ia,ib]=unique(RasterHPL_nW{n}(:,2));
    first_spike_times{n,1}=nan(size(P_trig{n},1),1);
    first_spike_times{n,1}(a)=RasterHPL_nW{n}(ia,1);
    [a,ia,ib]=unique(RasterHPL_W{n}(:,2));
    first_spike_times{n,2}=nan(size(T_trig{n},1),1);
    first_spike_times{n,2}(a)=RasterHPL_W{n}(ia,1);
if all(cellfun(@(x) any(~isnan(x)),first_spike_times(n,:)))
     [fsl_p(n)]=ranksum(first_spike_times{n,1},first_spike_times{n,2});
end
end
%
FSLavg=cellfun(@(x) mean(x,'omitnan'),first_spike_times);
FSLstd=cellfun(@(x) std(x,0,'omitnan'),first_spike_times);
FSLmed=cellfun(@(x) median(x,'omitnan'),first_spike_times);
FSLiqr=cellfun(@(x) iqr(x),first_spike_times);
FSLany=cellfun(@(x) mean(~isnan(x)),first_spike_times);
FSLsig=true(size(FSLavg));%repmat(fsl_p<.05,1,2);%
figure('Position',[ 910   176   500   588])
[pI,pC,mS]=plot_ratecomp2(FSLavg,FSLsig,{'QPuff','WPuff'},{vpmr;pomr},{'VPM';'POm'},'mean','first spike latency (ms)');title('mean first spike latency')
if write_figs;exportgraphics (gcf,[figdir Fig_no 'QWPuff_firstspikelatency.pdf'],'BackgroundColor','none','ContentType','vector');end

SumStats(SumStatLine).VPM_mean_fsl_cond1=mS(1,1,1);
SumStats(SumStatLine).VPM_sem_fsl_cond1=mS(1,1,2);
SumStats(SumStatLine).VPM_mean_fsl_cond2=mS(1,2,1);
SumStats(SumStatLine).VPM_sem_fsl_cond2=mS(1,2,2);
SumStats(SumStatLine).VPM_p_fsl_cond12=pI(1);

SumStats(SumStatLine).POm_mean_fsl_cond1=mS(2,1,1);
SumStats(SumStatLine).POm_sem_fsl_cond1=mS(2,1,2);
SumStats(SumStatLine).POm_mean_fsl_cond2=mS(2,2,1);
SumStats(SumStatLine).POm_sem_fsl_cond2=mS(2,2,2);
SumStats(SumStatLine).POm_p_fsl_cond12=pI(2);

SumStats(SumStatLine).VPMPOm_p_fsl_cond1=pC(1);
SumStats(SumStatLine).VPMPOm_p_fsl_cond2=pC(2);

XCond=5;
SumStatTable{3,XCond}=SumStats(SumStatLine).VPM_mean_cond1_base;
SumStatTable{4,XCond}=SumStats(SumStatLine).VPM_sem_cond1_base;
SumStatTable{6,XCond}=SumStats(SumStatLine).VPM_mean_cond1_resp;
SumStatTable{7,XCond}=SumStats(SumStatLine).VPM_sem_cond1_resp;
SumStatTable{9,XCond}=SumStats(SumStatLine).VPM_p_cond1;
SumStatTable{3,XCond+1}=SumStats(SumStatLine).VPM_mean_cond2_base;
SumStatTable{4,XCond+1}=SumStats(SumStatLine).VPM_sem_cond2_base;
SumStatTable{6,XCond+1}=SumStats(SumStatLine).VPM_mean_cond2_resp;
SumStatTable{7,XCond+1}=SumStats(SumStatLine).VPM_sem_cond2_resp;
SumStatTable{9,XCond+1}=SumStats(SumStatLine).VPM_p_cond2;
SumStatTable{5,XCond}=SumStats(SumStatLine).VPM_p_base_cond12;
SumStatTable{8,XCond}=SumStats(SumStatLine).VPM_p_resp_cond12;
SumStatTable{16+3,XCond}=SumStats(SumStatLine).POm_mean_cond1_base;
SumStatTable{16+4,XCond}=SumStats(SumStatLine).POm_sem_cond1_base;
SumStatTable{16+6,XCond}=SumStats(SumStatLine).POm_mean_cond1_resp;
SumStatTable{16+7,XCond}=SumStats(SumStatLine).POm_sem_cond1_resp;
SumStatTable{16+9,XCond}=SumStats(SumStatLine).POm_p_cond1;
SumStatTable{16+3,XCond+1}=SumStats(SumStatLine).POm_mean_cond2_base;
SumStatTable{16+4,XCond+1}=SumStats(SumStatLine).POm_sem_cond2_base;
SumStatTable{16+6,XCond+1}=SumStats(SumStatLine).POm_mean_cond2_resp;
SumStatTable{16+7,XCond+1}=SumStats(SumStatLine).POm_sem_cond2_resp;
SumStatTable{16+9,XCond+1}=SumStats(SumStatLine).POm_p_cond2;
SumStatTable{16+5,XCond}=SumStats(SumStatLine).POm_p_base_cond12;
SumStatTable{16+8,XCond}=SumStats(SumStatLine).POm_p_resp_cond12;
SumStatTable{10,XCond}=SumStats(SumStatLine).VPM_mean_modulation_cond1;
SumStatTable{11,XCond}=SumStats(SumStatLine).VPM_sem_modulation_cond1;
SumStatTable{12,XCond}=SumStats(SumStatLine).VPM_p_modulation_cond10;
SumStatTable{10,XCond+1}=SumStats(SumStatLine).VPM_mean_modulation_cond2;
SumStatTable{11,XCond+1}=SumStats(SumStatLine).VPM_sem_modulation_cond2;
SumStatTable{12,XCond+1}=SumStats(SumStatLine).VPM_p_modulation_cond20;
SumStatTable{13,XCond}=SumStats(SumStatLine).VPM_p_modulation_cond12;
SumStatTable{16+10,XCond}=SumStats(SumStatLine).POm_mean_modulation_cond1;
SumStatTable{16+11,XCond}=SumStats(SumStatLine).POm_sem_modulation_cond1;
SumStatTable{16+12,XCond}=SumStats(SumStatLine).POm_p_modulation_cond10;
SumStatTable{16+10,XCond+1}=SumStats(SumStatLine).POm_mean_modulation_cond2;
SumStatTable{16+11,XCond+1}=SumStats(SumStatLine).POm_sem_modulation_cond2;
SumStatTable{16+12,XCond+1}=SumStats(SumStatLine).POm_p_modulation_cond20;
SumStatTable{16+13,XCond}=SumStats(SumStatLine).POm_p_modulation_cond12;
SumStatTable{36,XCond}=SumStats(SumStatLine).VPMPOm_p_modulation_cond1;
SumStatTable{36,XCond+1}=SumStats(SumStatLine).VPMPOm_p_modulation_cond2;
SumStatTable{34,XCond}=SumStats(SumStatLine).VPMPOm_p_cond1_base;
SumStatTable{34,XCond+1}=SumStats(SumStatLine).VPMPOm_p_cond2_base;
SumStatTable{35,XCond}=SumStats(SumStatLine).VPMPOm_p_cond1_resp;
SumStatTable{35,XCond+1}=SumStats(SumStatLine).VPMPOm_p_cond2_resp;
SumStatTable{14,XCond}=SumStats(SumStatLine).VPM_mean_fsl_cond1;
SumStatTable{15,XCond}=SumStats(SumStatLine).VPM_sem_fsl_cond1;
SumStatTable{14,XCond+1}=SumStats(SumStatLine).VPM_mean_fsl_cond2;
SumStatTable{15,XCond+1}=SumStats(SumStatLine).VPM_sem_fsl_cond2;
SumStatTable{16,XCond}=SumStats(SumStatLine).VPM_p_fsl_cond12;
SumStatTable{16+14,XCond}=SumStats(SumStatLine).POm_mean_fsl_cond1;
SumStatTable{16+15,XCond}=SumStats(SumStatLine).POm_sem_fsl_cond1;
SumStatTable{16+14,XCond+1}=SumStats(SumStatLine).POm_mean_fsl_cond2;
SumStatTable{16+15,XCond+1}=SumStats(SumStatLine).POm_sem_fsl_cond2;
SumStatTable{16+16,XCond}=SumStats(SumStatLine).POm_p_fsl_cond12;
SumStatTable{37,XCond}=SumStats(SumStatLine).VPMPOm_p_fsl_cond1;
SumStatTable{37,XCond+1}=SumStats(SumStatLine).VPMPOm_p_fsl_cond2;
% 
% figure('Position',[ 910   176   500   588])
% plot_ratecomp2(FSLstd,FSLsig,{'Puff','Touch'},{vpmr;pomr},{'VPM';'POm'},'mean','std first spike latency (ms)');title('std first spike latency')
% if write_figs;exportgraphics (gcf,[figdir Fig_no 'stdfirstspikelatency.pdf'],'BackgroundColor','none','ContentType','vector');end
% 
% figure('Position',[ 910   176   500   588])
% plot_ratecomp2(FSLmed,FSLsig,{'Puff','Touch'},{vpmr;pomr},{'VPM';'POm'},'mean','first spike latency (ms)');title('median first spike latency')
% if write_figs;exportgraphics (gcf,[figdir Fig_no 'firstspikelatencymedian.pdf'],'BackgroundColor','none','ContentType','vector');end
% 
% figure('Position',[ 910   176   500   588])
% plot_ratecomp2(FSLiqr,FSLsig,{'Puff','Touch'},{vpmr;pomr},{'VPM';'POm'},'mean','std first spike latency (ms)');title('iqr first spike latency')
% if write_figs;exportgraphics (gcf,[figdir Fig_no 'iqrfirstspikelatency.pdf'],'BackgroundColor','none','ContentType','vector');end
% 
% figure('Position',[ 910   176   500   588])
% plot_ratecomp2(FSLany,FSLsig,{'Puff','Touch'},{vpmr;pomr},{'VPM';'POm'},'mean','response probability');title('spike in window')
% if write_figs;exportgraphics (gcf,[figdir Fig_no 'responseprobability.pdf'],'BackgroundColor','none','ContentType','vector');end


%% Figure 4
split_time_before_trig=50;%ms
split_time_after_trig=50;%ms
split_safety_margin=150;%time_after_trig;
Wdname='C:\Data\Whisker\Wdec\Wdec_by_field.mat';
split_binsize=1;%ms
smoothwidth=milliseconds(10);%ms
meth='gaussian';
binsR=[-50 0 50];
selc={vpmr;pomr};
splits=[0 100/3 200/3 100];
splitname={'Tertile'};
split_var_nameC={'Acceleration';'Curvature';'Interval';'Velocity';'Setp';'Amp';'Angle';'Phase'};
%split_var_nameC={'Curvature';'Velocity';'Acceleration';'Interval';'Amp'};
type_nameC={'rawbefore','rawafter','absbefore','absafter','rel','absrel'};
split_allR=true;%[true;false];%

[RasterHP,RasterHT,~,~,splitPlotValues]=split_kinematics(split_time_before_trig,split_time_after_trig,split_safety_margin,Wdname,split_binsize,smoothwidth,meth,binsR,...
    selc,splits,splitname,DiscreteData,ppms,split_var_nameC,type_nameC,split_allR,figdir,TrialWp,TrialWt,Amp,Phase,Angle,Setp);
splitPlotValuesNames={'kinematic param','normalization', 'tertile cutoffs' ,...
    'vpm puff means', 'vpm puff sem','vpm p(puff=0)' ,'vpm touch means', 'vpm touch sem','vpm p(touch=0)','vpm p puff vs touch','vpm puff/touch p 1 vs 3',...
'pom puff means', 'pom puff sem','pom p(puff=0)','pom touch means', 'pom touch sem','pom p(touch=0)','pom p puff vs touch','pom puff/touch p 1 vs 3'};
%%
SplitTable=cell2table(splitPlotValues,'VariableNames',splitPlotValuesNames);

writetable(SplitTable,'C:\Data\Whisker\splitstats.xlsx')
%% Figure 5
Fig_no='Fig5_';
vpmL=vpm(hasLight(vpm));
pomL=pom(hasLight(pom));
vpmrLa=vpmra(hasLight(vpmra));
pomrLa=pomra(hasLight(pomra));
SumStatLine=4;
SumStats(SumStatLine).condition_1='whisking no Light';
SumStats(SumStatLine).condition_2='whisking Light';
SumStats(SumStatLine).VPM_cell_IDs=vpmrLa';
SumStats(SumStatLine).POm_cell_IDs=pomrLa';
[N_animals]=get_N_animals(SumStats(SumStatLine).VPM_cell_IDs,SumStats(SumStatLine).POm_cell_IDs,RecDB);
SumStats(SumStatLine).VPM_N=numel(vpmrLa);
SumStats(SumStatLine).VPM_N_animals=N_animals(3);
SumStats(SumStatLine).POm_N=numel(pomrLa);
SumStats(SumStatLine).POm_N_animals=N_animals(4);
SumStats(SumStatLine).Total_N_animals=N_animals(1);
C=11;
SumStatTable{1,C}=SumStats(SumStatLine).VPM_N;
SumStatTable{2,C}=SumStats(SumStatLine).VPM_N_animals;
SumStatTable{17,C}=SumStats(SumStatLine).POm_N;
SumStatTable{18,C}=SumStats(SumStatLine).POm_N_animals;
SumStatTable{33,C}=SumStats(SumStatLine).Total_N_animals;

outsafety=200;%ms
insafety=25;
isi=10;
include={'Light'};
exclude={'NoWhisking';'Puff';'Touch';'Exclude';'Grooming'};
[RatesWL,RatesWLp]=get_DD_Rates2(DiscreteData,exclude,include, outsafety,insafety,ppms,shuffN);
RatesWLp(RatesWLp>.5)=RatesWLp(RatesWLp>.5)-1;
outsafety=200;%ms
insafety=25;
isi=10;
include={'Light'};
exclude={'Whisking';'Puff';'Touch';'Exclude';'Grooming'};
[RatesQL,RatesQLp]=get_DD_Rates2(DiscreteData,exclude,include, outsafety,insafety,ppms,shuffN);
RatesQLp(RatesQLp>.5)=RatesQLp(RatesQLp>.5)-1;

figure('Position',[ 910   176   500   588])
[pI,pC,mS]=plot_ratecomp3(cat(3,cat(2,RatesQ,RatesW),cat(2,RatesQL,RatesWL)),cat(2,abs(RatesWp)<.05,abs(RatesWLp)<.05),{'Q';'W';'QL';'WL'},vpmrLa,{'nL';'L'},'mean','Rate (Hz)');
title('VPM Whisking Rate nL/L, responsive cells w Light')
if write_figs;exportgraphics (gcf,[figdir Fig_no 'WhiskingRatesLightRespondingCellsWLight_vpm.pdf'], 'BackgroundColor','none','ContentType','vector');end

SumStats(SumStatLine).VPM_mean_cond1_base=mS(1,1,1);
SumStats(SumStatLine).VPM_sem_cond1_base=mS(1,1,2);
SumStats(SumStatLine).VPM_mean_cond1_resp=mS(1,2,1);
SumStats(SumStatLine).VPM_sem_cond1_resp=mS(1,2,2);
SumStats(SumStatLine).VPM_p_cond1=pI(1);
SumStats(SumStatLine).VPM_mean_cond2_base=mS(2,1,1);
SumStats(SumStatLine).VPM_sem_cond2_base=mS(2,1,2);
SumStats(SumStatLine).VPM_mean_cond2_resp=mS(2,2,1);
SumStats(SumStatLine).VPM_sem_cond2_resp=mS(2,2,2);
SumStats(SumStatLine).VPM_p_cond2=pI(2);
SumStats(SumStatLine).VPM_p_base_cond12=pC(1);
SumStats(SumStatLine).VPM_p_resp_cond12=pC(2);
%%
figure('Position',[ 910   176   500   588])
[pI,pC,mS]=plot_ratecomp3(cat(3,cat(2,RatesQ,RatesW),cat(2,RatesQL,RatesWL)),cat(2,abs(RatesWp)<.05,abs(RatesWLp)<.05),{'Q';'W';'QL';'WL'},pomrLa,{'nL';'L'},'mean','Rate (Hz)');
title('POm Whisking Rate nL/L, responsive cells w Light')
if write_figs;exportgraphics (gcf,[figdir Fig_no 'WhiskingRatesLightRespondingCellsWLight_pom.pdf'], 'BackgroundColor','none','ContentType','vector');end

SumStats(SumStatLine).POm_mean_cond1_base=mS(1,1,1);
SumStats(SumStatLine).POm_sem_cond1_base=mS(1,1,2);
SumStats(SumStatLine).POm_mean_cond1_resp=mS(1,2,1);
SumStats(SumStatLine).POm_sem_cond1_resp=mS(1,2,2);
SumStats(SumStatLine).POm_p_cond1=pI(1);
SumStats(SumStatLine).POm_mean_cond2_base=mS(2,1,1);
SumStats(SumStatLine).POm_sem_cond2_base=mS(2,1,2);
SumStats(SumStatLine).POm_mean_cond2_resp=mS(2,2,1);
SumStats(SumStatLine).POm_sem_cond2_resp=mS(2,2,2);
SumStats(SumStatLine).POm_p_cond2=pI(2);
SumStats(SumStatLine).POm_p_base_cond12=pC(1);
SumStats(SumStatLine).POm_p_resp_cond12=pC(2);

pRcL=nan(2,2);
HR_nL=cat(2,RatesQ,RatesW);
for c=1:2
pRcL(1,c)=ranksum(HR_nL(vpmrLa,c),HR_nL(pomrLa,c));
end
HR_L=cat(2,RatesQL,RatesWL);
for c=1:2
pRcL(2,c)=ranksum(HR_L(vpmrLa,c),HR_L(pomrLa,c));
end

SumStats(SumStatLine).VPMPOm_p_cond1_base=pRcL(1,1);
SumStats(SumStatLine).VPMPOm_p_cond2_base=pRcL(2,1);
SumStats(SumStatLine).VPMPOm_p_cond1_resp=pRcL(1,2);
SumStats(SumStatLine).VPMPOm_p_cond2_resp=pRcL(2,2);


MD_spontWhiskL=cat(2,mod_depth(RatesQ,RatesQL),mod_depth(RatesW,RatesWL));
figure('Position',[ 910   176   500   588])
[pI,pC,mS,p0]=plot_ratecomp2(MD_spontWhiskL,RatesWLp<.05,{'Q nL vs L','W nL vs L'},{vpmrLa;pomrLa},{'VPM';'POm'},'mean','mod depth');title('Whisking light modulation, responsive cells')

 if write_figs;exportgraphics (gcf,[figdir Fig_no 'modulationLightWhiskingRespondingCells.pdf'], 'BackgroundColor','none','ContentType','vector');end

MD_spontWhiskL2=cat(2,mod_depth(RatesQ,RatesW),mod_depth(RatesQL,RatesWL));
figure('Position',[ 910   176   500   588])
[pI,pC,mS,p0]=plot_ratecomp2(MD_spontWhiskL2,RatesWLp<.05,{'nL: Q vs W','L: Q vs W'},{vpmrLa;pomrLa},{'VPM';'POm'},'mean','mod depth');title('Whisking light modulation, responsive cells')

 if write_figs;exportgraphics (gcf,[figdir Fig_no 'modulationWhisking_respectiveLightRespondingCells.pdf'], 'BackgroundColor','none','ContentType','vector');end




SumStats(SumStatLine).VPM_mean_modulation_cond1=mS(1,1,1);
SumStats(SumStatLine).VPM_sem_modulation_cond1=mS(1,1,2);
SumStats(SumStatLine).VPM_p_modulation_cond10=p0(1,1);
SumStats(SumStatLine).VPM_mean_modulation_cond2=mS(1,2,1);
SumStats(SumStatLine).VPM_sem_modulation_cond2=mS(1,2,2);
SumStats(SumStatLine).VPM_p_modulation_cond20=p0(1,2);
SumStats(SumStatLine).VPM_p_modulation_cond12=pI(1);

SumStats(SumStatLine).POm_mean_modulation_cond1=mS(2,1,1);
SumStats(SumStatLine).POm_sem_modulation_cond1=mS(2,1,2);
SumStats(SumStatLine).POm_p_modulation_cond10=p0(2,1);
SumStats(SumStatLine).POm_mean_modulation_cond2=mS(2,2,1);
SumStats(SumStatLine).POm_sem_modulation_cond2=mS(2,2,2);
SumStats(SumStatLine).POm_p_modulation_cond20=p0(2,2);
SumStats(SumStatLine).POm_p_modulation_cond12=pI(2);

SumStats(SumStatLine).VPMPOm_p_modulation_cond1=pC(1);
SumStats(SumStatLine).VPMPOm_p_modulation_cond2=pC(2);

XCond=11;
SumStatTable{3,XCond}=SumStats(SumStatLine).VPM_mean_cond1_base;
SumStatTable{4,XCond}=SumStats(SumStatLine).VPM_sem_cond1_base;
SumStatTable{6,XCond}=SumStats(SumStatLine).VPM_mean_cond1_resp;
SumStatTable{7,XCond}=SumStats(SumStatLine).VPM_sem_cond1_resp;
SumStatTable{9,XCond}=SumStats(SumStatLine).VPM_p_cond1;
SumStatTable{3,XCond+1}=SumStats(SumStatLine).VPM_mean_cond2_base;
SumStatTable{4,XCond+1}=SumStats(SumStatLine).VPM_sem_cond2_base;
SumStatTable{6,XCond+1}=SumStats(SumStatLine).VPM_mean_cond2_resp;
SumStatTable{7,XCond+1}=SumStats(SumStatLine).VPM_sem_cond2_resp;
SumStatTable{9,XCond+1}=SumStats(SumStatLine).VPM_p_cond2;
SumStatTable{5,XCond}=SumStats(SumStatLine).VPM_p_base_cond12;
SumStatTable{8,XCond}=SumStats(SumStatLine).VPM_p_resp_cond12;
SumStatTable{16+3,XCond}=SumStats(SumStatLine).POm_mean_cond1_base;
SumStatTable{16+4,XCond}=SumStats(SumStatLine).POm_sem_cond1_base;
SumStatTable{16+6,XCond}=SumStats(SumStatLine).POm_mean_cond1_resp;
SumStatTable{16+7,XCond}=SumStats(SumStatLine).POm_sem_cond1_resp;
SumStatTable{16+9,XCond}=SumStats(SumStatLine).POm_p_cond1;
SumStatTable{16+3,XCond+1}=SumStats(SumStatLine).POm_mean_cond2_base;
SumStatTable{16+4,XCond+1}=SumStats(SumStatLine).POm_sem_cond2_base;
SumStatTable{16+6,XCond+1}=SumStats(SumStatLine).POm_mean_cond2_resp;
SumStatTable{16+7,XCond+1}=SumStats(SumStatLine).POm_sem_cond2_resp;
SumStatTable{16+9,XCond+1}=SumStats(SumStatLine).POm_p_cond2;
SumStatTable{16+5,XCond}=SumStats(SumStatLine).POm_p_base_cond12;
SumStatTable{16+8,XCond}=SumStats(SumStatLine).POm_p_resp_cond12;
SumStatTable{10,XCond}=SumStats(SumStatLine).VPM_mean_modulation_cond1;
SumStatTable{11,XCond}=SumStats(SumStatLine).VPM_sem_modulation_cond1;
SumStatTable{12,XCond}=SumStats(SumStatLine).VPM_p_modulation_cond10;
SumStatTable{10,XCond+1}=SumStats(SumStatLine).VPM_mean_modulation_cond2;
SumStatTable{11,XCond+1}=SumStats(SumStatLine).VPM_sem_modulation_cond2;
SumStatTable{12,XCond+1}=SumStats(SumStatLine).VPM_p_modulation_cond20;
SumStatTable{13,XCond}=SumStats(SumStatLine).VPM_p_modulation_cond12;
SumStatTable{16+10,XCond}=SumStats(SumStatLine).POm_mean_modulation_cond1;
SumStatTable{16+11,XCond}=SumStats(SumStatLine).POm_sem_modulation_cond1;
SumStatTable{16+12,XCond}=SumStats(SumStatLine).POm_p_modulation_cond10;
SumStatTable{16+10,XCond+1}=SumStats(SumStatLine).POm_mean_modulation_cond2;
SumStatTable{16+11,XCond+1}=SumStats(SumStatLine).POm_sem_modulation_cond2;
SumStatTable{16+12,XCond+1}=SumStats(SumStatLine).POm_p_modulation_cond20;
SumStatTable{16+13,XCond}=SumStats(SumStatLine).POm_p_modulation_cond12;
SumStatTable{36,XCond}=SumStats(SumStatLine).VPMPOm_p_modulation_cond1;
SumStatTable{36,XCond+1}=SumStats(SumStatLine).VPMPOm_p_modulation_cond2;
SumStatTable{34,XCond}=SumStats(SumStatLine).VPMPOm_p_cond1_base;
SumStatTable{34,XCond+1}=SumStats(SumStatLine).VPMPOm_p_cond2_base;
SumStatTable{35,XCond}=SumStats(SumStatLine).VPMPOm_p_cond1_resp;
SumStatTable{35,XCond+1}=SumStats(SumStatLine).VPMPOm_p_cond2_resp;


exclude={'Pole';'Exclude';'Grooming';'Puff'};
include={'Light'};
safety_marginEx=[5 75];
safety_marginIn=[10 0];
min_t=500;%minimal time in bin

phasebinN=16;
ampbins=0:3:65;
abssetpbins=0:3:40;
absangbins=0:3:65;
setpbins=-30:3:40;
angbins=-40:3:65;
wbins={phasebinN;ampbins;abssetpbins;absangbins;setpbins;angbins};
[Num_bin_occL,rates_to_plotL,plot_x,name_vars,PhaseModL,whisker_corrsL,name_whisker_corrs,neuron_correlations_with_paramsL]=get_whisking_in_air_modulation(wbins,min_t,DiscreteData,exclude,include,safety_marginEx,safety_marginIn,ppms,amp_thresh,Amp,Setp,Angle,Phase);

pref_phase_nLL=cat(2,PhaseMod.pref_phase',PhaseModL.pref_phase');
SNR_nLL=cat(2,PhaseMod.SNR',PhaseModL.SNR');
p_nLL=cat(2,PhaseMod.p',PhaseModL.p');
cell_select_nLL={vpmrLa,pomrLa};
s=200;
figure('Position',[ 222         460        1018         518],'Color','white')
tiledlayout (1,2)
nexttile;
Lx=1;Cx=1;
polarscatter(pref_phase_nLL(cell_select_nLL{Cx}(p_nLL(cell_select_nLL{Cx},Lx)>.05),Lx),...
                    SNR_nLL(cell_select_nLL{Cx}(p_nLL(cell_select_nLL{Cx},Lx)>.05),Lx),s,'m.');hold on
polarscatter(pref_phase_nLL(cell_select_nLL{Cx}(p_nLL(cell_select_nLL{Cx},Lx)<=.05),Lx),...
                    SNR_nLL(cell_select_nLL{Cx}(p_nLL(cell_select_nLL{Cx},Lx)<=.05),Lx),s,'r.');hold on
Lx=2;Cx=1;
polarscatter(pref_phase_nLL(cell_select_nLL{Cx}(p_nLL(cell_select_nLL{Cx},Lx)>.05),Lx),...
                    SNR_nLL(cell_select_nLL{Cx}(p_nLL(cell_select_nLL{Cx},Lx)>.05),Lx),25,'cd','filled');hold on
polarscatter(pref_phase_nLL(cell_select_nLL{Cx}(p_nLL(cell_select_nLL{Cx},Lx)<=.05),Lx),...
                    SNR_nLL(cell_select_nLL{Cx}(p_nLL(cell_select_nLL{Cx},Lx)<=.05),Lx),25,'bd','filled');hold on
polarplot(pref_phase_nLL(cell_select_nLL{Cx},:)',SNR_nLL(cell_select_nLL{Cx},:)' ,'k');hold on
title('VPM pref phase with light')
set(gca,'RLim', [0 2.5000])
nexttile;
Lx=1;Cx=2;
polarscatter(pref_phase_nLL(cell_select_nLL{Cx}(p_nLL(cell_select_nLL{Cx},Lx)>.05),Lx),...
                    SNR_nLL(cell_select_nLL{Cx}(p_nLL(cell_select_nLL{Cx},Lx)>.05),Lx),s,'m.');hold on
polarscatter(pref_phase_nLL(cell_select_nLL{Cx}(p_nLL(cell_select_nLL{Cx},Lx)<=.05),Lx),...
                    SNR_nLL(cell_select_nLL{Cx}(p_nLL(cell_select_nLL{Cx},Lx)<=.05),Lx),s,'r.');hold on
Lx=2;Cx=2;
polarscatter(pref_phase_nLL(cell_select_nLL{Cx}(p_nLL(cell_select_nLL{Cx},Lx)>.05),Lx),...
                    SNR_nLL(cell_select_nLL{Cx}(p_nLL(cell_select_nLL{Cx},Lx)>.05),Lx),25,'cd','filled');hold on
polarscatter(pref_phase_nLL(cell_select_nLL{Cx}(p_nLL(cell_select_nLL{Cx},Lx)<=.05),Lx),...
                    SNR_nLL(cell_select_nLL{Cx}(p_nLL(cell_select_nLL{Cx},Lx)<=.05),Lx),25,'bd','filled');hold on
polarplot(pref_phase_nLL(cell_select_nLL{Cx},:)',SNR_nLL(cell_select_nLL{Cx},:)','k');hold on
title('POm pref phase with light')
set(gca,'RLim', [0 2.5000])
if write_figs;exportgraphics (gcf,[figdir Fig_no 'phase_mod_change_whisking_in_air_light.pdf'], 'BackgroundColor','none','ContentType','vector');end

figure('Position',[ 910   176   500   588])
plot_ratecomp2(SNR_nLL,p_nLL(:,2)<=.05,{'no Light','Light'},{vpmrLa;pomrLa},{'VPM';'POm'},'mean','SNR');title('nL/L phase SNR')
if write_figs;exportgraphics (gcf,[figdir Fig_no 'SNRchange_nLL.pdf'], 'BackgroundColor','none','ContentType','vector');end

round(100.*[mean(p_nLL(vpmrLa,:)<=.05) mean(p_nLL(pomrLa,:)<=.05)])

parameter='Puff';
exclude={'Pole';'Exclude';'Grooming'};
include={{'Light'};[-5 35]};
contPlot={MeanWpL,MeanPL,MeanAL};

figure('Position',[38 260 1413 525],'Color','White');Name='PuffLPopAllR';
plot_Pop_PSTH_wrapper3(parameter, DiscreteData,puff_triglength_limit,ppms,psth_time_before_trig,psth_time_after_trig,exclude,psth_safety_margin,psth_binsize,[],Name,{vpmrLa;pomrLa},{'VPM';'POm'},include,'latency',100,contPlot,{'meanbefore','scaleminmax','mono'},'baseline corr. whisker angle');
exportgraphics(gcf,[figdir Fig_no 'Pop_puff_responsive_' parameter '_L_instRateResp_' datestr(now,'yymmdd') '.pdf'],'BackgroundColor','none')


contPlot={MeanWp,MeanP,MeanAp};
exclude={'Pole';'Light';'Exclude';'Grooming'};
include=[];

figure('Position',[38 260 1413 525],'Color','White');Name='PuffnLPopAllR';
plot_Pop_PSTH_wrapper3(parameter, DiscreteData,puff_triglength_limit,ppms,psth_time_before_trig,psth_time_after_trig,exclude,psth_safety_margin,psth_binsize,[],Name,{vpmrLa;pomrLa},{'VPM';'POm'},include,'latency',100,contPlot,{'meanbefore','scaleminmax','mono'},'baseline corr. whisker angle');
exportgraphics(gcf,[figdir Fig_no 'Pop_puff_responsive_' parameter '_nL_instRateResp_' datestr(now,'yymmdd') '.pdf'],'BackgroundColor','none')

trigLen=cellfun(@(x) ones(1,numel(x))*30*20,trigsPL,'UniformOutput',0);
parameter=[trigsPL trigLen];
[H_PuffL,~,p_PuffL]=get_PSTH_pop(parameter, DiscreteData,puff_triglength_limit,ppms,rate_time_before_trig,rate_time_after_trig,exclude,psth_safety_margin,rate_binsize,isi,'left',[]);
p_PuffL=p_PuffL{1}(:,2);
sig_PuffL=p_PuffL<sig;
H_PL=cat(2,H_PuffL{:})';


%
SumStatLine=5;
SumStats(SumStatLine).condition_1='Puff no Light';
SumStats(SumStatLine).condition_2='Puff Light';
SumStats(SumStatLine).VPM_cell_IDs=vpmrLa';
SumStats(SumStatLine).POm_cell_IDs=pomrLa';
SumStats(SumStatLine).VPM_N=numel(vpmrLa);
SumStats(SumStatLine).POm_N=numel(pomrLa);
[N_animals]=get_N_animals(SumStats(SumStatLine).VPM_cell_IDs,SumStats(SumStatLine).POm_cell_IDs,RecDB);
SumStats(SumStatLine).VPM_N=numel(vpmrLa);
SumStats(SumStatLine).VPM_N_animals=N_animals(3);
SumStats(SumStatLine).POm_N=numel(pomrLa);
SumStats(SumStatLine).POm_N_animals=N_animals(4);
SumStats(SumStatLine).Total_N_animals=N_animals(1);
C=9;
SumStatTable{1,C}=SumStats(SumStatLine).VPM_N;
SumStatTable{2,C}=SumStats(SumStatLine).VPM_N_animals;
SumStatTable{17,C}=SumStats(SumStatLine).POm_N;
SumStatTable{18,C}=SumStats(SumStatLine).POm_N_animals;
SumStatTable{33,C}=SumStats(SumStatLine).Total_N_animals;
% 
% 
% figure('Position',[38 260 1413 525],'Color','White');Name='PuffLPop';
% plot_Pop_PSTH_wrapper3(parameter, DiscreteData,puff_triglength_limit,ppms,psth_time_before_trig,psth_time_after_trig,exclude,psth_safety_margin,psth_binsize,[],Name,{vpm;pom},{'VPM';'POm'},include,'latency',100,contPlot,{'meanbefore','scaleminmax','mono'},'baseline corr. whisker angle');
% exportgraphics(gcf,[figdir Fig_no 'Pop_all_' parameter '_L_instRateResp_' datestr(now,'yymmdd') '.pdf'],'BackgroundColor','none')
% figure('Position',[38 260 1413 525],'Color','White');Name='PuffLPop';
% plot_Pop_PSTH_wrapper3(parameter, DiscreteData,puff_triglength_limit,ppms,psth_time_before_trig,psth_time_after_trig,exclude,psth_safety_margin,psth_binsize,[],Name,{vpmL;pomL},{'VPM';'POm'},include,'latency',100,contPlot,{'meanbefore','scaleminmax','mono'},'baseline corr. whisker angle');
% exportgraphics(gcf,[figdir Fig_no 'Pop_all_Puff_nL_instRateResp_' datestr(now,'yymmdd') '.pdf'],'BackgroundColor','none')

figure('Position',[910   176   500   588])
[pI,pC,mS]=plot_ratecomp3(cat(3,H_P,H_PL),cat(2,sig_Puff,sig_PuffL),{'0';'1';'0';'1'},vpmrLa,{'nL Puff','L Puff'},'mean','Rate [Hz]');
title('nL/L Puff Response VPM', 'resp cells w light')
if write_figs;exportgraphics(gcf,[figdir Fig_no 'RatesLPuff_VPM_RespondingCellsWithLight.pdf'], 'BackgroundColor','none','ContentType','vector');end


SumStats(SumStatLine).VPM_mean_cond1_base=mS(1,1,1);
SumStats(SumStatLine).VPM_sem_cond1_base=mS(1,1,2);
SumStats(SumStatLine).VPM_mean_cond1_resp=mS(1,2,1);
SumStats(SumStatLine).VPM_sem_cond1_resp=mS(1,2,2);
SumStats(SumStatLine).VPM_p_cond1=pI(1);
SumStats(SumStatLine).VPM_mean_cond2_base=mS(2,1,1);
SumStats(SumStatLine).VPM_sem_cond2_base=mS(2,1,2);
SumStats(SumStatLine).VPM_mean_cond2_resp=mS(2,2,1);
SumStats(SumStatLine).VPM_sem_cond2_resp=mS(2,2,2);
SumStats(SumStatLine).VPM_p_cond2=pI(2);
SumStats(SumStatLine).VPM_p_base_cond12=pC(1);
SumStats(SumStatLine).VPM_p_resp_cond12=pC(2);



figure('Position',[910   176   500   588])
[pI,pC,mS]=plot_ratecomp3(cat(3,H_P,H_PL),cat(2,sig_Puff,sig_PuffL),{'0';'1';'0';'1'},pomrLa,{'nL Puff','L Puff'},'mean','Rate [Hz]');
title('nL/L Puff Response POm', 'resp cells w light')
if write_figs;exportgraphics(gcf,[figdir Fig_no 'RatesLPuff_POm_RespondingCellsWithLight.pdf'], 'BackgroundColor','none','ContentType','vector');end


SumStats(SumStatLine).POm_mean_cond1_base=mS(1,1,1);
SumStats(SumStatLine).POm_sem_cond1_base=mS(1,1,2);
SumStats(SumStatLine).POm_mean_cond1_resp=mS(1,2,1);
SumStats(SumStatLine).POm_sem_cond1_resp=mS(1,2,2);
SumStats(SumStatLine).POm_p_cond1=pI(1);
SumStats(SumStatLine).POm_mean_cond2_base=mS(2,1,1);
SumStats(SumStatLine).POm_sem_cond2_base=mS(2,1,2);
SumStats(SumStatLine).POm_mean_cond2_resp=mS(2,2,1);
SumStats(SumStatLine).POm_sem_cond2_resp=mS(2,2,2);
SumStats(SumStatLine).POm_p_cond2=pI(2);
SumStats(SumStatLine).POm_p_base_cond12=pC(1);
SumStats(SumStatLine).POm_p_resp_cond12=pC(2);

pRc=nan(2,2);
for c=1:2
pRc(1,c)=ranksum(H_P(vpmrLa,c,1),H_P(pomrLa,c,1));
end
for c=1:2
pRc(2,c)=ranksum(H_PL(vpmrLa,c,1),H_PL(pomrLa,c,1));
end

SumStats(SumStatLine).VPMPOm_p_cond1_base=pRc(1,1);
SumStats(SumStatLine).VPMPOm_p_cond2_base=pRc(2,1);
SumStats(SumStatLine).VPMPOm_p_cond1_resp=pRc(1,2);
SumStats(SumStatLine).VPMPOm_p_cond2_resp=pRc(2,2);

%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p_PuffQWr=sig_test_comp_trig({trigsPnW;trigsPW},DiscreteData,rate_bins,2,'both',ppms);
sig_PuffQWr=p_PuffQWr<sig;
p_PuffQWb=sig_test_comp_trig({trigsPnW;trigsPW},DiscreteData,rate_bins,1,'both',ppms);
sig_PuffQWb=p_PuffQWb<sig;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tempD=cat(2,log2(H_P(:,2)./H_P(:,1)),log2(H_PL(:,2)./H_PL(:,1)));tempD(any(isnan(tempD),2),:)=nan;
figure('Position',[ 910   176   500   588])
plot_ratecomp2(tempD,sig_Touch,{'noLPuff','LPuff'},{vpmrLa;pomrLa},{'VPM';'POm'},'mean','Rate fold change');title('nL/L Puff fold change Response, responsive cells')
if write_figs;exportgraphics (gcf,[figdir Fig_no 'foldchangeRates_nLL_RespondingCells.pdf'], 'BackgroundColor','none','ContentType','vector');end


tempD=cat(2,mod_depth(H_P),mod_depth(H_PL));
figure('Position',[ 910   176   500   588])
[pI,pC,mS,p0]=plot_ratecomp2(tempD,true(size(sig_Touch)),{'noLPuff','LPuff'},{vpmrLa;pomrLa},{'VPM';'POm'},'mean','mod depth');title('nL/L Puff modulation, responsive cells')
 if write_figs;exportgraphics (gcf,[figdir Fig_no 'modulation_nLL_RespondingCellsLog.pdf'], 'BackgroundColor','none','ContentType','vector');end


tempD=cat(2,H_P(:,2)./H_P(:,1),H_PL(:,2)./H_PL(:,1));tempD(any(isnan(tempD),2),:)=nan;
figure('Position',[ 910   176   500   588])
plot_ratecomp2(tempD,true(size(sig_Touch)),{'noLPuff','LPuff'},{vpmrLa;pomrLa},{'VPM';'POm'},'mean','Rate fold change');title('nL/L Puff fold change Response, responsive cells')
 set(gca,'YScale','log')
 set(gca,'YTick',[0.05 0.1 0.25 0.5 1 2 4 8 16 32 50 100])
 ylim([0.05 120])

if write_figs;exportgraphics (gcf,[figdir Fig_no 'foldchangeRates_nLL_RespondingCellsLOG.pdf'], 'BackgroundColor','none','ContentType','vector');end


 
SumStats(SumStatLine).VPM_mean_modulation_cond1=mS(1,1,1);
SumStats(SumStatLine).VPM_sem_modulation_cond1=mS(1,1,2);
SumStats(SumStatLine).VPM_p_modulation_cond10=p0(1,1);
SumStats(SumStatLine).VPM_mean_modulation_cond2=mS(1,2,1);
SumStats(SumStatLine).VPM_sem_modulation_cond2=mS(1,2,2);
SumStats(SumStatLine).VPM_p_modulation_cond20=p0(1,2);
SumStats(SumStatLine).VPM_p_modulation_cond12=pI(1);

SumStats(SumStatLine).POm_mean_modulation_cond1=mS(2,1,1);
SumStats(SumStatLine).POm_sem_modulation_cond1=mS(2,1,2);
SumStats(SumStatLine).POm_p_modulation_cond10=p0(2,1);
SumStats(SumStatLine).POm_mean_modulation_cond2=mS(2,2,1);
SumStats(SumStatLine).POm_sem_modulation_cond2=mS(2,2,2);
SumStats(SumStatLine).POm_p_modulation_cond20=p0(2,2);
SumStats(SumStatLine).POm_p_modulation_cond12=pI(2);

SumStats(SumStatLine).VPMPOm_p_modulation_cond1=pC(1);
SumStats(SumStatLine).VPMPOm_p_modulation_cond2=pC(2);


% Supp Figure 5
Fig_no='EDFig5_';
lat_binsize=1;
lat_time_before_trig=0;
lat_time_after_trig=75;
lat_safety_margin=lat_time_after_trig;


puff_triglength_limit=[28 32];%ms
% no light
trigLen=cellfun(@(x) ones(1,numel(x))*30*20,trigsPL,'UniformOutput',0);
parameter=[trigsPL trigLen];
tempDD=struct2cell(DiscreteData);
spikes=squeeze(tempDD(1,1,:));clear temp*
[RasterHPL_nW,P_trig]=get_raster_pop_no_self_safety(parameter,spikes, DiscreteData,triglength_limit,ppms,lat_time_before_trig,lat_time_after_trig,exclude,lat_safety_margin,lat_binsize,isi,[]);
%light
trigLen=cellfun(@(x) ones(1,numel(x))*30*20,trigsP,'UniformOutput',0);
parameter=[trigsP trigLen];
[RasterHPL_W,T_trig]=get_raster_pop_no_self_safety(parameter,spikes, DiscreteData,triglength_limit,ppms,lat_time_before_trig,lat_time_after_trig,exclude,lat_safety_margin,lat_binsize,isi,[]);
first_spike_times=cell(numel(RasterHPL_nW),2);
fsl_p=ones(numel(RasterHPL_nW),1);
for n=1:numel(RasterHPL_nW)
    [a,ia,ib]=unique(RasterHPL_nW{n}(:,2));
    first_spike_times{n,1}=nan(size(P_trig{n},1),1);
    first_spike_times{n,1}(a)=RasterHPL_nW{n}(ia,1);
    [a,ia,ib]=unique(RasterHPL_W{n}(:,2));
    first_spike_times{n,2}=nan(size(T_trig{n},1),1);
    first_spike_times{n,2}(a)=RasterHPL_W{n}(ia,1);
if all(cellfun(@(x) any(~isnan(x)),first_spike_times(n,:)))
     [fsl_p(n)]=ranksum(first_spike_times{n,1},first_spike_times{n,2});
end
end
%
FSLavg=cellfun(@(x) mean(x,'omitnan'),first_spike_times);
FSLstd=cellfun(@(x) std(x,0,'omitnan'),first_spike_times);
FSLmed=cellfun(@(x) median(x,'omitnan'),first_spike_times);
FSLiqr=cellfun(@(x) iqr(x),first_spike_times);
FSLany=cellfun(@(x) mean(~isnan(x)),first_spike_times);
FSLsig=true(size(FSLavg));%repmat(fsl_p<.05,1,2);%
figure('Position',[ 910   176   500   588])
[pI,pC,mS]=plot_ratecomp2(FSLavg,FSLsig,{'nLPuff','LPuff'},{vpmrLa;pomrLa},{'VPM';'POm'},'mean','first spike latency (ms)');title('mean first spike latency')
if write_figs;exportgraphics (gcf,[figdir Fig_no 'nLLPuff_firstspikelatency.pdf'],'BackgroundColor','none','ContentType','vector');end

SumStats(SumStatLine).VPM_mean_fsl_cond1=mS(1,1,1);
SumStats(SumStatLine).VPM_sem_fsl_cond1=mS(1,1,2);
SumStats(SumStatLine).VPM_mean_fsl_cond2=mS(1,2,1);
SumStats(SumStatLine).VPM_sem_fsl_cond2=mS(1,2,2);
SumStats(SumStatLine).VPM_p_fsl_cond12=pI(1);

SumStats(SumStatLine).POm_mean_fsl_cond1=mS(2,1,1);
SumStats(SumStatLine).POm_sem_fsl_cond1=mS(2,1,2);
SumStats(SumStatLine).POm_mean_fsl_cond2=mS(2,2,1);
SumStats(SumStatLine).POm_sem_fsl_cond2=mS(2,2,2);
SumStats(SumStatLine).POm_p_fsl_cond12=pI(2);

SumStats(SumStatLine).VPMPOm_p_fsl_cond1=pC(1);
SumStats(SumStatLine).VPMPOm_p_fsl_cond2=pC(2);

XCond=9;
SumStatTable{3,XCond}=SumStats(SumStatLine).VPM_mean_cond1_base;
SumStatTable{4,XCond}=SumStats(SumStatLine).VPM_sem_cond1_base;
SumStatTable{6,XCond}=SumStats(SumStatLine).VPM_mean_cond1_resp;
SumStatTable{7,XCond}=SumStats(SumStatLine).VPM_sem_cond1_resp;
SumStatTable{9,XCond}=SumStats(SumStatLine).VPM_p_cond1;
SumStatTable{3,XCond+1}=SumStats(SumStatLine).VPM_mean_cond2_base;
SumStatTable{4,XCond+1}=SumStats(SumStatLine).VPM_sem_cond2_base;
SumStatTable{6,XCond+1}=SumStats(SumStatLine).VPM_mean_cond2_resp;
SumStatTable{7,XCond+1}=SumStats(SumStatLine).VPM_sem_cond2_resp;
SumStatTable{9,XCond+1}=SumStats(SumStatLine).VPM_p_cond2;
SumStatTable{5,XCond}=SumStats(SumStatLine).VPM_p_base_cond12;
SumStatTable{8,XCond}=SumStats(SumStatLine).VPM_p_resp_cond12;
SumStatTable{16+3,XCond}=SumStats(SumStatLine).POm_mean_cond1_base;
SumStatTable{16+4,XCond}=SumStats(SumStatLine).POm_sem_cond1_base;
SumStatTable{16+6,XCond}=SumStats(SumStatLine).POm_mean_cond1_resp;
SumStatTable{16+7,XCond}=SumStats(SumStatLine).POm_sem_cond1_resp;
SumStatTable{16+9,XCond}=SumStats(SumStatLine).POm_p_cond1;
SumStatTable{16+3,XCond+1}=SumStats(SumStatLine).POm_mean_cond2_base;
SumStatTable{16+4,XCond+1}=SumStats(SumStatLine).POm_sem_cond2_base;
SumStatTable{16+6,XCond+1}=SumStats(SumStatLine).POm_mean_cond2_resp;
SumStatTable{16+7,XCond+1}=SumStats(SumStatLine).POm_sem_cond2_resp;
SumStatTable{16+9,XCond+1}=SumStats(SumStatLine).POm_p_cond2;
SumStatTable{16+5,XCond}=SumStats(SumStatLine).POm_p_base_cond12;
SumStatTable{16+8,XCond}=SumStats(SumStatLine).POm_p_resp_cond12;
SumStatTable{10,XCond}=SumStats(SumStatLine).VPM_mean_modulation_cond1;
SumStatTable{11,XCond}=SumStats(SumStatLine).VPM_sem_modulation_cond1;
SumStatTable{12,XCond}=SumStats(SumStatLine).VPM_p_modulation_cond10;
SumStatTable{10,XCond+1}=SumStats(SumStatLine).VPM_mean_modulation_cond2;
SumStatTable{11,XCond+1}=SumStats(SumStatLine).VPM_sem_modulation_cond2;
SumStatTable{12,XCond+1}=SumStats(SumStatLine).VPM_p_modulation_cond20;
SumStatTable{13,XCond}=SumStats(SumStatLine).VPM_p_modulation_cond12;
SumStatTable{16+10,XCond}=SumStats(SumStatLine).POm_mean_modulation_cond1;
SumStatTable{16+11,XCond}=SumStats(SumStatLine).POm_sem_modulation_cond1;
SumStatTable{16+12,XCond}=SumStats(SumStatLine).POm_p_modulation_cond10;
SumStatTable{16+10,XCond+1}=SumStats(SumStatLine).POm_mean_modulation_cond2;
SumStatTable{16+11,XCond+1}=SumStats(SumStatLine).POm_sem_modulation_cond2;
SumStatTable{16+12,XCond+1}=SumStats(SumStatLine).POm_p_modulation_cond20;
SumStatTable{16+13,XCond}=SumStats(SumStatLine).POm_p_modulation_cond12;
SumStatTable{36,XCond}=SumStats(SumStatLine).VPMPOm_p_modulation_cond1;
SumStatTable{36,XCond+1}=SumStats(SumStatLine).VPMPOm_p_modulation_cond2;
SumStatTable{34,XCond}=SumStats(SumStatLine).VPMPOm_p_cond1_base;
SumStatTable{34,XCond+1}=SumStats(SumStatLine).VPMPOm_p_cond2_base;
SumStatTable{35,XCond}=SumStats(SumStatLine).VPMPOm_p_cond1_resp;
SumStatTable{35,XCond+1}=SumStats(SumStatLine).VPMPOm_p_cond2_resp;
SumStatTable{14,XCond}=SumStats(SumStatLine).VPM_mean_fsl_cond1;
SumStatTable{15,XCond}=SumStats(SumStatLine).VPM_sem_fsl_cond1;
SumStatTable{14,XCond+1}=SumStats(SumStatLine).VPM_mean_fsl_cond2;
SumStatTable{15,XCond+1}=SumStats(SumStatLine).VPM_sem_fsl_cond2;
SumStatTable{16,XCond}=SumStats(SumStatLine).VPM_p_fsl_cond12;
SumStatTable{16+14,XCond}=SumStats(SumStatLine).POm_mean_fsl_cond1;
SumStatTable{16+15,XCond}=SumStats(SumStatLine).POm_sem_fsl_cond1;
SumStatTable{16+14,XCond+1}=SumStats(SumStatLine).POm_mean_fsl_cond2;
SumStatTable{16+15,XCond+1}=SumStats(SumStatLine).POm_sem_fsl_cond2;
SumStatTable{16+16,XCond}=SumStats(SumStatLine).POm_p_fsl_cond12;
SumStatTable{37,XCond}=SumStats(SumStatLine).VPMPOm_p_fsl_cond1;
SumStatTable{37,XCond+1}=SumStats(SumStatLine).VPMPOm_p_fsl_cond2;

%%

writetable(SumStatTable,'C:\Data\Whisker\sumstats3.xlsx')
%%

singleCellTableNames={'NeuronID','Nucleus','AnimalID','VGAT','','','','','','','','','','','','','','','',''}


