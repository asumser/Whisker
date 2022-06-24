%% load stuff
DrivePath='G:\Meine Ablage\';
DataPath='G:\Meine Ablage\Arbeit\Projects\Database\EphysData\';
load([DrivePath 'Arbeit\Projects\Database\ephys_database.mat'])
RecDB.Properties.RowNames=strcat(RecDB.AnimalName,{'_'},datestr(RecDB.StartTime,'HH_MM_SS'));
load([DataPath 'DiscreteData.mat'])
ppms=20;
load('D:\Whisker\Wdec\Wdec_by_field.mat', 'Angle','Phase','Amp_extrema','Setpoint_extrema','Curvature')
Amp=Amp_extrema;clear Amp_extrema
Setp=Setpoint_extrema;clear Setpoint_extrema
minFillQ=4;
% vpm=find(RecDB.DefinedNucleus=='VPM' & RecDB.FillQuality<=minFillQ & RecDB.RecordingType=='juxta');
% pom=find(RecDB.DefinedNucleus=='POm' & RecDB.FillQuality<=minFillQ & RecDB.RecordingType=='juxta');
%% fig stuff
figdir=['C:\Data\Whisker\Paperfigs_matlab\' datestr(now,'YYmmDD') '\'];
set(0,'defaultfigurecolor',[1 1 1]);
set(0, 'DefaultAxesBox', 'off')
set(0,'defaultAxesFontName', 'Arial')
set(0,'defaultTextFontName', 'Arial')
write_figs=true;
if write_figs
    mkdir(figdir)
end
%% settings
psth_time_before_trig=25;%ms
psth_time_after_trig=75;%ms
psth_binsize=1;%ms
psth_safety_margin=psth_time_after_trig;%ms
rate_time_before_trig=50;%ms
rate_time_after_trig=50;%ms
rate_binsize=50;%ms
sig=.05;
isi=10;
puff_triglength_limit=[28 32];%ms
touch_triglength_limit=[0 inf];%ms
phasebinN=32;
amp_thresh=3;

rate_bins=-rate_time_before_trig:rate_binsize:rate_time_after_trig;
psth_bins=-psth_time_before_trig:psth_binsize:psth_time_after_trig;
psth_relativeTime=psth_bins(1:end-1)+diff(psth_bins);
%%
% exclude={'Pole';'Light';'Exclude';'Grooming'};
% [trigs_rep,trigs_rep_first]=get_reg_puff_trigs(DiscreteData,ppms,exclude,psth_safety_margin);
% [PuffV]=mk_vector_from_start_length(DiscreteData,'Puff');
% %%
% rep_time_before_trig=200;
% rep_time_after_trig=1200;
% rep_binsize=2.5;
% rep_bins=-rep_time_before_trig:rep_binsize:rep_time_after_trig;
% rep_relativeTime=rep_bins(1:end-1)+diff(rep_bins);
% [MeanWpr]=get_trig_cont_pop(trigs_rep_first,rep_relativeTime,Angle,20,DiscreteData);
% [MeanPr]=get_trig_cont_pop(trigs_rep_first,rep_relativeTime,PuffV,1,DiscreteData);
% %
%     figure('Position',[38 260 1413 525],'Color','White');
% contPlot={MeanWpr,MeanPr};
% Name='regular puff';
% [P_PSTH_instrate]=plot_Pop_PSTH_wrapper3(trigs_rep_first, DiscreteData,puff_triglength_limit,ppms,rep_time_before_trig,rep_time_after_trig,exclude,psth_safety_margin,rep_binsize,[],Name,{vpmra},{'VPM'},[],'latency',-inf,contPlot,{'meanbefore','scaleminmax'},'baseline corr. whisker angle');
% if write_figs;exportgraphics(gcf,[figdir 'RegularPuff.pdf'], 'BackgroundColor','none','ContentType','vector');end
% 

%% calculate supporting data
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

%% Figure 2
Fig_no='Fig2_';
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

% extFigure 2: puff psth incl notouch
Name='Puff_allResp';
figure('Position',[38 260 1413 525],'Color','White');
[P_PSTH_instrate]=plot_Pop_PSTH_wrapper3(parameter, DiscreteData,puff_triglength_limit,ppms,psth_time_before_trig,psth_time_after_trig,exclude,psth_safety_margin,psth_binsize,[],Name,{vpmra;pomra},{'VPM';'POm'},include,'latency',-inf,contPlot,{'meanbefore','scaleminmax','mono'},'baseline corr. whisker angle');
if write_figs;exportgraphics(gcf,[figdir 'sFig2_Pop_including_notouchneurons_' parameter '_instRateResp_' datestr(now,'yymmdd') '.pdf'],'BackgroundColor','none');end
%puff rates
parameter='Puff';
exclude={'Pole';'Light';'Exclude';'Grooming'};include=[];
[H_Puff,~,p_Puff,~,p_trigs,~,p_trigL]=get_PSTH_pop(parameter, DiscreteData,puff_triglength_limit,ppms,rate_time_before_trig,rate_time_after_trig,exclude,psth_safety_margin,rate_binsize,isi,'left',[]);
p_Puff=p_Puff{1}(:,2);
sig_Puff=p_Puff<sig;
H_P=cat(2,H_Puff{:})';
figure('Position',[910   176   500   588])
plot_ratecomp2(H_P,sig_Puff,{'before Puff';'after Puff'},{vpm;pom},{'VPM';'POm'},'meansd','Rate [Hz]');title('Puff Response, all cells')
if write_figs;exportgraphics(gcf,[figdir Fig_no 'RatesPuffAllCells.pdf'], 'BackgroundColor','none','ContentType','vector');end
figure('Position',[ 910   176   500   588])
plot_ratecomp2(H_P,sig_Puff,{'before Puff';'after Puff'},{vpmra;pomra},{'VPM';'POm'},'meansd','Rate [Hz]');title('Puff Response, responsive cells')
if write_figs;exportgraphics(gcf,[figdir Fig_no 'RatesPuffRespondingCells.pdf'], 'BackgroundColor','none','ContentType','vector');end
figure('Position',[ 910   176   500   588])
plot_ratecomp2(H_P,sig_Puff,{'before Puff';'after Puff'},{vpmr;pomr},{'VPM';'POm'},'meansd','Rate [Hz]');title('Puff Response, responsive cells(haveTouch)')
if write_figs;exportgraphics(gcf,[figdir Fig_no 'RatesPuffRespondingCellsWithTouch.pdf'], 'BackgroundColor','none','ContentType','vector');end

%
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


%touch rates
[~,~,~,~,ttrigs,zeta_T,t_trigL]=get_PSTH_pop(parameter, DiscreteData,touch_triglength_limit,ppms,rate_time_before_trig,rate_time_after_trig,exclude,psth_safety_margin,rate_binsize,isi,'left',[]);
[H_Touch,~,p_Touch]=get_PSTH_pop(parameter, DiscreteData,touch_triglength_limit,ppms,rate_time_before_trig,rate_time_after_trig,exclude,psth_safety_margin,rate_binsize,isi,'left',[]);
p_Touch=p_Touch{1}(:,2);
sig_Touch=p_Touch<sig;
H_T=cat(2,H_Touch{:})';
figure('Position',[ 910   176   500   588])
plot_ratecomp2(H_T,sig_Touch,{'before Touch';'after Touch'},{vpmr;pomr},{'VPM';'POm'},'meansd','Rate [Hz]');title('Touch Response, responsive cells')
if write_figs;exportgraphics(gcf,[figdir Fig_no 'RatesTouchRespondingCells.pdf'],'BackgroundColor','none','ContentType','vector');end
% sigT=nan(numel(zeta_T),1);
% for n=1:numel(zeta_T)
% sigT(n)=zeta_T{n}.dblP;
% end
%%
%touch puff comp
p_PuffPTr=sig_test_comp_trig({p_trigs;ttrigs},DiscreteData,rate_bins,2,'both',ppms);
p_PuffPTb=sig_test_comp_trig({p_trigs;ttrigs},DiscreteData,rate_bins,1,'both',ppms);
sig_PTr=p_PuffPTr<sig;
sig_PTb=p_PuffPTb<sig;
tempD=cat(2,diff(H_P,1,2),diff(H_T,1,2));
figure('Position',[ 910   176   500   588])
plot_ratecomp2(tempD,sig_PTr,{'Puff','Touch'},{vpmr;pomr},{'VPM';'POm'},'meansd','deltaRate [Hz]');title('diffTouch/Puff Response, responsive cells')
if write_figs;exportgraphics (gcf,[figdir Fig_no 'diffRatesTouchVsPuffRespondingCells.pdf'],'BackgroundColor','none','ContentType','vector');end
tempD=cat(2,H_P(:,1),H_T(:,1));
figure('Position',[ 910   176   500   588])
plot_ratecomp2(tempD,sig_PTb,{'Puff','Touch'},{vpmr;pomr},{'VPM';'POm'},'meansd','Rate [Hz]');title('baseline before Touch/Puff activity, responsive cells')
if write_figs;exportgraphics (gcf,[figdir Fig_no 'baseRatesTouchVsPuffRespondingCells.pdf'],'BackgroundColor','none','ContentType','vector');end
tempD=cat(2,H_P(:,2),H_T(:,2));
figure('Position',[ 910   176   500   588])
plot_ratecomp2(tempD,sig_PTr,{'Puff','Touch'},{vpmr;pomr},{'VPM';'POm'},'meansd','Rate [Hz]');title('Touch/Puff Response, responsive cells')
if write_figs;exportgraphics (gcf,[figdir Fig_no 'absRatesTouchVsPuffRespondingCells.pdf'], 'BackgroundColor','none','ContentType','vector');end
tempD=log2(cat(2,H_P(:,2)./H_P(:,1),H_T(:,2)./H_T(:,1)));
figure('Position',[ 910   176   500   588])
plot_ratecomp2(tempD,sig_PTr,{'Puff','Touch'},{vpmr;pomr},{'VPM';'POm'},'meansd','log2 Rate fold change');title('Touch/Puff fold change Response, responsive cells')
if write_figs;exportgraphics (gcf,[figdir Fig_no 'foldchangeRatesTouchVsPuffRespondingCells.pdf'], 'BackgroundColor','none','ContentType','vector');end


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
cell_example=pomr(9);
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

%%
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

% set(gca,'CameraViewAngle',8.2541,'CameraTarget',[1647       -1821       -3110])
% grid off
%set(gca,'XColor', 'none','YColor','none','ZColor','none')
%%
% raw_data_example=pomr(9);
raw_data_example=vpmr(4);
tE=(1:DiscreteData(raw_data_example).LengthInd)/(ppms*1000);
tW=downsample(tE,20);
load(['D:\AnalysisMatFiles\' RecDB.Properties.RowNames{raw_data_example} '\' RecDB.Properties.RowNames{raw_data_example} 'analysis.mat'],'filteredResponse')
%
figure
tiledlayout(2,1)
a(1)=nexttile(1);
plot(tE,filteredResponse.data,'k');box off;hold on
plot(DiscreteData(raw_data_example).PuffStart/(ppms*1000),max(filteredResponse.data),'b.');
plot(DiscreteData(raw_data_example).TouchStart/(ppms*1000),max(filteredResponse.data),'r.');
a(2)=nexttile(2);
plot(tW,Angle{raw_data_example},'Color',[.8 .8 .8]);box off;hold on
plot(DiscreteData(raw_data_example).PuffStart/(ppms*1000),max(Angle{raw_data_example}),'b.');
plot(DiscreteData(raw_data_example).TouchStart/(ppms*1000),max(Angle{raw_data_example}),'r.');

linkaxes(a,'x')
%%
exportgraphics(gcf,[figdir Fig_no 'raw_data_puff_cell' num2str(raw_data_example) '.pdf'],'BackgroundColor','none','ContentType','vector')
exportgraphics(gcf,[figdir Fig_no 'raw_data_touch_cell' num2str(raw_data_example) '.pdf'],'BackgroundColor','none','ContentType','vector')

%% Figure 3
Fig_no='Fig3_';
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
figure('Position',[38 260 1413 525],'Color','White');Name='Puff noWh allResp';
plot_Pop_PSTH_wrapper3(parameter, DiscreteData,puff_triglength_limit,ppms,psth_time_before_trig,psth_time_after_trig,exclude,psth_safety_margin,psth_binsize,[],Name,{vpmra;pomra},{'VPM';'POm'},include,'latency',-inf,contPlot,{'meanbefore','scaleminmax','mono'},'baseline corr. whisker angle');
if write_figs;exportgraphics(gcf,[figdir Fig_no 'Pop_including_notouchneurons_Puff_nW_instRateResp_' datestr(now,'yymmdd') '.pdf'],'BackgroundColor','none','ContentType','vector');end
fig2a2=figure('Position',[-1731         836        1322         171],'Color','White');
Name='Puff nW trig whisking';
plot_wh_trace_with_sem(MeanWpnw,SEM_Wpnw,{vpmr;pomr},{'VPM';'POm'},'wh angle',psth_relativeTime,Name)
if write_figs;exportgraphics(fig2a2,[figdir Fig_no 'Pop_Puffnw_wh_with_sem_' datestr(now,'yymmdd') '.pdf'],'BackgroundColor','none');end



[H_PuffnW,~,p_PuffnW]=get_PSTH_pop(parameter, DiscreteData,puff_triglength_limit,ppms,rate_time_before_trig,rate_time_after_trig,exclude,psth_safety_margin,rate_binsize,isi,'left',[]);
H_PnW=cat(2,H_PuffnW{:})';
p_PuffnW=p_PuffnW{1}(:,2);
sig_PuffnW=p_PuffnW<sig;
figure('Position',[ 910   176   500   588])
plot_ratecomp2(H_PnW,sig_PuffnW,{'before Puff';'after Puff'},{vpmra;pomra},{'VPM';'POm'},'meansd','Rate [Hz]');title('Q Puff Response, responsive cells')
if write_figs;exportgraphics (gcf, [figdir Fig_no 'RatesPuff_nW_RespondingCells.pdf'],'BackgroundColor','none','ContentType','vector');end
figure('Position',[ 910   176   500   588])
plot_ratecomp2(H_PnW,sig_PuffnW,{'before Puff';'after Puff'},{vpmr;pomr},{'VPM';'POm'},'meansd','Rate [Hz]');title('Q Puff Response, responsive cells(haveTouch)')
if write_figs;exportgraphics (gcf, [figdir Fig_no 'RatesPuff_nW_RespondingCellsWithTouch.pdf'],'BackgroundColor','none','ContentType','vector');end

%
% Whisking
contPlot={MeanWpw,MeanPW,MeanApw};
trigLen=cellfun(@(x) ones(1,numel(x))*30*20,trigsPW,'UniformOutput',0);
parameter=[trigsPW trigLen];
figure('Position',[38 260 1413 525],'Color','White');Name='Puff Wh';
plot_Pop_PSTH_wrapper3(parameter, DiscreteData,puff_triglength_limit,ppms,psth_time_before_trig,psth_time_after_trig,exclude,psth_safety_margin,psth_binsize,[],Name,{vpmr;pomr},{'VPM';'POm'},include,'latency',-inf,contPlot,{'meanbefore','scaleminmax','mono'},'baseline corr. whisker angle');
if write_figs;exportgraphics(gcf,[figdir Fig_no 'Pop_Puff_W_instRateResp_' datestr(now,'yymmdd') '.pdf'],'BackgroundColor','none');end
figure('Position',[38 260 1413 525],'Color','White');Name='Puff Wh_allResp';
plot_Pop_PSTH_wrapper3(parameter, DiscreteData,puff_triglength_limit,ppms,psth_time_before_trig,psth_time_after_trig,exclude,psth_safety_margin,psth_binsize,[],Name,{vpmra;pomra},{'VPM';'POm'},include,'latency',-inf,contPlot,{'meanbefore','scaleminmax','mono'},'baseline corr. whisker angle');
if write_figs;exportgraphics(gcf,[figdir Fig_no 'Pop_including_notouchneurons_Puff_W_instRateResp_' datestr(now,'yymmdd') '.pdf'],'BackgroundColor','none');end
fig2a2=figure('Position',[-1731         836        1322         171],'Color','White');
Name='Puff W trig whisking';
plot_wh_trace_with_sem(MeanWpw,SEM_Wpw,{vpmr;pomr},{'VPM';'POm'},'wh angle',psth_relativeTime,Name)
if write_figs;exportgraphics(fig2a2,[figdir Fig_no 'Pop_Puffw_wh_with_sem_' datestr(now,'yymmdd') '.pdf'],'BackgroundColor','none');end


[H_PuffW,~,p_PuffW]=get_PSTH_pop(parameter, DiscreteData,puff_triglength_limit,ppms,rate_time_before_trig,rate_time_after_trig,exclude,psth_safety_margin,rate_binsize,isi,'left',[]);
H_PW=cat(2,H_PuffW{:})';
p_PuffW=p_PuffW{1}(:,2);
sig_PuffW=p_PuffW<sig;
figure('Position',[ 910   176   500   588])
plot_ratecomp2(H_PW,sig_PuffW,{'before Puff';'after Puff'},{vpmra;pomra},{'VPM';'POm'},'meansd','Rate [Hz]');title('W Puff Response, responsive cells')
if write_figs;exportgraphics (gcf, [figdir Fig_no 'RatesPuff_W_RespondingCells.pdf'],'BackgroundColor','none','ContentType','vector');end
figure('Position',[ 910   176   500   588])
plot_ratecomp2(H_PW,sig_PuffW,{'before Puff';'after Puff'},{vpmr;pomr},{'VPM';'POm'},'meansd','Rate [Hz]');title('W Puff Response, responsive cells(haveTouch)')
if write_figs;exportgraphics (gcf, [figdir Fig_no 'RatesPuff_W_RespondingCellsWithTouch.pdf'],'BackgroundColor','none','ContentType','vector');end
%
%
p_PuffQWr=sig_test_comp_trig({trigsPnW;trigsPW},DiscreteData,rate_bins,2,'both',ppms);
sig_PuffQWr=p_PuffQWr<sig;
p_PuffQWb=sig_test_comp_trig({trigsPnW;trigsPW},DiscreteData,rate_bins,1,'both',ppms);
sig_PuffQWb=p_PuffQWb<sig;
tempD=cat(2,diff(H_PnW,1,2),diff(H_PW,1,2));
figure('Position',[ 910   176   500   588])
plot_ratecomp2(tempD,sig_PuffQWr,{'QPuff','WPuff'},{vpmr;pomr},{'VPM';'POm'},'meansd','deltaRate [Hz]');title('diffQ/W Puff Response, responsive cells')
if write_figs;exportgraphics (gcf,[figdir Fig_no 'diffRatesQWPuffRespondingCells.pdf'],'BackgroundColor','none','ContentType','vector');end
tempD=cat(2,H_PnW(:,1),H_PW(:,1));
figure('Position',[ 910   176   500   588])
plot_ratecomp2(tempD,sig_PuffQWb,{'QPuff','WPuff'},{vpmr;pomr},{'VPM';'POm'},'meansd','Rate [Hz]');title('baseline activity before Q/W Puff, responsive cells')
if write_figs;exportgraphics (gcf,[figdir Fig_no 'baseRatesQWPuffRespondingCells.pdf'],'BackgroundColor','none','ContentType','vector');end
tempD=cat(2,H_PnW(:,2),H_PW(:,2));
figure('Position',[ 910   176   500   588])
plot_ratecomp2(tempD,sig_PuffQWr,{'QPuff','WPuff'},{vpmr;pomr},{'VPM';'POm'},'meansd','Rate [Hz]');title('Q/W Puff Response, responsive cells')
if write_figs;exportgraphics (gcf,[figdir Fig_no 'absRatesQWPuffRespondingCells.pdf'], 'BackgroundColor','none','ContentType','vector');end
tempD=cat(3,cat(2,H_PnW(:,2)./H_PnW(:,1),H_PW(:,2)./H_PW(:,1)),log2(cat(2,H_PnW(:,2)./H_PnW(:,1),H_PW(:,2)./H_PW(:,1))));
figure('Position',[ 910   176   500   588])
plot_ratecomp2(tempD,sig_PuffQWr,{'QPuff','WPuff'},{vpmr;pomr},{'VPM';'POm'},'meansd','log2 Rate fold change');title('Q/W Puff fold change Response, responsive cells')
if write_figs;exportgraphics (gcf,[figdir Fig_no 'foldchangeRatesQWRespondingCells.pdf'], 'BackgroundColor','none','ContentType','vector');end

%% Ratecomp Whisking
whisking_time_before_trig=300;%ms
whisking_time_after_trig=300;%ms
whisking_binsize=300;%ms
whsiking_safety_margin=whisking_time_after_trig;
isi=10;
sig=.05;
parameter='Whisking';
exclude={'Puff';'Pole';'Light';'Exclude';'Grooming'};
triglength_limit=[whisking_time_after_trig inf];%ms
[H_W,~,p_W]=get_PSTH_pop(parameter, DiscreteData,triglength_limit,ppms,whisking_time_before_trig,whisking_time_after_trig,exclude,whsiking_safety_margin,whisking_binsize,isi,'left',[]);
p_W=p_W{1}(:,2);
sig_W=p_W<sig;
H_Wc=cat(2,H_W{:})';
figure('Position',[ 910   176   500   588])
plot_ratecomp2(H_Wc,sig_W,{'preWhisk','Whisking'},{vpmr;pomr},{'VPM';'POm'},'meansd','Rate (Hz)');title('Whisking Rate, responsive cells')
if write_figs;exportgraphics (gcf,[figdir Fig_no 'WhiskingRatesRespondingCells.pdf'], 'BackgroundColor','none','ContentType','vector');end

%%
Fig_no='Fig4_';
exclude={'Pole';'Light';'Exclude';'Grooming';'Puff'};
include=[];
safety_margin=[100 200];
min_t=500;%minimal time in bin

phasebinN=16;
ampbins=[0:3:30 inf];
abssetpbins=[0:3:30 inf];
absangbins=[0:3:40 inf];
setpbins=[-90 -25:3:30 inf];
angbins=[-90 -30:3:45 inf];
wbins={phasebinN;ampbins;abssetpbins;absangbins;setpbins;angbins};
[Num_bin_occ,rates_to_plot,plot_x,name_vars,PhaseMod,whisker_corrs,name_whisker_corrs,neuron_correlations_with_params]=get_whisking_in_air_modulation(wbins,min_t,DiscreteData,exclude,include,safety_margin,ppms,amp_thresh,Amp,Setp,Angle,Phase);
% %%
% min_t=100;%minimal time in bin
% ampbins=0:.05:1;
% abssetpbins=0:.05:1;
% absangbins=0:.05:1;
% setpbins=0:.05:1;
% angbins=0:.05:1;
% wbins={phasebinN;ampbins;abssetpbins;absangbins;setpbins;angbins};
% [Num_bin_occ,rates_to_plot,plot_x,name_vars,PhaseMod,whisker_corrs,name_whisker_corrs,neuron_correlations_with_params]=get_whisking_in_air_modulation_norm(wbins,min_t,DiscreteData,exclude,include,safety_margin,ppms,amp_thresh,Amp,Setp,Angle,Phase);

%%
cell_select={vpmr;pomr};
cell_select_names={'VPM';'POm'};
cell_colors={'r';'b'};
min_pointsN=4;
figure('Position',[148         198        1593         688])
tt= tiledlayout(3,6);
title(tt,'mean whisking in air correlations')
for v=1:numel(rates_to_plot)
    nexttile(v);
    for c=1:numel(cell_select)
        sel=sum(~isnan(rates_to_plot{v}(cell_select{c},:)),1)>min_pointsN;
        stA=std(rates_to_plot{v}(cell_select{c},sel),0,1,'omitnan')./sqrt(sum(~isnan(rates_to_plot{v}(cell_select{c},sel)),1));
        mA=mean(rates_to_plot{v}(cell_select{c},sel),1,'omitnan');
        SA=cat(2,mA-stA,fliplr(mA+stA));
        plotSA=cat(2,plot_x{v}(sel),fliplr(plot_x{v}(sel)));
        fill(plotSA,SA,cell_colors{c},'FaceAlpha',.25,'EdgeColor','none');hold on
        plot(plot_x{v}(sel),mA,cell_colors{c});hold on;box off
    end
    ylabel('mean Rate (Hz)');
    xlabel(name_vars{v});
    title(name_vars{v})
    xlim(plot_x{v}([1 end]));
    ms=prctile(rates_to_plot{v}(cat(1,cell_select{:}),:),99,'all');
    nexttile(v+numel(rates_to_plot));
    imagesc(plot_x{v},1:numel(cell_select{1}),rates_to_plot{v}(cell_select{1},:),[0 ms]);title([cell_select_names{1} ' ' name_vars{v}]);box off;colormap(CustomColormap);
    xlim(plot_x{v}([1 end]));
    colorbar;
    nexttile(v+2*numel(rates_to_plot));
    imagesc(plot_x{v},1:numel(cell_select{2}),rates_to_plot{v}(cell_select{2},:),[0 ms]);title([cell_select_names{2} ' ' name_vars{v}]);box off;colormap(CustomColormap);
    xlim(plot_x{v}([1 end]));
    colorbar;
end
%if write_figs;exportgraphics (gcf,[figdir Fig_no 'response_relative_to_whisking_in_air.pdf'], 'BackgroundColor','none','ContentType','vector');end
%%
figure
plot(plot_x{1},rates_to_plot{1}(vpmr,:),'Color',[.75 .75 .75]);hold on
%%
figure('Position',[148         198        1593         688])
tt= tiledlayout(2,6);
title(tt,'mean whisking in air correlations')
for v=1:numel(rates_to_plot)
    
    for c=1:numel(cell_select)
        nexttile(v+numel(rates_to_plot)*(c-1));
       % sel=sum(~isnan(rates_to_plot{v}(cell_select{c},:)),1)>min_pointsN;
       plot(plot_x{v},rates_to_plot{v}(cell_select{c},:),'Color',[.75 .75 .75]);hold on
       % stA=std(rates_to_plot{v}(cell_select{c},sel),0,1,'omitnan')./sqrt(sum(~isnan(rates_to_plot{v}(cell_select{c},sel)),1));
        mA=mean(rates_to_plot{v}(cell_select{c},:),1,'omitnan');
%         SA=cat(2,mA-stA,fliplr(mA+stA));
%         plotSA=cat(2,plot_x{v}(sel),fliplr(plot_x{v}(sel)));
%         fill(plotSA,SA,cell_colors{c},'FaceAlpha',.25,'EdgeColor','none');hold on
        plot(plot_x{v},mA,cell_colors{c});hold off;box off
         ylabel('mean Rate (Hz)');
    xlabel(name_vars{v});
    title(name_vars{v})
    xlim(plot_x{v}([1 end]));
    end
   
%     ms=prctile(rates_to_plot{v}(cat(1,cell_select{:}),:),99,'all');
%     nexttile(v+numel(rates_to_plot));
%     imagesc(plot_x{v},1:numel(cell_select{1}),rates_to_plot{v}(cell_select{1},:),[0 ms]);title([cell_select_names{1} ' ' name_vars{v}]);box off;colormap(CustomColormap);
%     xlim(plot_x{v}([1 end]));
%     colorbar;
%     nexttile(v+2*numel(rates_to_plot));
%     imagesc(plot_x{v},1:numel(cell_select{2}),rates_to_plot{v}(cell_select{2},:),[0 ms]);title([cell_select_names{2} ' ' name_vars{v}]);box off;colormap(CustomColormap);
%     xlim(plot_x{v}([1 end]));
%     colorbar;
end
if write_figs;exportgraphics (gcf,[figdir Fig_no 'response_relative_to_whisking_in_air.pdf'], 'BackgroundColor','none','ContentType','vector');end
%%
cell_select={vpmra;pomra};
figure('Position',[680   558   560   420],'Color','white')
tiledlayout(2,3);
nexttile(1,[2 2]);
polarscatter(PhaseMod.pref_phase(cell_select{1}),PhaseMod.SNR(cell_select{1}),75,'m.');hold on
polarscatter(PhaseMod.pref_phase(cell_select{1}(PhaseMod.p(cell_select{1})<.05)),PhaseMod.SNR(cell_select{1}(PhaseMod.p(cell_select{1})<.05)),75,'r.');hold on
polarscatter(PhaseMod.pref_phase(cell_select{2}),PhaseMod.SNR(cell_select{2}),75,'c.');hold on
polarscatter(PhaseMod.pref_phase(cell_select{2}(PhaseMod.p(cell_select{2})<.05)),PhaseMod.SNR(cell_select{2}(PhaseMod.p(cell_select{2})<.05)),75,'b.');hold off

title('phase SNR','preferred phase');
nexttile(3);
histogram(PhaseMod.SNR(cell_select{1}),0:.25:max(PhaseMod.SNR(cat(1,cell_select{1},cell_select{2})))+.2,'FaceColor','r');hold on;box off
title('phase SNR');
nexttile(6);
histogram(PhaseMod.SNR(cell_select{2}),0:.25:max(PhaseMod.SNR(cat(1,cell_select{1},cell_select{2})))+.2,'FaceColor','b');hold on;box off
if write_figs;exportgraphics (gcf,[figdir Fig_no 'phase_mod_whisking_in_air.pdf'], 'BackgroundColor','none','ContentType','vector');end

%%



%
%
figure('Position',[213         387        1587         401],'Color','white')
tiledlayout(2,numel(rates_to_plot)+2)
for v=1:numel(rates_to_plot)
    nexttile(v);
    histogram(neuron_correlations_with_params{v}(vpmr,1),-1:.2:1,'FaceColor','r');
    title(['VPM corr. with ' name_vars{v}]);box off;
    nexttile(v+numel(rates_to_plot)+2);
    histogram(neuron_correlations_with_params{v}(pomr,1),-1:.2:1,'FaceColor','b');
    title(['POm corr. with ' name_vars{v}]);box off;
end
nexttile(numel(rates_to_plot)+1,[2 2]);
rcmap=get_blue_red_cmap;
imagesc(1:numel(name_whisker_corrs),1:numel(name_whisker_corrs),whisker_corrs,[-1 1]);colormap(rcmap);
set(gca,'YTickLabel',name_whisker_corrs,'XTickLabel',name_whisker_corrs,'XTickLabelRotation',30,'YTickLabelRotation',30)
cb=colorbar;cb.Label.String='correlation coefficient';
title('Correlation between whisker variables')
if write_figs;exportgraphics (gcf,[figdir Fig_no 'correlation_whisking_in_air.pdf'], 'BackgroundColor','none','ContentType','vector');end


%% Figure 5
Fig_no='Fig5_';
parameter='Puff';
exclude={'Pole';'Exclude';'Grooming'};
include={{'Light'};[-5 35]};
contPlot={MeanWpL,MeanPL,MeanL};
figure('Position',[38 260 1413 525],'Color','White');Name='PuffLPop';
plot_Pop_PSTH_wrapper3(parameter, DiscreteData,puff_triglength_limit,ppms,psth_time_before_trig,psth_time_after_trig,exclude,psth_safety_margin,psth_binsize,[],Name,{vpm;pom},{'VPM';'POm'},include,'latency',100,contPlot,{'meanbefore','scaleminmax','scaleminmax'},'baseline corr. whisker angle');
exportgraphics(gcf,[figdir Fig_no 'Pop_all_' parameter '_L_instRateResp_' datestr(now,'yymmdd') '.pdf'],'BackgroundColor','none')
figure('Position',[38 260 1413 525],'Color','White');Name='PuffLPopAllR';
plot_Pop_PSTH_wrapper3(parameter, DiscreteData,puff_triglength_limit,ppms,psth_time_before_trig,psth_time_after_trig,exclude,psth_safety_margin,psth_binsize,[],Name,{vpmra;pomra},{'VPM';'POm'},include,'latency',100,contPlot,{'meanbefore','scaleminmax','scaleminmax'},'baseline corr. whisker angle');
exportgraphics(gcf,[figdir Fig_no 'Pop_puff_responsive_' parameter '_L_instRateResp_' datestr(now,'yymmdd') '.pdf'],'BackgroundColor','none')
vpmL=vpm(hasLight(vpm));
pomL=pom(hasLight(pom));
vpmrLa=vpmra(hasLight(vpmra));
pomrLa=pomra(hasLight(pomra));
contPlot={MeanWp,MeanP};
exclude={'Pole';'Light';'Exclude';'Grooming'};
include=[];
figure('Position',[38 260 1413 525],'Color','White');Name='PuffLPop';
plot_Pop_PSTH_wrapper3(parameter, DiscreteData,puff_triglength_limit,ppms,psth_time_before_trig,psth_time_after_trig,exclude,psth_safety_margin,psth_binsize,[],Name,{vpmL;pomL},{'VPM';'POm'},include,'latency',100,contPlot,{'meanbefore','scaleminmax','scaleminmax'},'baseline corr. whisker angle');
exportgraphics(gcf,[figdir Fig_no 'Pop_all_' parameter '_nL_instRateResp_' datestr(now,'yymmdd') '.pdf'],'BackgroundColor','none')
figure('Position',[38 260 1413 525],'Color','White');Name='PuffLPopAllR';
plot_Pop_PSTH_wrapper3(parameter, DiscreteData,puff_triglength_limit,ppms,psth_time_before_trig,psth_time_after_trig,exclude,psth_safety_margin,psth_binsize,[],Name,{vpmrLa;pomrLa},{'VPM';'POm'},include,'latency',100,contPlot,{'meanbefore','scaleminmax','scaleminmax'},'baseline corr. whisker angle');
exportgraphics(gcf,[figdir Fig_no 'Pop_puff_responsive_' parameter '_nL_instRateResp_' datestr(now,'yymmdd') '.pdf'],'BackgroundColor','none')

trigLen=cellfun(@(x) ones(1,numel(x))*30*20,trigsPL,'UniformOutput',0);
parameter=[trigsPL trigLen];
[H_PuffL,~,p_PuffL]=get_PSTH_pop(parameter, DiscreteData,puff_triglength_limit,ppms,rate_time_before_trig,rate_time_after_trig,exclude,psth_safety_margin,rate_binsize,isi,'left',[]);
p_PuffL=p_PuffL{1}(:,2);
sig_PuffL=p_PuffL<sig;
H_PL=cat(2,H_PuffL{:})';

figure('Position',[910   176   500   588])
plot_ratecomp2(H_PL,sig_PuffL,{'before PuffL';'after PuffL'},{vpm;pom},{'VPM';'POm'},'meansd','Rate [Hz]');title('Puff Light Response, all cells')
if write_figs;exportgraphics (gcf,[figdir Fig_no 'RatesPuffLAllCells.pdf'], 'BackgroundColor','none','ContentType','vector');end
figure('Position',[ 910   176   500   588])
plot_ratecomp2(H_PL,sig_PuffL,{'before PuffL';'after PuffL'},{vpmra;pomra},{'VPM';'POm'},'meansd','Rate [Hz]');title('Puff Light Response, responsive cells')
if write_figs;exportgraphics (gcf,[figdir Fig_no 'RatesPuffLRespondingCells.pdf'], 'BackgroundColor','none','ContentType','vector');end


figure('Position',[910   176   500   588])
plot_ratecomp2(H_P,sig_Puff,{'before PuffL';'after PuffL'},{vpmL;pomL},{'VPM';'POm'},'meansd','Rate [Hz]');title('Puff no Light Response, all cells')
if write_figs;exportgraphics (gcf,[figdir Fig_no 'RatesPuffnLAllCells.pdf'], 'BackgroundColor','none','ContentType','vector');end
figure('Position',[ 910   176   500   588])
plot_ratecomp2(H_P,sig_Puff,{'before PuffL';'after PuffL'},{vpmrLa;pomrLa},{'VPM';'POm'},'meansd','Rate [Hz]');title('Puff no Light Response, responsive cells')
if write_figs;exportgraphics (gcf,[figdir Fig_no 'RatesPuffnLRespondingCells.pdf'], 'BackgroundColor','none','ContentType','vector');end


%LnL puff comp
%significance: compare rates directly?
tempD=cat(2,diff(H_P,1,2),diff(H_PL,1,2));tempD(any(isnan(tempD),2),:)=nan;
figure('Position',[ 910   176   500   588])
plot_ratecomp2(tempD,sig_PuffW,{'noLPuff','LPuff'},{vpmra;pomra},{'VPM';'POm'},'meansd','deltaRate [Hz]');title('diffnL/L Puff Response, responsive cells')
if write_figs;exportgraphics (gcf,[figdir Fig_no 'diffRates_nLL_PuffRespondingCells.pdf'],'BackgroundColor','none','ContentType','vector');end
tempD=cat(2,H_P(:,1),H_PL(:,1));tempD(any(isnan(tempD),2),:)=nan;
figure('Position',[ 910   176   500   588])
plot_ratecomp2(tempD,sig_PuffW,{'noLPuff','LPuff'},{vpmra;pomra},{'VPM';'POm'},'meansd','Rate [Hz]');title('baseline activity before nL/L Puff, responsive cells')
if write_figs;exportgraphics (gcf,[figdir Fig_no 'baseRates_nLL_PuffRespondingCells.pdf'],'BackgroundColor','none','ContentType','vector');end
tempD=cat(2,H_P(:,2),H_PL(:,2));tempD(any(isnan(tempD),2),:)=nan;
figure('Position',[ 910   176   500   588])
plot_ratecomp2(tempD,sig_PuffW,{'noLPuff','LPuff'},{vpmra;pomra},{'VPM';'POm'},'meansd','Rate [Hz]');title('nL/L Puff Response, responsive cells')
if write_figs;exportgraphics (gcf,[figdir Fig_no 'absRates_nLL_PuffRespondingCells.pdf'], 'BackgroundColor','none','ContentType','vector');end
tempD=cat(2,H_P(:,2)./H_P(:,1),H_PL(:,2)./H_PL(:,1));tempD(any(isnan(tempD),2),:)=nan;
figure('Position',[ 910   176   500   588])
plot_ratecomp2(tempD,sig_Touch,{'noLPuff','LPuff'},{vpmra;pomra},{'VPM';'POm'},'meansd','Rate fold change');title('nL/L Puff fold change Response, responsive cells')
if write_figs;exportgraphics (gcf,[figdir Fig_no 'foldchangeRates_nLL_RespondingCells.pdf'], 'BackgroundColor','none','ContentType','vector');end

%%

exclude={'Pole';'Exclude';'Grooming';'Puff'};
include={'Light'};
safety_margin=[25 50];
min_t=1000;%minimal time in bin

phasebinN=16;
ampbins=0:3:65;
abssetpbins=0:3:40;
absangbins=0:3:65;
setpbins=-30:3:40;
angbins=-40:3:65;
wbins={phasebinN;ampbins;abssetpbins;absangbins;setpbins;angbins};
[Num_bin_occL,rates_to_plotL,plot_x,name_vars,PhaseModL,whisker_corrsL,name_whisker_corrs,neuron_correlations_with_paramsL]=get_whisking_in_air_modulation(wbins,min_t,DiscreteData,exclude,include,safety_margin,ppms,amp_thresh,Amp,Setp,Angle,Phase);



cell_select={vpmra;pomra};
cell_select_names={'VPM';'POm'};
cell_colors={'r';'b'};
min_pointsN=4;
figure('Position',[148         198        1593         688])
tt= tiledlayout(3,6);
title(tt,'mean whisking in air correlations')
for v=1:numel(rates_to_plotL)
    nexttile(v);
    for c=1:numel(cell_select)
        sel=sum(~isnan(rates_to_plotL{v}(cell_select{c},:)),1)>min_pointsN;
        stA=std(rates_to_plotL{v}(cell_select{c},sel),0,1,'omitnan')./sqrt(sum(~isnan(rates_to_plotL{v}(cell_select{c},sel)),1));
        mA=mean(rates_to_plotL{v}(cell_select{c},sel),1,'omitnan');
        SA=cat(2,mA-stA,fliplr(mA+stA));
        plotSA=cat(2,plot_x{v}(sel),fliplr(plot_x{v}(sel)));
        fill(plotSA,SA,cell_colors{c},'FaceAlpha',.25,'EdgeColor','none');hold on
        plot(plot_x{v}(sel),mA,cell_colors{c});hold on;box off
    end
    ylabel('mean Rate (Hz)');
    xlabel(name_vars{v});
    title(name_vars{v})
    xlim(plot_x{v}([1 end]));
    ms=prctile(rates_to_plotL{v}(cat(1,cell_select{:}),:),99,'all');
    nexttile(v+numel(rates_to_plotL));
    imagesc(plot_x{v},1:numel(cell_select{1}),rates_to_plotL{v}(cell_select{1},:),[0 ms]);title([cell_select_names{1} ' ' name_vars{v}]);box off;colormap(CustomColormap);
    xlim(plot_x{v}([1 end]));
    colorbar;
    nexttile(v+2*numel(rates_to_plotL));
    imagesc(plot_x{v},1:numel(cell_select{2}),rates_to_plotL{v}(cell_select{2},:),[0 ms]);title([cell_select_names{2} ' ' name_vars{v}]);box off;colormap(CustomColormap);
    xlim(plot_x{v}([1 end]));
    colorbar;
end
if write_figs;exportgraphics (gcf,[figdir Fig_no 'response_relative_to_whisking_in_air_light.pdf'], 'BackgroundColor','none','ContentType','vector');end

%
figure('Position',[680   558   560   420],'Color','white')
tiledlayout(2,3);
nexttile(1,[2 2]);
polarscatter(PhaseModL.pref_phase(vpmra),PhaseModL.MD(vpmra),'r.');hold on
polarscatter(PhaseModL.pref_phase(pomra),PhaseModL.MD(pomra),'b.');hold off
title('phase modulation depth','preferred phase');
nexttile(3);
histogram(PhaseModL.MD(vpmra),0:.2:max(PhaseModL.MD(cat(1,vpmra,pomra)))+.2,'FaceColor','r');hold on;box off
title('phase modulation depth');
nexttile(6);
histogram(PhaseModL.MD(pomra),0:.2:max(PhaseModL.MD(cat(1,vpmra,pomra)))+.2,'FaceColor','b');hold on;box off
if write_figs;exportgraphics (gcf,[figdir Fig_no 'phase_mod_whisking_in_air_light.pdf'], 'BackgroundColor','none','ContentType','vector');end

%



%
%
figure('Position',[213         387        1587         401],'Color','white')
tiledlayout(2,numel(rates_to_plot)+2)
for v=1:numel(rates_to_plot)
    nexttile(v);
    histogram(neuron_correlations_with_paramsL{v}(vpmra,1),-1:.1:1,'FaceColor','r');
    title(['VPM corr. with ' name_vars{v}]);box off;
    nexttile(v+numel(rates_to_plot)+2);
    histogram(neuron_correlations_with_paramsL{v}(pomra,1),-1:.1:1,'FaceColor','b');
    title(['POm corr. with ' name_vars{v}]);box off;
end
nexttile(numel(rates_to_plot)+1,[2 2]);
rcmap=get_blue_red_cmap;
imagesc(1:numel(name_whisker_corrs),1:numel(name_whisker_corrs),whisker_corrsL,[-1 1]);colormap(rcmap);
set(gca,'YTickLabel',name_whisker_corrs,'XTickLabel',name_whisker_corrs,'XTickLabelRotation',30,'YTickLabelRotation',30)
cb=colorbar;cb.Label.String='correlation coefficient';
title('Correlation between whisker variables')
if write_figs;exportgraphics (gcf,[figdir Fig_no 'correlation_whisking_in_air_light.pdf'], 'BackgroundColor','none','ContentType','vector');end


%% Figure 6
%% Figure 7



















