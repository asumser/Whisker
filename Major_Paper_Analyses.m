%% Puff/Touch Analyses
DrivePath='C:\Google Drive\';
DataPath='C:\Google Drive\Arbeit\Projects\Database\EphysData\';
load([DrivePath 'Arbeit\Projects\Database\ephys_database.mat'])
RecDB.Properties.RowNames=strcat(RecDB.AnimalName,{'_'},datestr(RecDB.StartTime,'HH_MM_SS'));
load([DataPath 'DiscreteData.mat'])
Wdir='C:\Google Drive\Arbeit\Projects\Database\EphysData\Whisker\';

ppms=20;
figdir=['C:\Data\Whisker\Paperfigs_matlab\' datestr(now,'YYmmDD') '\'];
set(0,'defaultfigurecolor',[1 1 1]);
set(0, 'DefaultAxesBox', 'off')
write_figs=true;
if write_figs
    mkdir(figdir)
end
%
minFillQ=4;
%response=RecDB.FastPuffResponse | RecDB.EarlyPuffResponse;
vpm=find(RecDB.DefinedNucleus=='VPM' & RecDB.FillQuality<=minFillQ & RecDB.RecordingType=='juxta');
pom=find(RecDB.DefinedNucleus=='POm' & RecDB.FillQuality<=minFillQ & RecDB.RecordingType=='juxta');
response=zeta_output.p_Puff<.05;
vpmr=find(RecDB.DefinedNucleus=='VPM' & RecDB.FillQuality<=minFillQ & RecDB.RecordingType=='juxta' & response);
pomr=find(RecDB.DefinedNucleus=='POm' & RecDB.FillQuality<=minFillQ & RecDB.RecordingType=='juxta' & response);
%%
load('F:\Whisker\Wdec\Wdec_by_field.mat', 'Angle','Phase','Amp')
%% Ratecomp Puff Diag
time_before_trig=35;%ms
time_after_trig=35;%ms
binsize=35;%ms
safety_margin=time_after_trig;
isi=10;
sig=.05;
%Puff
parameter='Puff';
exclude={'Pole';'Light';'Exclude';'Grooming'};
triglength_limit=[28 32];%ms
[H_Puff,~,p_Puff,~,ptrigs,lat]=get_PSTH_pop(parameter, DiscreteData,triglength_limit,ppms,time_before_trig,time_after_trig,exclude,safety_margin,binsize,isi,[]);
p_Puff=p_Puff{1}(:,2);
sig_Puff=p_Puff<sig;
%p=find_test_bin(parameter, DiscreteData,triglength_limit,ppms,time_before_trig,time_after_trig,exclude,safety_margin,[]);

H_Pc=cat(2,H_Puff{:})';
vpmr=vpm(sig_Puff(vpm));
pomr=pom(sig_Puff(pom));
%
for n=1:numel(lat)
    temp(n)=lat{n};
end
%
zeta_output.Latency_Puff=cat(1,temp(:).dblPeakT);
zeta_output.p_Puff=cat(1,temp(:).dblP);
zeta_output.Magnitude_Puff=cat(1,temp(:).dblZETA);
zeta_output.MagnitudeNeg_Puff=cat(1,temp(:).dblD_InvSign);
zeta_output.LatencyNeg_Puff=cat(1,temp(:).dblPeakT_InvSign);
figure,scatter(zeta_output.p_Puff,p_Puff)
%%
figure, scatter(zeta_output.Latency_Puff(vpmr),zeta_output.Magnitude_Puff(vpmr),15,'k','filled');hold on;
scatter(zeta_output.Latency_Puff(pomr),zeta_output.Magnitude_Puff(pomr),15,'r','filled');hold on;
legend('vpm','pom')
%%
dL=cell(numel(zeta_P),2);
for n=1:numel(zeta_P)
    dL{n,1}=zeta_P{1, 3}.vecD
end

%%
minFillQ=4;
vpm=find(RecDB.DefinedNucleus=='VPM' & RecDB.FillQuality<=minFillQ & RecDB.RecordingType=='juxta');
pom=find(RecDB.DefinedNucleus=='POm' & RecDB.FillQuality<=minFillQ & RecDB.RecordingType=='juxta');
vpmr=vpm(sig_Puff(vpm));
pomr=pom(sig_Puff(pom));
%% Ratecomp Puff
time_before_trig=35;%ms
time_after_trig=35;%ms
binsize=35;%ms
safety_margin=time_after_trig;
isi=10;
sig=.05;
%Puff
parameter='Puff';
exclude={'Pole';'Light';'Exclude';'Grooming'};
triglength_limit=[28 32];%ms
[H_Puff,~,p_Puff,~,ptrigs,lat]=get_PSTH_pop(parameter, DiscreteData,triglength_limit,ppms,time_before_trig,time_after_trig,exclude,safety_margin,binsize,isi,[]);
% p_Puff=p_Puff{1}(:,2);
% sig_Puff=p_Puff<sig;
H_Pc=cat(2,H_Puff{:})';
figure('Position',[910   176   500   588])
plot_ratecomp2(H_Pc,sig_Puff,{'before Puff';'after Puff'},{vpm;pom},{'VPM';'POm'},'meansd','Rate [Hz]');title('Puff Response, all cells')
if write_figs;export_fig ([figdir 'RatesPuffAllCells.pdf'], '-q95', '-pdf' ,'-transparent' ,'-native');end
figure('Position',[ 910   176   500   588])
plot_ratecomp2(H_Pc,sig_Puff,{'before Puff';'after Puff'},{vpmr;pomr},{'VPM';'POm'},'meansd','Rate [Hz]');title('Puff Response, responsive cells')
if write_figs;export_fig ([figdir 'RatesPuffRespondingCells.pdf'], '-q95', '-pdf' ,'-transparent' ,'-native');end

% figure('Position',[223         257        1013         642],'Name','Puff all VPM&POm cells')
% plot_ratecomp_diag(H_Pc,sig_Puff,{vpm,pom},'before Puff','after Puff',{'VPM','POm'},{'b','r'},'meansd')
% if write_figs;export_fig ([figdir 'RatesPuffDiagAllCells.pdf'], '-q95', '-pdf' ,'-transparent' ,'-native');end
% figure('Position',[223         257        1013         642],'Name','Puff Resp VPM&POm cells')
% plot_ratecomp_diag(H_Pc,sig_Puff,{vpmr,pomr},'before Puff','after Puff',{'VPM','POm'},{'b','r'},'meansd')
% if write_figs;export_fig ([figdir 'RatesPuffDiagRespCells.pdf'], '-q95', '-pdf' ,'-transparent' ,'-native');end

%
%% Ratecomp Touch
time_before_trig=35;%ms
time_after_trig=35;%ms
binsize=35;%ms
safety_margin=time_after_trig;
isi=10;
sig=.05;
%Puff
parameter='Touch';
exclude={'Puff';'Light';'Exclude';'Grooming'};
triglength_limit=[0 inf];%ms
[H_Touch,~,p_Touch,~,ttrigs,latT]=get_PSTH_pop(parameter, DiscreteData,triglength_limit,ppms,time_before_trig,time_after_trig,exclude,safety_margin,binsize,isi,[]);
p_Touch=p_Touch{1}(:,2);
sig_Touch=p_Touch<sig;
%p=find_test_bin(parameter, DiscreteData,triglength_limit,ppms,time_before_trig,time_after_trig,exclude,safety_margin,[]);

H_Tc=cat(2,H_Touch{:})';
%%
%
for n=1:numel(lat)
    tempT(n)=latT{n};
end
%
zeta_output.Latency_Touch=cat(1,tempT(:).dblPeakT);
zeta_output.p_Touch=cat(1,tempT(:).dblP);
zeta_output.Magnitude_Touch=cat(1,tempT(:).dblZETA);
zeta_output.MagnitudeNeg_Touch=cat(1,tempT(:).dblD_InvSign);
zeta_output.LatencyNeg_Touch=cat(1,tempT(:).dblPeakT_InvSign);
zeta_output.Magnitude_Touch(isnan(zeta_output.Latency_Touch))=nan;
%%
figure
scatter(zeta_output.Magnitude(vpm),zeta_output.MagnitudeT(vpm),15,'k','filled');axis equal;hold on
scatter(zeta_output.Magnitude(pom),zeta_output.MagnitudeT(pom),15,'r','filled');axis equal;hold on
line([0 8],[0 8],'Color','k')
line([0 8],[2 2],'Color','k','LineStyle','--')
line([2 2],[0 8],'Color','k','LineStyle','--')
xlabel('Puff zeta');ylabel('Touch zeta')
xlim([0 8]);ylim([0 8]);
box off;
legend('vpm','pom','equi line','significance cutoff');legend('boxoff')
%%
time_before_trig=35;%ms
time_after_trig=35;%ms
binsize=35;%ms
safety_margin=time_after_trig;
isi=10;
sig=.05;
%Puff
parameter='Touch';
exclude={'Puff';'Light';'Exclude';'Grooming'};
triglength_limit=[0 inf];%ms
[H_Touch,~,p_Touch,~,ttrigs,latT]=get_PSTH_pop(parameter, DiscreteData,triglength_limit,ppms,time_before_trig,time_after_trig,exclude,safety_margin,binsize,isi,[]);
% p_Touch=p_Touch{1}(:,2);
% sig_Touch=p_Touch<sig;
%p=find_test_bin(parameter, DiscreteData,triglength_limit,ppms,time_before_trig,time_after_trig,exclude,safety_margin,[]);

H_Tc=cat(2,H_Touch{:})';

figure('Position',[ 910   176   500   588])
plot_ratecomp2(H_Tc,sig_Touch,{'before Touch';'after Touch'},{vpm;pom},{'VPM';'POm'},'meansd','Rate [Hz]');title('Touch Response, all cells')
if write_figs;export_fig ([figdir 'RatesTouchAllCells.pdf'], '-q95', '-pdf' ,'-transparent' ,'-native');end
figure('Position',[ 910   176   500   588])
plot_ratecomp2(H_Tc,sig_Touch,{'before Touch';'after Touch'},{vpmr;pomr},{'VPM';'POm'},'meansd','Rate [Hz]');title('Touch Response, responsive cells')
if write_figs;export_fig ([figdir 'RatesTouchRespondingCells.pdf'], '-q95', '-pdf' ,'-transparent' ,'-native');end

% figure('Position',[223         257        1013         642],'Name','Touch all VPM&POm cells')
% plot_ratecomp_diag(H_Tc,sig_Touch,{vpm,pom},'before Touch','after Touch',{'VPM','POm'},{'b','r'},'meansd')
% if write_figs;export_fig ([figdir 'RatesTouchDiagAllCells.pdf'], '-q95', '-pdf' ,'-transparent' ,'-native');end
% 
% figure('Position',[223         257        1013         642],'Name','Touch PuffResp VPM&POm cells')
% plot_ratecomp_diag(H_Tc,sig_Touch,{vpmr,pomr},'before Touch','after Touch',{'VPM','POm'},{'b','r'},'meansd')
% if write_figs;export_fig ([figdir 'RatesTouchDiagPuffRespCells.pdf'], '-q95', '-pdf' ,'-transparent' ,'-native');end


%% Ratecomp Puff Touch

H_PT=cat(2,diff(H_Pc,1,2),diff(H_Tc,1,2));
H_PT=cat(2,zeta_output.Magnitude_Puff,zeta_output.Magnitude_Touch);
% H_PT=cat(2,log2(H_Pc(:,2)./H_Pc(:,1)),log2(H_Tc(:,2)./H_Tc(:,1)));
% figure('Position',[223         257        1013         642],'Name','Touch_Puff all VPM&POm cells')
% plot_ratecomp_diag(H_PT,any(cat(2,sig_Touch,sig_Puff),2),{vpm,pom},'zeta Puff','zeta Touch',{'VPM','POm'},{'b','r'},'meansd')
% if write_figs;export_fig ([figdir 'RatesPuffTouchDiagAllCells.pdf'], '-q95', '-pdf' ,'-transparent' ,'-native');end
% 
% figure('Position',[223         257        1013         642],'Name','Touch PuffResp VPM&POm cells')
% plot_ratecomp_diag(H_PT,any(cat(2,sig_Touch,sig_Puff),2),{vpmr,pomr},'zeta Puff','zeta Touch',{'VPM','POm'},{'b','r'},'meansd')
% if write_figs;export_fig ([figdir 'RatesPuffTouchDiagPuffRespCells.pdf'], '-q95', '-pdf' ,'-transparent' ,'-native');end

%
figure('Position',[ 910   176   500   588])
plot_ratecomp2(H_PT,any(cat(2,sig_Touch,sig_Puff),2),{'Puff','Touch'},{vpm;pom},{'VPM';'POm'},'meansd','zeta');title('Touch/Puff Response, all cells')
if write_figs;export_fig ([figdir 'RatesTouchVsPuffAllCells.pdf'], '-q95', '-pdf' ,'-transparent' ,'-native');end
figure('Position',[ 910   176   500   588])
plot_ratecomp2(H_PT,any(cat(2,sig_Touch,sig_Puff),2),{'Puff','Touch'},{vpmr;pomr},{'VPM';'POm'},'meansd','zeta');title('Touch/Puff Response, responsive cells')
if write_figs;export_fig ([figdir 'RatesTouchVsPuffRespondingCells.pdf'], '-q95', '-pdf' ,'-transparent' ,'-native');end


%% Ratecomp Whisk Diag
time_before_trig=300;%ms
time_after_trig=300;%ms
binsize=300;%ms
safety_margin=time_after_trig;
isi=10;
sig=.05;
%Puff
parameter='Whisking';
exclude={'Puff';'Pole''Light';'Exclude';'Grooming'};
triglength_limit=[time_after_trig inf];%ms
[H_W,~,p_W]=get_PSTH_pop(parameter, DiscreteData,triglength_limit,ppms,time_before_trig,time_after_trig,exclude,safety_margin,binsize,isi,[]);
p_W=p_W{1}(:,2);
sig_W=p_W<sig;
%p=find_test_bin(parameter, DiscreteData,triglength_limit,ppms,time_before_trig,time_after_trig,exclude,safety_margin,[]);

H_Wc=cat(2,H_W{:})';

figure('Position',[223         257        1013         642],'Name','Whisk all VPM&POm cells')
plot_ratecomp_diag(H_Wc,sig_W,{vpm,pom},'before Whisking','while Whisking',{'VPM','POm'},{'b','r'},'meansd')
export_fig ([figdir 'RatesWhiskDiagAllCells.pdf'], '-q95', '-pdf' ,'-transparent' ,'-native')

figure('Position',[223         257        1013         642],'Name','Whisk PuffResp VPM&POm cells')
plot_ratecomp_diag(H_Wc,sig_W,{vpmr,pomr},'before Whisking','while Whisking',{'VPM','POm'},{'b','r'},'meansd')
export_fig ([figdir 'RatesWhiskDiagPuffRespCells.pdf'], '-q95', '-pdf' ,'-transparent' ,'-native')

%
figure('Position',[ 910   176   500   588])
plot_ratecomp2(H_Wc,sig_W,{'preWhisking','Whisking'},{vpm;pom},{'VPM';'POm'},'meansd','Rate [Hz]');title('Whisking Response, all cells')
if write_figs;export_fig ([figdir 'RatesWhiskingAllCells.pdf'], '-q95', '-pdf' ,'-transparent' ,'-native');end
figure('Position',[ 910   176   500   588])
plot_ratecomp2(H_Wc,sig_W,{'preWhisking','Whisking'},{vpmr;pomr},{'VPM';'POm'},'meansd','Rate [Hz]');title('Whisking Response, responsive cells')
if write_figs;export_fig ([figdir 'RatesWhiskingRespondingCells.pdf'], '-q95', '-pdf' ,'-transparent' ,'-native');end

%% Puff when (non)Whisking
%general params
time_before_trig=35;%ms
time_after_trig=35;%ms
binsize=35;%ms
safety_margin=time_after_trig;
isi=10;
%Puff
include_margin=[50 0];
triglength_limit=[28 32];%ms
parameter='Puff';

exclude={'Pole';'Light';'Exclude';'Grooming'};
include_forced={'Whisking'};
[H_PuffW,~,p_PuffW,~,p_trigs_W,zeta_W]=get_PSTH_pop(parameter, DiscreteData,triglength_limit,ppms,time_before_trig,time_after_trig,exclude,safety_margin,binsize,isi,[],include_forced,include_margin);
%
H_PW=cat(2,H_PuffW{:})';
p_PuffW=p_PuffW{1}(:,2);
sig_PuffW=p_PuffW<sig;

figure('Position',[910   176   500   588])
plot_ratecomp2(H_PW,sig_PuffW,{'before PuffW';'after PuffW'},{vpm;pom},{'VPM';'POm'},'meansd','Rate [Hz]');title('Puff Response Whisking, all cells')
if write_figs;export_fig ([figdir 'RatesPuffWhiskingAllCells.pdf'], '-q95', '-pdf' ,'-transparent' ,'-native');end
figure('Position',[ 910   176   500   588])
plot_ratecomp2(H_PW,sig_PuffW,{'before PuffW';'after PuffW'},{vpmr;pomr},{'VPM';'POm'},'meansd','Rate [Hz]');title('Puff Response Whisking, responsive cells')
if write_figs;export_fig ([figdir 'RatesPuffWhiskingRespondingCells.pdf'], '-q95', '-pdf' ,'-transparent' ,'-native');end
% 
% figure('Position',[223         257        1013         642],'Name','Puff all VPM&POm cells')
% plot_ratecomp_diag(H_PW,sig_PuffW,{vpm,pom},'before PuffW','after PuffW',{'VPM','POm'},{'b','r'},'meansd')
% if write_figs;export_fig ([figdir 'RatesPuffDiagWhiskingDiagAllCells.pdf'], '-q95', '-pdf' ,'-transparent' ,'-native');end
% figure('Position',[223         257        1013         642],'Name','Puff Resp VPM&POm cells')
% plot_ratecomp_diag(H_PW,sig_PuffW,{vpmr,pomr},'before PuffW','after PuffW',{'VPM','POm'},{'b','r'},'meansd')
% if write_figs;export_fig ([figdir 'RatesPuffDiagWhiskingDiagRespCells.pdf'], '-q95', '-pdf' ,'-transparent' ,'-native');end
% %

safety_margin=[time_after_trig time_after_trig];
exclude={'Pole';'Light';'Exclude';'Grooming';'Whisking'};
[H_PuffnW,~,p_PuffnW,~,p_trigs_nW,zeta_nW]=get_PSTH_pop(parameter, DiscreteData,triglength_limit,ppms,time_before_trig,time_after_trig,exclude,safety_margin,binsize,isi,[]);

H_PnW=cat(2,H_PuffnW{:})';
p_PuffnW=p_PuffnW{1}(:,2);
sig_PuffnW=p_PuffnW<sig;

figure('Position',[910   176   500   588])
plot_ratecomp2(H_PnW,sig_PuffnW,{'before PuffQ';'after PuffQ'},{vpm;pom},{'VPM';'POm'},'meansd','Rate [Hz]');title('Puff Response nonWhisking, all cells')
if write_figs;export_fig ([figdir 'RatesPuffQAllCells.pdf'], '-q95', '-pdf' ,'-transparent' ,'-native');end
figure('Position',[ 910   176   500   588])
plot_ratecomp2(H_PnW,sig_PuffnW,{'before PuffQ';'after PuffQ'},{vpmr;pomr},{'VPM';'POm'},'meansd','Rate [Hz]');title('Puff Response nonWhisking, responsive cells')
if write_figs;export_fig ([figdir 'RatesPuffQRespondingCells.pdf'], '-q95', '-pdf' ,'-transparent' ,'-native');end

% figure('Position',[223         257        1013         642],'Name','Puff all VPM&POm cells')
% plot_ratecomp_diag(H_PnW,sig_PuffnW,{vpm,pom},'before PuffQ','after PuffQ',{'VPM','POm'},{'b','r'},'meansd')
% if write_figs;export_fig ([figdir 'RatesPuffDiagQDiagAllCells.pdf'], '-q95', '-pdf' ,'-transparent' ,'-native');end
% figure('Position',[223         257        1013         642],'Name','Puff Resp VPM&POm cells')
% plot_ratecomp_diag(H_PnW,sig_PuffnW,{vpmr,pomr},'before PuffQ','after PuffQ',{'VPM','POm'},{'b','r'},'meansd')
% if write_figs;export_fig ([figdir 'RatesPuffDiagQDiagRespCells.pdf'], '-q95', '-pdf' ,'-transparent' ,'-native');end
% %
%%
H_PWQ=zeros(numel(zeta_W),2);
for n=1:numel(zeta_W)
    H_PWQ(n,1)=zeta_nW{n}.dblZETA;
    H_PWQ(n,2)=zeta_W{n}.dblZETA;
end


%H_PWQ=cat(2,diff(H_PnW,1,2),diff(H_PW,1,2));



% H_PT=cat(2,log2(H_Pc(:,2)./H_Pc(:,1)),log2(H_Tc(:,2)./H_Tc(:,1)));
% figure('Position',[223         257        1013         642],'Name','PuffQW all VPM&POm cells')
% plot_ratecomp_diag(H_PWQ,any(cat(2,sig_PuffnW,sig_PuffW),2),{vpm,pom},'diff QPuff','diff WPuff',{'VPM','POm'},{'b','r'},'meansd')
% if write_figs;export_fig ([figdir 'RatesPuffQWDiagAllCells.pdf'], '-q95', '-pdf' ,'-transparent' ,'-native');end
% 
% figure('Position',[223         257        1013         642],'Name','PuffQW resp VPM&POm cells')
% plot_ratecomp_diag(H_PWQ,any(cat(2,sig_PuffnW,sig_PuffW),2),{vpmr,pomr},'diff QPuff','diff WPuff',{'VPM','POm'},{'b','r'},'meansd')
% if write_figs;export_fig ([figdir 'RatesPuffQWDiagPuffRespCells.pdf'], '-q95', '-pdf' ,'-transparent' ,'-native');end
% 
% %
figure('Position',[ 910   176   500   588])
plot_ratecomp2(H_PWQ,any(cat(2,sig_PuffnW,sig_PuffW),2),{'diff QPuff','diff WPuff'},{vpm;pom},{'VPM';'POm'},'meansd','Rate change [dHz]');title('Puff Q/W Response, all cells')
if write_figs;export_fig ([figdir 'RatesPuffQWAllCells.pdf'], '-q95', '-pdf' ,'-transparent' ,'-native');end
figure('Position',[ 910   176   500   588])
plot_ratecomp2(H_PWQ,any(cat(2,sig_PuffnW,sig_PuffW),2),{'diff QPuff','diff WPuff'},{vpmr;pomr},{'VPM';'POm'},'meansd','Rate change [dHz]');title('Puff Q/W Response, responsive cells')
if write_figs;export_fig ([figdir 'RatesPuffQWRespondingCells.pdf'], '-q95', '-pdf' ,'-transparent' ,'-native');end





%% Ratecomp Puff Light
time_before_trig=35;%ms
time_after_trig=35;%ms
binsize=35;%ms
safety_margin=50;
isi=10;
sig=.05;
%Puff
parameter='Puff';
exclude={'Pole';'Light';'Exclude';'Grooming'};
triglength_limit=[28 32];%ms
[H_Puff,~,p_Puff]=get_PSTH_pop(parameter, DiscreteData,triglength_limit,ppms,time_before_trig,time_after_trig,exclude,safety_margin,binsize,isi,[]);
p_Puff=p_Puff{1}(:,2);
sig_Puff=p_Puff<sig;
H_Pc=cat(2,H_Puff{:})';

parameter='Puff';
exclude={'Pole';'Exclude';'Grooming'};
include={'Light'};
include_safety_margin=[50 50];
triglength_limit=[28 32];%ms
[H_PuffL,~,p_PuffL]=get_PSTH_pop(parameter, DiscreteData,triglength_limit,ppms,time_before_trig,time_after_trig,exclude,safety_margin,binsize,isi,[],include,include_safety_margin);
p_PuffL=p_PuffL{1}(:,2);
sig_PuffL=p_PuffL<sig;
H_PcL=cat(2,H_PuffL{:})';
%%

figure('Position',[910   176   500   588])
plot_ratecomp2(H_PcL,sig_PuffL,{'before PuffL';'after PuffL'},{vpm;pom},{'VPM';'POm'},'meansd','Rate [Hz]');title('Puff Light Response, all cells')
if write_figs;export_fig ([figdir 'RatesPuffLAllCells.pdf'], '-q95', '-pdf' ,'-transparent' ,'-native');end
figure('Position',[ 910   176   500   588])
plot_ratecomp2(H_PcL,sig_PuffL,{'before PuffL';'after PuffL'},{vpmr;pomr},{'VPM';'POm'},'meansd','Rate [Hz]');title('Puff Light Response, responsive cells')
if write_figs;export_fig ([figdir 'RatesPuffLRespondingCells.pdf'], '-q95', '-pdf' ,'-transparent' ,'-native');end
% 
% figure('Position',[223         257        1013         642],'Name','Puff all VPM&POm cells')
% plot_ratecomp_diag(H_PcL,sig_PuffL,{vpm,pom},'before PuffL','after PuffL',{'VPM','POm'},{'b','r'},'meansd')
% if write_figs;export_fig ([figdir 'RatesPuffLDiagAllCells.pdf'], '-q95', '-pdf' ,'-transparent' ,'-native');end
% figure('Position',[223         257        1013         642],'Name','Puff Resp VPM&POm cells')
% plot_ratecomp_diag(H_PcL,sig_PuffL,{vpmr,pomr},'before PuffL','after PuffL',{'VPM','POm'},{'b','r'},'meansd')
% if write_figs;export_fig ([figdir 'RatesPuffLDiagRespCells.pdf'], '-q95', '-pdf' ,'-transparent' ,'-native');end
%% diff
H_PnLL=cat(2,diff(H_Pc,1,2),diff(H_PcL,1,2));
% H_PT=cat(2,log2(H_Pc(:,2)./H_Pc(:,1)),log2(H_Tc(:,2)./H_Tc(:,1)));
figure('Position',[223         257        1013         642],'Name','PuffnLL all VPM&POm cells')
plot_ratecomp_diag(H_PnLL,any(cat(2,sig_Puff,sig_PuffL),2),{vpm,pom},'diff nLPuff','diff LPuff',{'VPM','POm'},{'b','r'},'meansd')
if write_figs;export_fig ([figdir 'RatesPuffnLLDiagAllCells.pdf'], '-q95', '-pdf' ,'-transparent' ,'-native');end

figure('Position',[223         257        1013         642],'Name','PuffnLL resp VPM&POm cells')
plot_ratecomp_diag(H_PnLL,any(cat(2,sig_Puff,sig_PuffL),2),{vpmr,pomr},'diff nLPuff','diff LPuff',{'VPM','POm'},{'b','r'},'meansd')
if write_figs;export_fig ([figdir 'RatesPuffnLLDiagPuffRespCells.pdf'], '-q95', '-pdf' ,'-transparent' ,'-native');end

%
figure('Position',[ 910   176   500   588])
plot_ratecomp2(H_PnLL,any(cat(2,sig_Puff,sig_PuffL),2),{'diff nLPuff','diff LPuff'},{vpm;pom},{'VPM';'POm'},'meansd','Rate change [dHz]');title('Puff nL/L Response, all cells')
if write_figs;export_fig ([figdir 'RatesPuffnLLAllCells.pdf'], '-q95', '-pdf' ,'-transparent' ,'-native');end
figure('Position',[ 910   176   500   588])
plot_ratecomp2(H_PnLL,any(cat(2,sig_Puff,sig_PuffL),2),{'diff nLPuff','diff LPuff'},{vpmr;pomr},{'VPM';'POm'},'meansd','Rate change [dHz]');title('Puff nL/L Response, responsive cells')
if write_figs;export_fig ([figdir 'RatesPuffnLLRespondingCells.pdf'], '-q95', '-pdf' ,'-transparent' ,'-native');end

%% Ratecomp Touch Light
time_before_trig=35;%ms
time_after_trig=35;%ms
binsize=35;%ms
safety_margin=time_after_trig;
isi=10;
sig=.05;
%Puff
parameter='Touch';
exclude={'Puff';'Light';'Exclude';'Grooming'};
triglength_limit=[28 32];%ms
[H_Touch,~,p_Touch]=get_PSTH_pop(parameter, DiscreteData,triglength_limit,ppms,time_before_trig,time_after_trig,exclude,safety_margin,binsize,isi,[]);
p_Touch=p_Touch{1}(:,2);
sig_Touch=p_Touch<sig;
H_Tc=cat(2,H_Touch{:})';

parameter='Touch';
exclude={'Puff';'Exclude';'Grooming'};
include={'Light'};
include_safety_margin=[35 35];
triglength_limit=[28 32];%ms
[H_TouchL,~,p_TouchL]=get_PSTH_pop(parameter, DiscreteData,triglength_limit,ppms,time_before_trig,time_after_trig,exclude,safety_margin,binsize,isi,[],include,include_safety_margin);
p_TouchL=p_TouchL{1}(:,2);
sig_TouchL=p_TouchL<sig;
H_TcL=cat(2,H_TouchL{:})';
%%

figure('Position',[910   176   500   588])
plot_ratecomp2(H_TcL,sig_TouchL,{'before TouchL';'after TouchL'},{vpm;pom},{'VPM';'POm'},'meansd','Rate [Hz]');title('Touch Light Response, all cells')
if write_figs;export_fig ([figdir 'RatesTouchLAllCells.pdf'], '-q95', '-pdf' ,'-transparent' ,'-native');end
figure('Position',[ 910   176   500   588])
plot_ratecomp2(H_TcL,sig_TouchL,{'before TouchL';'after TouchL'},{vpmr;pomr},{'VPM';'POm'},'meansd','Rate [Hz]');title('Touch Light Response, responsive cells')
if write_figs;export_fig ([figdir 'RatesTouchLRespondingCells.pdf'], '-q95', '-pdf' ,'-transparent' ,'-native');end

figure('Position',[223         257        1013         642],'Name','Touch all VPM&POm cells')
plot_ratecomp_diag(H_TcL,sig_TouchL,{vpm,pom},'before TouchL','after TouchL',{'VPM','POm'},{'b','r'},'meansd')
if write_figs;export_fig ([figdir 'RatesTouchLDiagAllCells.pdf'], '-q95', '-pdf' ,'-transparent' ,'-native');end
figure('Position',[223         257        1013         642],'Name','Touch Resp VPM&POm cells')
plot_ratecomp_diag(H_TcL,sig_TouchL,{vpmr,pomr},'before TouchL','after TouchL',{'VPM','POm'},{'b','r'},'meansd')
if write_figs;export_fig ([figdir 'RatesTouchLDiagRespCells.pdf'], '-q95', '-pdf' ,'-transparent' ,'-native');end
%% diff
H_PnLL=cat(2,diff(H_Pc,1,2),diff(H_PcL,1,2));
% H_PT=cat(2,log2(H_Pc(:,2)./H_Pc(:,1)),log2(H_Tc(:,2)./H_Tc(:,1)));
figure('Position',[223         257        1013         642],'Name','PuffnLL all VPM&POm cells')
plot_ratecomp_diag(H_PnLL,any(cat(2,sig_Puff,sig_PuffL),2),{vpm,pom},'diff nLPuff','diff LPuff',{'VPM','POm'},{'b','r'},'meansd')
if write_figs;export_fig ([figdir 'RatesPuffnLLDiagAllCells.pdf'], '-q95', '-pdf' ,'-transparent' ,'-native');end

figure('Position',[223         257        1013         642],'Name','PuffnLL resp VPM&POm cells')
plot_ratecomp_diag(H_PnLL,any(cat(2,sig_Puff,sig_PuffL),2),{vpmr,pomr},'diff nLPuff','diff LPuff',{'VPM','POm'},{'b','r'},'meansd')
if write_figs;export_fig ([figdir 'RatesPuffnLLDiagPuffRespCells.pdf'], '-q95', '-pdf' ,'-transparent' ,'-native');end

%
figure('Position',[ 910   176   500   588])
plot_ratecomp2(H_PnLL,any(cat(2,sig_Puff,sig_PuffL),2),{'diff nLPuff','diff LPuff'},{vpm;pom},{'VPM';'POm'},'meansd','Rate change [dHz]');title('Puff nL/L Response, all cells')
if write_figs;export_fig ([figdir 'RatesPuffnLLAllCells.pdf'], '-q95', '-pdf' ,'-transparent' ,'-native');end
figure('Position',[ 910   176   500   588])
plot_ratecomp2(H_PnLL,any(cat(2,sig_Puff,sig_PuffL),2),{'diff nLPuff','diff LPuff'},{vpmr;pomr},{'VPM';'POm'},'meansd','Rate change [dHz]');title('Puff nL/L Response, responsive cells')
if write_figs;export_fig ([figdir 'RatesPuffnLLRespondingCells.pdf'], '-q95', '-pdf' ,'-transparent' ,'-native');end



%% single cell examples
vpm_example=vpmr(1);
pom_example=pomr(6);

time_before_trig=50;%ms
time_after_trig=150;%ms
binsize=2.5;%ms
safety_margin=[time_before_trig+15 time_after_trig];
isi=10;
%Puff
parameter='Puff';
exclude={'Pole';'Light';'Exclude';'Grooming'};
triglength_limit=[28 32];%ms
[RasterHP,P_trig]=get_raster_pop(parameter, DiscreteData,triglength_limit,ppms,time_before_trig,time_after_trig,exclude,safety_margin,binsize,isi,[]);
[H_Puff]=get_PSTH_pop(parameter, DiscreteData,triglength_limit,ppms,time_before_trig,time_after_trig,exclude,safety_margin,binsize,isi,[]);

parameter='Touch';
exclude={'Puff';'Light';'Exclude';'Grooming'};
triglength_limit=[-inf inf];%ms
[RasterHT,T_trig]=get_raster_pop(parameter, DiscreteData,triglength_limit,ppms,time_before_trig,time_after_trig,exclude,safety_margin,binsize,isi,[]);
[H_Touch]=get_PSTH_pop(parameter, DiscreteData,triglength_limit,ppms,time_before_trig,time_after_trig,exclude,safety_margin,binsize,isi,[]);

%%
bins=-time_before_trig:binsize:time_after_trig;
pbins=bins(1:end-1)+diff(bins);

[WP_vpm_example]=get_trig_whisker(vpm_example,P_trig,pbins,Wdir,RecDB,DiscreteData);
[WT_vpm_example]=get_trig_whisker(vpm_example,T_trig,pbins,Wdir,RecDB,DiscreteData);


[WP_pom_example]=get_trig_whisker(pom_example,P_trig,pbins,Wdir,RecDB,DiscreteData);
[WT_pom_example]=get_trig_whisker(pom_example,T_trig,pbins,Wdir,RecDB,DiscreteData);

RastN=50;
figure('Position',[273   110   880   762])
subplot(4,2,1);
n=vpm_example;
ax(1)=scatter(RasterHP{n}(RasterHP{n}(:,2)<=RastN,1),RasterHP{n}(RasterHP{n}(:,2)<=RastN,2),'k.');box off;
xlim([-time_before_trig time_after_trig]);ylim([1 50]);
ylabel('trials')
title(['VPM Puff ' num2str(vpm_example)])

subplot(4,2,2);
n=pom_example;
ax(2)=scatter(RasterHP{n}(RasterHP{n}(:,2)<=RastN,1),RasterHP{n}(RasterHP{n}(:,2)<=RastN,2),'k.');box off;
xlim([-time_before_trig time_after_trig]);ylim([1 50]);
ylabel('trials')
title(['POm Puff ' num2str(pom_example)])


subplot(4,2,3);
n=vpm_example; 
ax(3)=bar(pbins,H_Puff{n},1,'FaceColor',ones(1,3)/2);box off;
xlim([-time_before_trig time_after_trig]);ylim([0 400]);
ylabel('Rate (Hz)')
xlabel('Time relative to Puff (ms)')
yyaxis right;
plot(pbins,WP_vpm_example,'r');box off;
ylabel('avg. whisker (°)')


subplot(4,2,4);
n=pom_example;
ax(4)=bar(pbins,H_Puff{n},1,'FaceColor',ones(1,3)/2);box off;
xlim([-time_before_trig time_after_trig]);ylim([0 400]);
ylabel('Rate (Hz)')
xlabel('Time relative to Puff (ms)')
yyaxis right;
plot(pbins,WP_pom_example,'r');box off;
ylabel('avg. whisker (°)')

subplot(4,2,5);
n=vpm_example;
ax(1)=scatter(RasterHT{n}(RasterHT{n}(:,2)<=RastN,1),RasterHT{n}(RasterHT{n}(:,2)<=RastN,2),'k.');box off;
xlim([-time_before_trig time_after_trig]);ylim([1 50]);
ylabel('trials')
title('VPM Touch')


subplot(4,2,6);
n=pom_example;
ax(2)=scatter(RasterHT{n}(RasterHT{n}(:,2)<=RastN,1),RasterHT{n}(RasterHT{n}(:,2)<=RastN,2),'k.');box off;
xlim([-time_before_trig time_after_trig]);ylim([1 50]);
ylabel('trials')
title('POm Touch')

subplot(4,2,7);
n=vpm_example; 
ax(3)=bar(pbins,H_Touch{n},1,'FaceColor',ones(1,3)/2);box off;
xlim([-time_before_trig time_after_trig]);ylim([0 400]);
ylabel('Rate (Hz)')
xlabel('Time relative to Touch (ms)')
yyaxis right;
plot(pbins,WT_vpm_example,'r');box off;
ylabel('avg. whisker (°)')

subplot(4,2,8);
n=pom_example;
ax(4)=bar(pbins,H_Touch{n},1,'FaceColor',ones(1,3)/2);box off;
xlim([-time_before_trig time_after_trig]);ylim([0 400]);
ylabel('Rate (Hz)')
xlabel('Time relative to Touch (ms)')
yyaxis right;
plot(pbins,WT_pom_example,'r');box off;
ylabel('avg. whisker (°)')

%export_fig ([figdir 'Fig_1_example_VPMPom_PuffTOuch' datestr(now,'yymmdd') '.pdf'], '-q95', '-pdf' ,'-transparent' ,'-native')
exportgraphics(gcf,[figdir 'Fig_1_example_VPMPom_PuffTOuch' datestr(now,'yymmdd') '.pdf'],'BackgroundColor','none','ContentType','vector')


%% Ratecomp general

outsafety=200;%ms
insafety=100;
isi=10;
include={'Whisking'};
exclude={'Puff';'Touch';'Light';'Exclude';'Grooming'};
[RatesW,BurstsW]=get_DD_Rates(DiscreteData,exclude,include, outsafety,insafety,ppms,isi);


outsafety=200;%ms
insafety=100;
isi=10;
include={};
exclude={'Whisking';'Puff';'Touch';'Light';'Exclude';'Grooming'};
[RatesQ,BurstsQ]=get_DD_Rates(DiscreteData,exclude,include, outsafety,insafety,ppms,isi);

outsafety=200;%ms
insafety=25;
isi=10;
include={'Whisking';'Light'};
exclude={'Puff';'Touch';'Exclude';'Grooming'};
[RatesWL,BurstsWL]=get_DD_Rates(DiscreteData,exclude,include, outsafety,insafety,ppms,isi);

outsafety=200;%ms
insafety=25;
isi=10;
include={'Light'};
exclude={'Whisking';'Puff';'Touch';'Exclude';'Grooming'};
[RatesQL,BurstsQL]=get_DD_Rates(DiscreteData,exclude,include, outsafety,insafety,ppms,isi);

outsafety=200;%ms
insafety=25;
isi=10;
include={'Light'};
exclude={'Puff';'Touch';'Exclude';'Grooming'};
[RatesL,BurstsL]=get_DD_Rates(DiscreteData,exclude,include, outsafety,insafety,ppms,isi);

outsafety=100;%ms
insafety=0;
isi=10;
include={'Puff'};
exclude={'Light';'Touch';'Exclude';'Grooming'};
[RatesP,BurstsP]=get_DD_Rates(DiscreteData,exclude,include, outsafety,insafety,ppms,isi);

outsafety=100;%ms
insafety=0;
isi=10;
include={'Light';'Puff'};
exclude={'Touch';'Exclude';'Grooming'};
[RatesLP,BurstsLP]=get_DD_Rates(DiscreteData,exclude,include, outsafety,insafety,ppms,isi);

outsafety=100;%ms
insafety=0;
isi=10;
include={'Touch'};
exclude={'Light';'Puff';'Exclude';'Grooming'};
[RatesT,BurstsT]=get_DD_Rates(DiscreteData,exclude,include, outsafety,insafety,ppms,isi);

outsafety=100;%ms
insafety=0;
isi=10;
include={'Touch';'Light'};
exclude={'Puff';'Exclude';'Grooming'};
[RatesLT,BurstsLT]=get_DD_Rates(DiscreteData,exclude,include, outsafety,insafety,ppms,isi);
%%
figpos=[ 910   176   500   588];sep_names={'VPM';'POm'};sel_all={vpm;pom};sel_resp={vpmr;pomr};sigV=true(size(RatesQ));
%QW
comp_label={'Q','W'};
values=[RatesQ RatesW];ylab='Average Rate [Hz]';figname='RatesGeneralQW';tit='Quiet & Whisking Rates';
ratecomp_wrapper(figpos,figdir,figname,write_figs,values,sigV,comp_label,sel_all,sel_resp,sep_names,ylab,tit)
values=[BurstsQ BurstsW];ylab='Spikes/Burst';figname='BurstsGeneralQW';tit='Quiet & Whisking Bursts';
ratecomp_wrapper(figpos,figdir,figname,write_figs,values,sigV,comp_label,sel_all,sel_resp,sep_names,ylab,tit)
%% QL
comp_label={'QnL','QL'};
values=[RatesQ RatesQL];ylab='Average Rate [Hz]';figname='RatesGeneralQL';tit='Quiet & Quiet+Light Rates';
ratecomp_wrapper(figpos,figdir,figname,write_figs,values,sigV,comp_label,sel_all,sel_resp,sep_names,ylab,tit)
values=[BurstsQ BurstsQL];ylab='Spikes/Burst';figname='BurstsGeneralQL';tit='Quiet & Quiet+Light Bursts';
ratecomp_wrapper(figpos,figdir,figname,write_figs,values,sigV,comp_label,sel_all,sel_resp,sep_names,ylab,tit)
%% WL
comp_label={'WnL','WL'};
values=[RatesW RatesWL];ylab='Average Rate [Hz]';figname='RatesGeneralWL';tit='Whisk & Whisk+Light Rates';
ratecomp_wrapper(figpos,figdir,figname,write_figs,values,sigV,comp_label,sel_all,sel_resp,sep_names,ylab,tit)
values=[BurstsW BurstsWL];ylab='Spikes/Burst';figname='BurstsGeneralWL';tit='Whisk & Whisk+Light Bursts';
ratecomp_wrapper(figpos,figdir,figname,write_figs,values,sigV,comp_label,sel_all,sel_resp,sep_names,ylab,tit)


%% population psths
%general params
time_before_trig=25;%ms
time_after_trig=75;%ms
binsize=2.5;%ms
safety_margin=time_after_trig;
isi=10;
mintrig=10;
binnum=(time_before_trig+time_after_trig)/binsize;
%Puff
parameter='Puff';
exclude={'Pole';'Light';'Exclude';'Grooming'};
triglength_limit=[28 32];%ms
[H_Puff,B_Puff,~,ntrigs,trigs]=get_PSTH_pop(parameter, DiscreteData,triglength_limit,ppms,time_before_trig,time_after_trig,exclude,safety_margin,binsize,isi,[]);
H_Puff(ntrigs<mintrig,:)={nan(binnum,1)};
B_Puff(ntrigs<mintrig,:)={nan(binnum,1)};


%
bins=-time_before_trig:binsize:time_after_trig;
relativeTime=bins(1:end-1)+diff(bins);
[MeanWp]=get_trig_cont_pop(trigs,relativeTime,Angle,20,DiscreteData);
%
[PuffV]=mk_vector_from_start_length(DiscreteData,'Puff');
[MeanP]=get_trig_cont_pop(trigs,relativeTime,PuffV,1,DiscreteData);
parameter='Touch';
exclude={'Puff';'Light';'Exclude';'Grooming'};
triglength_limit=[-inf inf];%ms
[H_Touch,B_Touch,~,ntrigsT,trigsT]=get_PSTH_pop(parameter, DiscreteData,triglength_limit,ppms,time_before_trig,time_after_trig,exclude,safety_margin,binsize,isi,[]);
H_Touch(ntrigs<mintrig,:)={nan(binnum,1)};
B_Touch(ntrigs<mintrig,:)={nan(binnum,1)};

bins=-time_before_trig:binsize:time_after_trig;
relativeTime=bins(1:end-1)+diff(bins);
[MeanWt]=get_trig_cont_pop(trigsT,relativeTime,Angle,20,DiscreteData);
[MeanT]=get_trig_cont_pop(trigsT,relativeTime,'Touch',1,DiscreteData);
%%
%plot
% figure('Position',[38 260 1413 525],'Color','White');Name='PuffBS';
% plot_Pop_PSTH_wrapper2(B_Puff(:,3),Name,{vpmr;pomr},{'VPM','POm'},time_before_trig,binsize,time_after_trig,{MeanW,MeanP},{'meanbefore','scaleminmax'},'baseline corr. whisker angle')
figure('Position',[38 260 1413 525],'Color','White');Name='PuffSS';
plot_Pop_PSTH_wrapper2(H_Puff,Name,{vpmr;pomr},{'VPM','POm'},time_before_trig,binsize,time_after_trig,{MeanW,MeanP},{'meanbefore','scaleminmax'},'baseline corr. whisker angle')

%
exportgraphics(gcf,[figdir 'Pop_' Name '_response_' datestr(now,'yymmdd') '.pdf'],'BackgroundColor','none')

%plot_Pop_PSTH_wrapper2(B_Puff(:,3),Name,{vpmr;pomr},{'VPM','POm'},time_before_trig,binsize,time_after_trig,H_Puff,true,'spike probability')
%
%export_fig ([figdir 'Pop_' Name '_response_' datestr(now,'yymmdd') '.pdf'], '-q95', '-pdf' ,'-transparent' ,'-native')
%%
time_before_trig=25;%ms
time_after_trig=75;%ms
binsize=2.5;%ms
safety_margin=time_after_trig;
isi=10;
mintrig=10;
binnum=(time_before_trig+time_after_trig)/binsize;
%Puff
parameter='Puff';
exclude={'Pole';'Light';'Exclude';'Grooming'};
include=[];
triglength_limit=[28 32];%ms
figure('Position',[38 260 1413 525],'Color','White');Name='PuffPopR';
plot_Pop_PSTH_wrapper_zeta(parameter, DiscreteData,triglength_limit,ppms,time_before_trig,time_after_trig,exclude,safety_margin,binsize,[],parameter,{vpmr;pomr},{'VPM';'POm'},include,{MeanWp,MeanP},{'meanbefore','scaleminmax'},'baseline corr. whisker angle')
exportgraphics(gcf,[figdir 'Pop_' parameter '_instRateResp_' datestr(now,'yymmdd') '.pdf'],'BackgroundColor','none')
%
% time_before_trig=50;%ms
% time_after_trig=100;%ms
% binsize=2.5;%ms
safety_margin=[time_before_trig time_after_trig];
isi=10;
mintrig=10;
binnum=(time_before_trig+time_after_trig)/binsize;
parameter='Touch';
exclude={'Puff';'Light';'Exclude';'Grooming'};
figure('Position',[38 260 1413 525],'Color','White');Name='TouchPopR';
triglength_limit=[-inf inf];%ms
plot_Pop_PSTH_wrapper_zeta(parameter, DiscreteData,triglength_limit,ppms,time_before_trig,time_after_trig,exclude,safety_margin,binsize,[],parameter,{vpmr;pomr},{'VPM';'POm'},include,{MeanWt,MeanT},{'meanbefore','scaleminmax'},'baseline corr. whisker angle')
exportgraphics(gcf,[figdir 'Pop_' parameter '_instRateResp_' datestr(now,'yymmdd') '.pdf'],'BackgroundColor','none')


%% Puff light
parameter='Puff';
exclude={'Pole';'Exclude';'Grooming'};
triglength_limit=[28 32];%ms
include={'Light'};
include_safety_margin=[0 40];
[H_PuffL,B_PuffL,~,ntrigsL,trigs]=get_PSTH_pop(parameter, DiscreteData,triglength_limit,ppms,time_before_trig,time_after_trig,exclude,safety_margin,binsize,isi,[],include,include_safety_margin);
%
H_PuffL(ntrigsL<mintrig,:)={nan(binnum,1)};
B_PuffL(ntrigsL<mintrig,:)={nan(binnum,1)};
%plot
bins=-time_before_trig:binsize:time_after_trig;
relativeTime=bins(1:end-1)+diff(bins);
[MeanWl]=get_trig_cont_pop(trigs,relativeTime,Angle,20,DiscreteData);
[MeanPl]=get_trig_cont_pop(trigs,relativeTime,'Puff',1,DiscreteData);
[MeanL]=get_trig_cont_pop(trigs,relativeTime,'Light',1,DiscreteData);
figure('Position',[38 260 1413 525],'Color','White');Name='Puff Light';
%plot_Pop_PSTH_wrapper2(H_PuffL,Name,{vpmr;pomr},{'VPM','POm'},time_before_trig,binsize,time_after_trig,{MeanW,MeanP,MeanL},{'meanbefore','scaleminmax','scaleminmax'},'baseline corr. whisker angle')
plot_Pop_PSTH_wrapper_zeta(parameter, DiscreteData,triglength_limit,ppms,time_before_trig,time_after_trig,exclude,safety_margin,binsize,[],parameter,{vpmr;pomr},{'VPM';'POm'},{include; include_safety_margin},{MeanWl,MeanPl,MeanL},{'meanbefore','scaleminmax','scaleminmax'},'baseline corr. whisker angle')

%plot_Pop_PSTH_wrapper2(B_Puff(:,3),Name,{vpmr;pomr},{'VPM','POm'},time_before_trig,binsize,time_after_trig,H_Puff,true,'spike probability')
%
%export_fig ([figdir 'Pop_' Name '_response_' datestr(now,'yymmdd') '.pdf'], '-q95', '-pdf' ,'-transparent' ,'-native')
exportgraphics(gcf,[figdir 'Pop_' Name '_response_' datestr(now,'yymmdd') '.pdf'],'BackgroundColor','none')

%%
vpmN=sum(ntrigsL(vpmr)>10);
pomN=sum(ntrigsL(pomr)>10);
Nmax=max(vpmN,pomN);
bins=-time_before_trig:binsize:time_after_trig;
pbins=bins(1:end-1)+diff(bins);
%% comp PSTH puff light
figure
tiledlayout(Nmax,2);
tx=-1;
for n=1:numel(vpmr)
    if ntrigsL(vpmr(n))>10
        tx=tx+2;
        nexttile(tx);
        bar(pbins,H_Puff{vpmr(n)},1,'k');hold on;
        bar(pbins,H_PuffL{vpmr(n)},1,'c','EdgeColor','none','FaceAlpha',.75);hold off;box off;
    end
end
tx=0;
for n=1:numel(pomr)
    if ntrigsL(pomr(n))>10
        tx=tx+2;
        nexttile(tx);
        bar(pbins,H_Puff{pomr(n)},1,'k');hold on;
        bar(pbins,H_PuffL{pomr(n)},1,'c','EdgeColor','none','FaceAlpha',.75);hold off;box off;
    end
end
%%
%Touch
parameter='Touch';
exclude={'Puff';'Light';'Exclude';'Grooming'};
triglength_limit=[-inf inf];%ms
[H_Touch,B_Touch,~,ntrigsT,trigsT]=get_PSTH_pop(parameter, DiscreteData,triglength_limit,ppms,time_before_trig,time_after_trig,exclude,safety_margin,binsize,isi,[]);
H_Touch(ntrigs<mintrig,:)={nan(binnum,1)};
B_Touch(ntrigs<mintrig,:)={nan(binnum,1)};

bins=-time_before_trig:binsize:time_after_trig;
relativeTime=bins(1:end-1)+diff(bins);
[MeanW]=get_trig_cont_pop(trigsT,relativeTime,Angle,20,DiscreteData);
[MeanT]=get_trig_cont_pop(trigsT,relativeTime,'Touch',1,DiscreteData);


%plot_Pop_PSTH_wrapper2(B_Puff(:,3),Name,{vpmr;pomr},{'VPM','POm'},time_before_trig,binsize,time_after_trig,H_Puff,true,'spike probability')
%
%export_fig ([figdir 'Pop_' Name '_response_' datestr(now,'yymmdd') '.pdf'], '-q95', '-pdf' ,'-transparent' ,'-native')

%plot
figure('Position',[38 260 1413 525]);Name='Touch';
plot_Pop_PSTH_wrapper2(H_Touch,Name,{vpmr;pomr},{'VPM','POm'},time_before_trig,binsize,time_after_trig,{MeanW,MeanT},{'meanbefore','scaleminmax'},'baseline corr. whisker angle')
exportgraphics(gcf,[figdir 'Pop_' Name '_response_' datestr(now,'yymmdd') '.pdf'],'BackgroundColor','none')
%export_fig ([figdir 'Pop_' Name '_response_' datestr(now,'yymmdd') '.pdf'], '-q95', '-pdf' ,'-transparent' ,'-native')
%%
%Touch
parameter='Touch';
exclude={'Light';'Exclude';'Grooming'};
triglength_limit=[-inf inf];%ms
[H_TouchEnd,B_TouchEnd,~,ntrigs]=get_PSTH_pop(parameter, DiscreteData,triglength_limit,ppms,time_before_trig,time_after_trig,exclude,safety_margin,binsize,isi,'trigend');
H_TouchEnd(ntrigs<mintrig,:)={nan(binnum,1)};
B_TouchEnd(ntrigs<mintrig,:)={nan(binnum,1)};


%plot
figure('Position',[38 260 1413 525]);Name='TouchEnd';
plot_Pop_PSTH_wrapper2(B_TouchEnd(:,3),Name,{vpmr;pomr},{'VPM','POm'},time_before_trig,binsize,time_after_trig)
export_fig ([figdir 'Pop_' Name '_End_response_' datestr(now,'yymmdd') '.pdf'], '-q95', '-pdf' ,'-transparent' ,'-native')

%%

%general params
time_before_trig=50;%ms
time_after_trig=100;%ms
binsize=2.5;%ms
safety_margin=time_after_trig;
isi=10;
%Puff
parameter='Puff';
exclude={'Whisking';'Pole';'Light';'Exclude';'Grooming'};
triglength_limit=[28 32];%ms
[H_Puff,B_Puff]=get_PSTH_pop(parameter, DiscreteData,triglength_limit,ppms,time_before_trig,time_after_trig,exclude,safety_margin,binsize,isi,[]);
%plot
figure('Position',[38 260 1413 525],'Color','White');Name='Puff non-Whisking';
%plot_Pop_PSTH_wrapper2(H_Puff,Name,{vpmr;pomr},{'VPM','POm'},time_before_trig,binsize,time_after_trig)
plot_Pop_PSTH_wrapper2(B_Puff(:,3),Name,{vpmr;pomr},{'VPM','POm'},time_before_trig,binsize,time_after_trig,H_Puff,true,'spike probability')

%export_fig ([figdir 'Pop_' Name '_PnW_response_' datestr(now,'yymmdd') '.pdf'], '-q95', '-pdf' ,'-transparent' ,'-native')
exportgraphics(gcf,[figdir 'Pop_' Name '_PnW_response_' datestr(now,'yymmdd') '.pdf'],'BackgroundColor','none','ContentType','vector')
%% Puff when Whisking
%general params
time_before_trig=50;%ms
time_after_trig=100;%ms
binsize=2.5;%ms
safety_margin=time_after_trig;
isi=10;
%Puff
parameter='Puff';
exclude={'Pole';'Light';'Exclude';'Grooming'};
include_forced='Whisking';
include_margin=[100 100];
triglength_limit=[28 32];%ms
[H_Puff,B_Puff]=get_PSTH_pop(parameter, DiscreteData,triglength_limit,ppms,time_before_trig,time_after_trig,exclude,safety_margin,binsize,isi,[],include_forced,include_margin);
%plot
figure('Position',[38 260 1413 525],'Color','White');Name='Puff Whisking';
plot_Pop_PSTH_wrapper2(H_Puff,Name,{vpmr;pomr},{'VPM','POm'},time_before_trig,binsize,time_after_trig)
%plot_Pop_PSTH_wrapper2(B_Puff(:,3),Name,{vpmr;pomr},{'VPM','POm'},time_before_trig,binsize,time_after_trig)
%
%export_fig ([figdir 'Pop_' Name '_PW_response_' datestr(now,'yymmdd') '.pdf'], '-q95', '-pdf' ,'-transparent' ,'-native')
exportgraphics(gcf,[figdir 'Pop_' Name '_PW_response_' datestr(now,'yymmdd') '.pdf'],'BackgroundColor','none','ContentType','vector')
% Puff when not Whisking
%general params
time_before_trig=50;%ms
time_after_trig=100;%ms
binsize=2.5;%ms
safety_margin=time_after_trig;
isi=10;
%Puff
parameter='Puff';
exclude={'Pole';'Whisking';'Light';'Exclude';'Grooming'};

triglength_limit=[28 32];%ms
[H_Puff,B_Puff]=get_PSTH_pop(parameter, DiscreteData,triglength_limit,ppms,time_before_trig,time_after_trig,exclude,safety_margin,binsize,isi,[]);
%plot
figure('Position',[38 260 1413 525],'Color','White');Name='Puff nonWhisking';
plot_Pop_PSTH_wrapper2(H_Puff,Name,{vpmr;pomr},{'VPM','POm'},time_before_trig,binsize,time_after_trig)
%plot_Pop_PSTH_wrapper2(B_Puff(:,3),Name,{vpmr;pomr},{'VPM','POm'},time_before_trig,binsize,time_after_trig)
% export_fig ([figdir 'Pop_' Name '_PnW_response_' datestr(now,'yymmdd') '.pdf'], '-q95', '-pdf' ,'-transparent' ,'-native')
exportgraphics(gcf,[figdir 'Pop_' Name '_PnW_response_' datestr(now,'yymmdd') '.pdf'],'BackgroundColor','none','ContentType','vector')

%% Whisking PSTH
time_before_trig=200;%ms
time_after_trig=500;%ms
binsize=10;%ms
safety_margin=time_after_trig;
isi=10;
%Whisker
parameter='Whisking';
exclude={'Pole';'Puff';'Light';'Exclude';'Grooming'};
triglength_limit=[300 inf];%ms
[H_Whisking,B_Whisking]=get_PSTH_pop(parameter, DiscreteData,triglength_limit,ppms,time_before_trig,time_after_trig,exclude,safety_margin,binsize,isi,[]);

figure('Position',[ 490   42  690  740])
Name='Whisking';
plot_Pop_PSTH_wrapper2(H_Whisking,Name,{vpmr,pomr},{'VPM','POm'},time_before_trig,binsize,time_after_trig)
export_fig ([figdir 'Pop_' Name '_response.pdf'], '-q95', '-pdf' ,'-transparent' ,'-native')



