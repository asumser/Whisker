function [NumSD,rates_to_plot,plot_x,name_vars,PhaseMod,whisker_corrs,name_whisker_corrs,neuron_correlations_with_params]=get_whisking_in_air_modulation(bins,min_t,DiscreteData,exclude,include,safety_marginEx,safety_marginIn,ppms,amp_thresh,Amp,Setp,Angle,Phase)


DDc=struct2cell(DiscreteData);
AllSpikes=squeeze(DDc(1,:,:));
SpikesSelected=cleanup_sections_exclusive(DiscreteData,exclude,safety_marginEx*ppms,AllSpikes);
if ~isempty(include)
    [SpikesSelected]=cleanup_sections_inclusive(DiscreteData,include,safety_marginIn*ppms,SpikesSelected);
end
[spont_idx]=get_all_spont_idx(DiscreteData,exclude,safety_marginEx*ppms,safety_marginIn*ppms,include);

phasebins=linspace(-pi,pi,bins{1}+1);
ampbins=bins{2};%0:3:65;
ampbinsp=ampbins(1:end-1)+1;
abssetpbins=bins{3};%0:3:40;
abssetpbinsp=abssetpbins(1:end-1)+1;
absangbins=bins{4};%0:3:65;
absangbinsp=absangbins(1:end-1)+1;
setpbins=bins{5};%-30:3:40;
setpbinsp=setpbins(1:end-1)+1;
angbins=bins{6};%-40:3:65;
angbinsp=angbins(1:end-1)+1;
wh_vars=cell(size(bins));
Rateamp=nan(size(DiscreteData,2),numel(ampbins)-1);
Ratesetp=nan(size(DiscreteData,2),numel(setpbins)-1);
Rateangle=nan(size(DiscreteData,2),numel(angbins)-1);
Rateabssetp=nan(size(DiscreteData,2),numel(abssetpbins)-1);
Rateabsangle=nan(size(DiscreteData,2),numel(absangbins)-1);
phasebinplot=phasebins(1:end-1)+mean(diff(phasebins))/2;
Ratephase=nan(size(DiscreteData,2),numel(phasebinplot));
ft =  fittype('mean_r + amp*cos(phase-pref_phase)', 'dependent',{'lam'},'independent',{'phase'},'coefficients',{'mean_r','amp','pref_phase'});
fo=fitoptions(ft);
PhaseMod.fitresult=cell(size(DiscreteData));
PhaseMod.ftgof=PhaseMod.fitresult;
Phase_D=PhaseMod.fitresult;
PhaseMod.mean_r=zeros(size(DiscreteData));
PhaseMod.amp=PhaseMod.mean_r;
PhaseMod.pref_phase=PhaseMod.mean_r;
PhaseMod.MD=PhaseMod.mean_r;
PhaseMod.SNR=PhaseMod.mean_r;
PhaseMod.p=PhaseMod.mean_r;
NumSD=cell(size(bins));
NumSD{1}=nan(size(DiscreteData,2),numel(ampbinsp));
NumSD{2}=nan(size(DiscreteData,2),numel(setpbinsp));
NumSD{3}=nan(size(DiscreteData,2),numel(abssetpbinsp));
NumSD{4}=nan(size(DiscreteData,2),numel(angbinsp));
NumSD{5}=nan(size(DiscreteData,2),numel(absangbinsp));
NumSD{6}=nan(size(DiscreteData,2),numel(phasebinplot));
name_vars={'Amplitude';'Setpoint';'|Setpoint|';'Angle';'|Angle|';'Phase'};
%
w=waitbar(0);
for r=1:size(DiscreteData,2)
    if r<80 || r>81
        tl=linspace(0,DiscreteData(r).LengthTime,numel(Angle{r}));
        tlE=linspace(0,DiscreteData(r).LengthTime,DiscreteData(r).LengthInd);
        
        spont_sp=SpikesSelected{r};
        spikesIndW=interp1(tl,1:numel(tl),tlE(spont_sp),'nearest');
        randspontIndW=unique(interp1(tl,1:numel(tl),tlE(spont_idx{r}),'nearest'));
        if numel(randspontIndW)>0
            amp=Amp{r};
            setp=Setp{r};setp=setp-median(setp,'omitnan');
            abssetp=abs(setp);
            Wangle=Angle{r};Wangle=Wangle-median(Wangle,'omitnan');
            absWangle=abs(Wangle);
            
            [Amp_SD]=histcounts(amp(randspontIndW),ampbins);
            NumSD{1}(r,:)=Amp_SD;
            Amp_SD(Amp_SD<min_t)=nan;
            [Amp_XD]=histcounts(amp(spikesIndW),ampbins);
            Amp_D=1000*Amp_XD./Amp_SD;
            Rateamp(r,:)=Amp_D;
            
            [absSetp_SD]=histcounts(abssetp(randspontIndW),abssetpbins);
            NumSD{3}(r,:)=absSetp_SD;
            absSetp_SD(absSetp_SD<min_t)=nan;
            [absSetp_XD]=histcounts(abssetp(spikesIndW),abssetpbins);
            absSetp_D=1000*absSetp_XD./absSetp_SD;
            Rateabssetp(r,:)=absSetp_D;
            
            [Setp_SD]=histcounts(setp(randspontIndW),setpbins);
            NumSD{2}(r,:)=Setp_SD;
            Setp_SD(Setp_SD<min_t)=nan;
            [Setp_XD]=histcounts(setp(spikesIndW),setpbins);
            Setp_D=1000*Setp_XD./Setp_SD;
            Ratesetp(r,:)=Setp_D;
            
            [absAngle_SD]=histcounts(absWangle(randspontIndW),absangbins);
            NumSD{5}(r,:)=absAngle_SD;
            absAngle_SD(absAngle_SD<min_t)=nan;
            [absAngle_XD]=histcounts(absWangle(spikesIndW),absangbins);
            absAngle_D=1000*absAngle_XD./absAngle_SD;
            Rateabsangle(r,:)=absAngle_D;
            
            [Angle_SD]=histcounts(Wangle(randspontIndW),angbins);
            NumSD{4}(r,:)=Angle_SD;
            Angle_SD(Angle_SD<min_t)=nan;
            [Angle_XD]=histcounts(Wangle(spikesIndW),angbins);
            Angle_D=1000*Angle_XD./Angle_SD;
            Rateangle(r,:)=Angle_D;
            
            Wphase_spont=Phase{r}(randspontIndW);
            
            Phase_SD=histcounts(Wphase_spont(amp(randspontIndW)>amp_thresh),phasebins);
            NumSD{6}(r,:)=Phase_SD;
            spW=spikesIndW(amp(spikesIndW)>amp_thresh);
            Phase_XD=histcounts(Phase{r}(spW),phasebins);
            Phase_XD(Phase_SD<min_t/10)=nan;
            Phase_D{r}=1000*Phase_XD./Phase_SD;
            try
                pval= circ_kuipertest(Phase{r}(spW), Wphase_spont(amp(randspontIndW)>amp_thresh), 200, false);
            catch
                pval=nan;
            end
            Ratephase(r,:)=Phase_D{r};
            if sum(~isnan(Phase_D{r}))>=numel(phasebins)/2
                fo.StartPoint=[mean(Phase_D{r},'omitnan') range(Phase_D{r})/2 0];
                fo.Lower=[0 0 -pi];
                fo.Upper=[inf inf pi];
                
                [PhaseMod.fitresult{r}, PhaseMod.ftgof{r}] = fit( reshape(phasebinplot(~isnan(Phase_D{r})),[],1), reshape(Phase_D{r}(~isnan(Phase_D{r})),[],1), ft, fo);
                PhaseMod.mean_r(r)=PhaseMod.fitresult{r}.mean_r;
                PhaseMod.amp(r)=PhaseMod.fitresult{r}.amp;
                PhaseMod.pref_phase(r)=PhaseMod.fitresult{r}.pref_phase;
                PhaseMod.MD(r)=2*PhaseMod.amp(r)/PhaseMod.mean_r(r);
                PhaseMod.SNR(r)= PhaseMod.MD(r).*sqrt(PhaseMod.mean_r(r) * .111);
                PhaseMod.p(r)=pval;
            end
            wh_vars{1}=cat(1,wh_vars{1},amp(randspontIndW)');
            wh_vars{2}=cat(1,wh_vars{2},setp(randspontIndW)');
            wh_vars{3}=cat(1,wh_vars{3},abssetp(randspontIndW)');
            wh_vars{4}=cat(1,wh_vars{4},Wangle(randspontIndW)');
            wh_vars{5}=cat(1,wh_vars{5},absWangle(randspontIndW)');
            wh_vars{6}=cat(1,wh_vars{6},Phase{r}(randspontIndW)');
        end
        waitbar(r/size(DiscreteData,2),w)
    end
end
close(w)

rates_to_plot={Rateamp;Ratesetp;Rateabssetp;Rateangle;Rateabsangle;Ratephase};
plot_x={ampbinsp;setpbinsp;abssetpbinsp;angbinsp;absangbinsp;phasebinplot};


wh_varsc=cat(2,wh_vars{:});
wh_varsc(any(isnan(wh_varsc),2),:)=[];
whisker_corrs=corr(wh_varsc);
name_whisker_corrs=categorical(name_vars,name_vars,name_vars);

neuron_correlations_with_params=cell(size(rates_to_plot));
for v=1:numel(rates_to_plot)
    corrA=nan(size(DiscreteData,2),2);
    for r=1:size(DiscreteData,2)
        if any(~isnan(rates_to_plot{v}(r,:)))
            [corrA(r,1),corrA(r,2)]=corr(plot_x{v}(~isnan(rates_to_plot{v}(r,:)))',rates_to_plot{v}(r,~isnan(rates_to_plot{v}(r,:)))');
        end
    end
    neuron_correlations_with_params{v}=corrA;
end
