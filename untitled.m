ad=discretize(amp(randspontIndW),ampbins);
spV=false(size(amp));
spikesWind=interp1(tl,1:numel(tl),tlE(DiscreteData(r).Spikes),'nearest');
spV(spikesWind)=true;
spVs=spV(randspontIndW);
ratex=nan(1,numel(ampbinsp));
for ab=1:numel(ampbinsp)
    sel=ad==ab;
    ratex(ab)=1000*sum(spVs(sel))/sum(sel);

end
%%
sel=amp(randspontIndW)>amp_thresh;
Nq=quantileranks(Wphase_spont(sel),12);
selPhase=Wphase_spont(sel);
pp=nan(1,12);
spV=false(size(amp));
spikesWind=interp1(tl,1:numel(tl),tlE(DiscreteData(r).Spikes),'nearest');
spV(spikesWind)=true;
spV=spV(randspontIndW);
spVs=spV(sel);
ratex=nan(1,12);
for x=1:12
    pp(x)=mean(selPhase(Nq==x))/pi;
ratex(x)=1000*sum(spVs((Nq==x)))./sum(Nq==x);

end
figure
plot(pp,ratex);hold on
plot(phasebinplot/pi,Phase_D{r})
%%
% figure
meanfreq=nan(95,1);
for n=1:95
instfrq = 1000/(2*pi)*diff(unwrap(Phase{n}));
meanfreq(n)=mean(instfrq(Amp{n}>3));
end
%%
1/mean(meanfreq([vpmr;pomr]))
