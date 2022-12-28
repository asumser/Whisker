function DiscreteData=update_whisking_times(DiscreteData,Amp,res,amp_thresh,min_dur)

for n=1:numel(Amp)
camp=Amp{n};
camp=camp>=amp_thresh;
wstarts=find(diff([0 camp])==1);
wends=find(diff([camp 0])==-1);
wdur=wends-wstarts;
wincl=wdur>=min_dur;
wstarts=wstarts(wincl)*res;
wdur=wdur(wincl)*res;
DiscreteData(n).WhiskingStart=wstarts;
DiscreteData(n).WhiskingLength=wdur;

end