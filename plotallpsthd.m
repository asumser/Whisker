figure
tiledlayout(max(numel(vpmr),numel(pomr)),2)
seln={vpmr,pomr};
psths={P_PSTH_instrate, T_PSTH_instrate};
c=[.25 .25 .25; 1 .25 .25];
for s=1:numel(seln)
for n=1:numel(seln{s})
nt(s,n)=nexttile(2*n+s-2);
for t=1:numel(psths)
    temp=psths{t}(seln{s}(n),:);
    temp=log2(temp./mean(temp(psth_relativeTime<0)));
bar(psth_relativeTime,temp,1,'FaceColor',c(t,:));hold on;
end


end
end