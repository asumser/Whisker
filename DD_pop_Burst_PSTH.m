function [H]=DD_pop_Burst_PSTH(trig,events,tb,ta,bs,ppms,isi)

bins=(-tb:bs:ta);
H=cell(size(trig,1),3);
w=waitbar(0,'compute rates');
for R=1:numel(trig)
    [~,burstI]=returnBursts2(events{R},isi*ppms);%SinB(n)=mean(cellfun(@numel,Bursts));
    [dtimes,~,eventI]=trig_diff_burst_v1(trig{R},burstI,tb*ppms,ta*ppms,ppms);% relative times to trigs
    dtimes=dtimes/ppms;% convert 2 ms
    [h,~,B]=histcounts(dtimes,bins);% mk hist
    H{R,3}=((1000./diff(bins)).*h./numel(trig{R}))';
    SpikesPerEvent=burstI(2,eventI);
    H{R,1}=nan(numel(bins)-1,1);
    H{R,2}=nan(numel(bins)-1,1);
    for b=1:numel(bins)-1
        temp=SpikesPerEvent(B==b);
        H{R,1}(b)=mean(temp);
        H{R,2}(b)=sum(temp(temp>1))/sum(temp);
    end
    %
    % [sig,p]=rates_ks(trig,spikes,bins,baselinebin,'smaller',ppms);
    % Rates(R).PuffRatesSig=sig;
    % Rates(R).PuffRatesP=p;
    waitbar(R/numel(trig),w)
end
close(w)