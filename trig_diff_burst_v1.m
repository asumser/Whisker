function [dtimes,trial,eventI]=trig_diff_burst_v1(trig,events,tb,ta,ppms)

trig=trig(:);

if diff(size(events))<0;events=events';end


T=repmat(trig,1,numel(events(1,:)));
E=repmat(events(1,:),numel(trig),1);
M=E-T;
dtimes=M(M>-tb & M<ta);
[trial,eventI]=find(M>-tb & M<ta);
