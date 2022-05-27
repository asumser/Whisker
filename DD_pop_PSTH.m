function [H,p,lat]=DD_pop_PSTH(trig,spikes,tb,ta,bs,ppms,varargin)
%scale time 2 indices
% tb=tb*ppms;
% ta=ta*ppms;
% bs=bs*ppms;
% mk bins
bins=(-tb:bs:ta);
H=cell(size(trig));
%w=waitbar(0,'compute rates');


    p=cell(1,size(trig,2));
    p(:)={nan(size(trig,1),numel(bins)-1)}; 
    lat=cell(size(trig,2),size(trig,2));

for R=1:size(trig,1)
for S=1:size(trig,2)
% [dtimes]=trig_diff_v3(trig{R,S},spikes{R},tb*ppms,ta*ppms,ppms,0);% relative times to trigs
% dtimes=dtimes/ppms;% convert 2 ms
% [h]=histcounts(dtimes,bins);% mk hist
h=histcounts((reshape(spikes{R},1,[])-reshape(trig{R,S},[],1))/ppms,bins);


H{R,S}=((1000./diff(bins)).*h./numel(trig{R,S}))';
if ~isempty(varargin)
 if numel(varargin)<=2
     varargin{3}=[];
     varargin{4}=[];
 end

[~,p{S}(R,:),lat{S,R}]=rates_ks(trig{R,S},spikes{R},bins,varargin{1},varargin{2},ppms,varargin{3},varargin{4}{R});

% 
 %[sig,p]=rates_ks(trig,spikes,bins,baselinebin,'smaller',ppms);
end
% Rates(R).PuffRatesSig=sig;
% Rates(R).PuffRatesP=p;
end
% waitbar(R/size(trig,1),w)
end
% close(w)