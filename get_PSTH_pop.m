
function [H,B,p,ntrigs,trig,lat,trigL]=get_PSTH_pop(parameter, DiscreteData,triglength_limit,ppms,time_before_trig,time_after_trig,exclude,safety_margin,binsize,isi,tail,trigstartend,varargin)

if ischar(parameter)
    
    [trig,trigL]=select_trigs(parameter,DiscreteData,triglength_limit,ppms,trigstartend);
    trig1=trig;
    for R=1:size(DiscreteData,2);[trig{R}]=exclude_reg_puff(trig{R});end
 %   trig=cleanup_self_safety(time_before_trig*ppms,time_after_trig*ppms,trig); %only isolated
    trig=cleanup_sections_exclusive(DiscreteData,exclude,safety_margin*ppms,trig);
    if ~isempty(varargin)
        [trig]=cleanup_sections_inclusive(DiscreteData,varargin{1},varargin{2}*ppms,trig);
    end
    for R=1:size(DiscreteData,2)
        td=trig{R}(:)-trig1{R}(:)';
        [a,b]=find(td==0);
        trigL{R}=trigL{R}(b);
    end
else
    trig=parameter(:,1);
%     for R=1:size(DiscreteData,2)
%         trig{R}=trig{R}-time_before_trig*ppms;
%     end
    trigL=parameter(:,2);
end
tempDD=struct2cell(DiscreteData);spikes=squeeze(tempDD(1,1,:));clear temp*
maxT=time_after_trig/1000;
[H,p,lat]=DD_pop_PSTH(trig,spikes,time_before_trig,time_after_trig,binsize,ppms,1,tail,maxT,trigL);
ntrigs=cellfun(@numel,trig);
% H=cellfun(@norm_psth,H,'UniformOutput',0);
if nargout>1;[B]=DD_pop_Burst_PSTH(trig,spikes,time_before_trig,time_after_trig,binsize,ppms,isi);end
