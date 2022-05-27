
function [trig,trigL]=select_trigs(parameter,DiscreteData,triglength_limit,ppms,varargin)
%select trigs from DD
DDcell=struct2cell(DiscreteData);
DDnames=fieldnames(DiscreteData);
trig=squeeze(DDcell(strcmp([parameter 'Start'],DDnames),1,:));
triglength=squeeze(DDcell(strcmp([parameter 'Length'],DDnames),1,:));
%convert lenth 2 ind
triglength_limit=triglength_limit*ppms;
trigL=cell(size(trig));
for R=1:size(DiscreteData,2);
    t=trig{R};
    tl=triglength{R};
    if strcmp('trigend',varargin{1})
        t=t+tl;
    end
    t=t(tl>=triglength_limit(1) & tl<=triglength_limit(2)); % select appropriate trigs
    trig{R}=t;
    trigL{R}=tl(tl>=triglength_limit(1) & tl<=triglength_limit(2));
end


end
