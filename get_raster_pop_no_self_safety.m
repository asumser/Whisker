function [H,trig]=get_raster_pop_no_self_safety(parameter,events, DiscreteData,triglength_limit,ppms,time_before_trig,time_after_trig,exclude,safety_margin,binsize,isi,trigstartend,varargin)

if ischar(parameter)
    trig=select_trigs(parameter,DiscreteData,triglength_limit,ppms,trigstartend);
    for R=1:size(DiscreteData,2)
        [trig{R}]=exclude_reg_puff(trig{R});
    end
    %trig=cleanup_self_safety(time_before_trig*ppms,time_after_trig*ppms,trig); %only isolated
    trig=cleanup_sections_exclusive(DiscreteData,exclude,safety_margin*ppms,trig);
    if ~isempty(varargin)
        [trig]=cleanup_sections_inclusive(DiscreteData,varargin{1},varargin{2}*ppms,trig);
    end
else
    trig=parameter(:,1);
  
end


%tempDD=struct2cell(DiscreteData);spikes=squeeze(tempDD(1,1,:));clear temp*
IPI=cell(size(trig));
for x=1:numel(trig)
    if ~isempty(trig{x})
        IPI{x}=cat(1,trig{x}(1),diff(trig{x}))/ppms;
    end
end
bins=(-time_before_trig:binsize:time_after_trig);
H=cell(size(trig));

for R=1:size(trig,1)
    for S=1:size(trig,2)
%          [IPIs,IPIsortIdx]=sort(IPI{R,S});
%         trig_sorted=trig{R,S}(IPIsortIdx);
        %h=(reshape(events{R},1,[])-reshape(trig_sorted,[],1))/ppms;
         h=(reshape(events{R},1,[])-reshape(trig{R,S},[],1))/ppms;
       
        a=find(h>=-time_before_trig & h<=time_after_trig);
        time=h(a);
        [trial,~]=ind2sub(size(h),a);
%         [trial,time]=find(h>=-time_before_trig & h<=time_after_trig);
%         time2=h(h>=-time_before_trig & h<=time_after_trig);
        H{R,S}=[time(:),trial(:)];
    end
    
end