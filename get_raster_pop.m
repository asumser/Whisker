function [H,trig]=get_raster_pop(parameter, DiscreteData,triglength_limit,ppms,time_before_trig,time_after_trig,exclude,safety_margin,binsize,isi,trigstartend,varargin)


trig=select_trigs(parameter,DiscreteData,triglength_limit,ppms,trigstartend);
for R=1:size(DiscreteData,2);[trig{R}]=exclude_reg_puff(trig{R});end
%trig=cleanup_self_safety(time_before_trig*ppms,time_after_trig*ppms,trig); %only isolated
trig=cleanup_sections_exclusive(DiscreteData,exclude,safety_margin*ppms,trig);

if ~isempty(varargin)
    [trig]=cleanup_sections_inclusive(DiscreteData,varargin{1},varargin{2}*ppms,trig);
end
tempDD=struct2cell(DiscreteData);spikes=squeeze(tempDD(1,1,:));clear temp*


bins=(-time_before_trig:binsize:time_after_trig);
H=cell(size(trig));

for R=1:size(trig,1)
    for S=1:size(trig,2)
        
        h=(reshape(spikes{R},1,[])-reshape(trig{R,S},[],1))/ppms;
        [trial,time]=find(h>=-time_before_trig & h<=time_after_trig);
        time2=h(h>=-time_before_trig & h<=time_after_trig);
        H{R,S}=[time2,trial];
    end
    
end