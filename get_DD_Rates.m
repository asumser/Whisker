function [Rates,SinB]=get_DD_Rates(DiscreteData,exclude,include, outsafety,insafety,ppms,isi)
DDcell=struct2cell(DiscreteData);
DDnames=fieldnames(DiscreteData);
[exclude_sections]=select_sections(DiscreteData,exclude,outsafety);
if ~isempty(include)
[include_sections]=select_sections(DiscreteData,include,insafety);
else
     include_sections=cell(size(DiscreteData,2),1);
    for n=1:size(DiscreteData,2); include_sections{n}=[1 DiscreteData(n).LengthInd]; end
end
[~,time_in_sections]=split_insections(include_sections,exclude_sections);

[Spikes]=cleanup_sections_exclusive(DiscreteData,exclude,outsafety*ppms,squeeze(DDcell(1,:,:)));
if ~isempty(include);[Spikes]=cleanup_sections_inclusive(DiscreteData,include,0,Spikes);end
if nargout>1
SinB=nan(size(Spikes));
for n=1:numel(Spikes);[Bursts]=returnBursts2(Spikes{n},isi*ppms);SinB(n)=mean(cellfun(@numel,Bursts));end
end
N_Spikes=cellfun(@numel,Spikes);
Rates=(N_Spikes./cell2mat(time_in_sections))*(ppms*1000);

