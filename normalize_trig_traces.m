function [normData]=normalize_trig_traces(inputData,t_dim,baseline_indx)

normData=cell(size(inputData));
for d=1:numel(inputData)
    if ~isempty(inputData{d})
temp=inputData{d};
if t_dim==1;temp=temp';end
temp=temp-median(temp(:,baseline_indx),2);
temp=temp./std(temp,0,2);
if t_dim==1;temp=temp';end
normData{d}=temp;
    end
end