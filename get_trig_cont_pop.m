function [mean_triggered,trial_triggered]=get_trig_cont_pop(trig,relativeTime,Name,DynamicVariableDS,DiscreteData)

if iscell(Name)
    DynamicVariable=Name;
    dvinput=true;
else
    dvinput=false;
end
trial_triggered=cell(size(trig));
mean_triggered=cell(size(trig));
for n=1:numel(trig)
    if ~isempty(trig{n})
        if dvinput && isempty(DynamicVariable{n});continue;end
        eT=linspace(0,DiscreteData(n).LengthTime,DiscreteData(n).LengthInd)';
        trigT=eT(trig{n});
        trigTest=trigT(:)+(relativeTime(:)/1000)';
        if DynamicVariableDS~=1
            DynamicVariableT=downsample(eT,DynamicVariableDS);
        else
            DynamicVariableT=eT;
        end
        if ~dvinput
            temp=zeros(size(eT));
            trigs=DiscreteData(n).([Name 'Start']);
            lengths=DiscreteData(n).([Name 'Length']);
            for t=1:numel(trigs)
                temp(trigs(t):trigs(t)+lengths(t))=1;
            end
        else
            temp=DynamicVariable{n};
        end
        trial_triggered{n}=interp1(DynamicVariableT,temp,trigTest,'linear');
        mean_triggered{n}=mean( trial_triggered{n},1,'omitnan');
    end
end

