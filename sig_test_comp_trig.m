function [p]=sig_test_comp_trig(trigs,DiscreteData,bins,testbin,tail,ppms)
p=nan(size(DiscreteData,2),1);
for r=1:size(DiscreteData,2)
    events=DiscreteData(r).Spikes(:)';
    Nxt=cell(2,1);
    for t=1:2
        trig=trigs{t}{r}(:);

        T=repmat(trig,1,numel(events));
        E=repmat(events,numel(trig),1);
        M=E-T;
        M(M<=min(bins)*ppms | M>=max(bins)*ppms)=nan;
        N=discretize(M,bins*ppms);
        Nx=[];for n=1:numel(bins)-1;Nx=[Nx sum(N==n,2)];end
        Nx=1000*Nx./repmat(diff(bins),size(Nx,1),1);
        if ~ischar(testbin)
            Nxt{t}=Nx(:,testbin);
        elseif strcmp(testbin,'diff')
            Nxt{t}=diff(Nx,1,2);
        elseif strcmp(testbin,'quot')
            Nxt{t}=Nx(:,2)./Nx(:,1);
        end
    end
    if all(~cellfun(@isempty,Nxt))
    [p(r)]=ranksum(Nxt{1},Nxt{2},'Tail',tail);
    end
end
