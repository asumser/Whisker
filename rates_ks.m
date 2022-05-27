function [sig,p,Nx]=rates_ks(trig,events,bins,baselinebin,tail,ppms,varargin)
if strcmp(tail,'none')
    sig=nan;
    p=nan;
    Nx=nan;
    
elseif strcmp(tail,'zeta')
    maxT=varargin{1};
    trigL=varargin{2};
    [p,sZ,sR,Nx] = zetatest(events(:)/(ppms*1000),[trig(:)/(ppms*1000) trig(:)/(ppms*1000)+trigL(:)/(ppms*1000)],maxT,100,0,4,[-inf inf],true);
    %[p,sZ,Nx,sR] = getZeta(events(:)/(ppms*1000),trig(:)/(ppms*1000),maxT,100,4,4,[-maxT maxT],false);
    %Nx.peakRelative=log2(sR.dblPeakRate./mean(sR.vecRate));
    if ~isempty(sR)
        Nx.peakRate=sR.dblPeakRate;
        Nx.meanRate=sR.dblMeanRate;
        Nx.peakTime=sR.dblPeakTime;
        Nx.peakOnsetTime=sR.dblOnset;
        Nx.vecT=sR.vecT;
        Nx.vecRate=sR.vecRate;
        
    else
        Nx.peakRate=[];
        Nx.meanRate=[];
        Nx.peakTime=[];
        Nx.peakOnsetTime=[];
        Nx.vecT=[];
        Nx.vecRate=[];
    end
    sig=p<.05;
    
else
    trig=trig(:);
    events=events(:)';
    T=repmat(trig,1,numel(events));
    E=repmat(events,numel(trig),1);
    M=E-T;
    M(M<=min(bins)*ppms | M>=max(bins)*ppms)=nan;
    N=discretize(M,bins*ppms);
    Nx=[];for n=1:numel(bins)-1;Nx=[Nx sum(N==n,2)];end
    Nx=1000*Nx./repmat(diff(bins),size(Nx,1),1);
    sig=nan(1,numel(bins)-1);
    p=sig;
    for n=1:numel(bins)-1
        if any(~isnan(Nx(:,n)))
            if strcmp(tail,'larger') || strcmp(tail,'smaller') || strcmp(tail,'unequal')
                [sig(n),p(n)]=kstest2(Nx(:,baselinebin),Nx(:,n),'Tail',tail);
            elseif strcmp(tail,'right') || strcmp(tail,'left') || strcmp(tail,'both')
                [p(n),sig(n)]=signrank(Nx(:,baselinebin),Nx(:,n),'Tail',tail);
            end
        end
    end
end
