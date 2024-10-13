function [varargout]=cleanup_sections_exclusive(DD,exclude,safety,varargin)

varargout=cell(size(varargin));
DDcell=struct2cell(DD);
DDnames=fieldnames(DD);

for a=1:numel(varargin)
    E=varargin{a};
    for b=1:numel(E)
        Ex=zeros(1,DD(b).LengthInd);
        for N=1:numel(exclude)
            if any(strcmp([exclude{N} 'Start'],DDnames))
                D=DDcell{strcmp([exclude{N} 'Start'],DDnames),1,b};
                Dl=DDcell{strcmp([exclude{N} 'Length'],DDnames),1,b};
                for t=1:numel(D)
                    Ex(D(t):D(t)+Dl(t))=1;
                end
            else
                Stemp=DD(b).Sections.(exclude{N});
                for t=1:size(Stemp,1)
                    Ex(Stemp(t,1):Stemp(t,2))=1;
                end
            end
        end
        trigs=E{b}(:)+[-safety(1):safety(end)];
        trigs(trigs<1)=1;
         trigs(trigs>DD(b).LengthInd)=DD(b).LengthInd;
         varargout{a}{b}=E{b}(~any(Ex(trigs)>0,2));
    end
    varargout{a}=varargout{a}';
end
