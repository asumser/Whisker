function [varargout]=cleanup_sections_exclusive(DD,exclude,safety,varargin)

varargout=cell(size(varargin));
DDcell=struct2cell(DD);

DDnames=fieldnames(DD);
Snames=fieldnames(DD(1).Sections);
S=zeros(size(Snames));
D=zeros(size(DDcell,1),1);
Dl=zeros(size(DDcell,1),1);
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





%%
% 
% varargout=cell(size(varargin));
% DDcell=struct2cell(DD);
% 
% DDnames=fieldnames(DD);
% Snames=fieldnames(DD(1).Sections);
% S=zeros(size(Snames));
% D=zeros(size(DDcell,1),1);
% Dl=zeros(size(DDcell,1),1);
% for N=1:numel(exclude)
%     if any(strcmp([exclude{N} 'Start'],DDnames))
%         D=D+strcmp([exclude{N} 'Start'],DDnames);
%         Dl=Dl+strcmp([exclude{N} 'Length'],DDnames); 
%     else
%        S=S+strcmp(exclude{N},Snames);
%     end
% end
% 
% for a=1:numel(varargin)
%     E=varargin{a};
%     for b=1:numel(E)
%         in=E{b}(:)';
%         indexin1=[cell2mat(DDcell(find(D),1,b)')' cell2mat(DDcell(find(D),1,b)')'+cell2mat(DDcell(find(Dl),1,b)')'];
%         Sections=struct2cell(DD(b).Sections); 
%         indexin2=cell2mat(Sections(find(S)));
%         indexin=[indexin1; indexin2];
%         if ~isempty(indexin)
%         indexin=indexin+[-safety(1) safety(end)];
%         ins=~any(in>indexin(:,1) & in<indexin(:,2),1);
%         
% %         decline=zeros(size(in));
% %         for n=1:size(indexin,1)
% %             decline=decline+(in>=indexin(n,1)-safety(1) & in<=indexin(n,2)+safety(end));
% %         end
%         %varargout{a}{b}=in(decline==0);
%         varargout{a}{b}=in(ins);
%         else
%             varargout{a}{b}=in;
%         end
%     end
%     varargout{a}=varargout{a}';
% end