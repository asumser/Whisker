function [selected_sections]=select_sections(DD,exclude,safety)

DDcell=struct2cell(DD);
selected_sections=cell(size(DDcell,3),1);
DDnames=fieldnames(DD);
Snames=fieldnames(DD(1).Sections);
S=zeros(size(Snames));
D=zeros(size(DDcell,1),1);
Dl=zeros(size(DDcell,1),1);
for N=1:numel(exclude)
    if any(strcmp([exclude{N} 'Start'],DDnames))
        D=D+strcmp([exclude{N} 'Start'],DDnames);
        Dl=Dl+strcmp([exclude{N} 'Length'],DDnames); 
    else
       S=S+strcmp(exclude{N},Snames);
    end
end

    for b=1:size(DDcell,3)
   
        indexin1=[cell2mat(DDcell(find(D),1,b)')' cell2mat(DDcell(find(D),1,b)')'+cell2mat(DDcell(find(Dl),1,b)')'];
        Sections=struct2cell(DD(b).Sections); 
        indexin2=cell2mat(Sections(find(S)));
        indexin=[indexin1; indexin2];
        if ~isempty(indexin)
       indexin(1,:)=indexin(1,:)-safety(1);
        indexin(2,:)=indexin(1,:)+safety(end);
        end
        selected_sections{b}=indexin;
    end
