function [sections,t]=split_insections(include_sections,exclude_sections)
sections=cell(size(include_sections));
t=sections;
 %w=waitbar(0,'sections');
parfor n=1:numel(include_sections)
   % waitbar(n/numel(include_sections),w)
    in=unique(include_sections{n},'rows');
    ex=unique(exclude_sections{n},'rows');
    s=false(max(in(:)),1);
    for m=1:size(in,1)
        s(in(m,1):in(m,2))=1;
    end
    for m=1:size(ex,1)
        s(ex(m,1):ex(m,2))=0;
    end
    [a,b]=squarewave2ind_v2(s,'==',1,0);
    sections{n}=[a b];
    t{n}=sum(s);
end
%close(w)
%     ex(:,1)=ex(:,1)+1;
%     ex(:,2)=ex(:,2)-1;
%     T=[ones(size(in,1),1); 2*ones(size(in,1),1);-1*ones(size(ex,1),1);-2*ones(size(ex,1),1)];
%     [~,a]=sort(I);
%     Ts=T(a)';
%     strfind(Ts,[-1 1 2 -2])
%     
%    in2=[];
%     for m=1:size(ex,1)
%    
%         
%        [~,a]=find_ind_in_start_stop_v1(ex(m,:),in);
%        if numel(a)==2
%        
%        
%        
%        if ~isempty(a)
%            if numel(a)==2
%                in(m,:)=nan;
%            elseif a==1
%                  in(m,1)=ex(find(ex(:,2)>in(m,1),1,'first'),2)+1;
%            elseif a==2
%                in(m,2)=ex(find(ex(:,1)<in(m,2),1,'last'),2)+1;
%            end
%            
%        end
%         
%         
%     end
%     sections{n}=include_sections;
% end