function [st,en,l]=squarewave2ind_v2(sqwave,relation,threshold,minLength)

d=eval(['diff(sqwave' relation 'threshold)']);
 if eval(['sqwave(1)' relation 'threshold']);d(1)=1;end
  if eval(['sqwave(end)' relation 'threshold']);d(end)=-1;end
st=find(d==1);
en=find(d==-1);
if ~isempty(en)
   
    while en(1)<st(1); en(1)=[]; end
    if numel(st)>numel(en)
        st(end)=[];
    elseif numel(st)<numel(en)
        en(1)=[];
    end
    l=en-st;
    ex=find(l<minLength);
    st(ex)=[];
    en(ex)=[];
    l(ex)=[];
else
    st=[];
    en=[];
    l=[];
    
end