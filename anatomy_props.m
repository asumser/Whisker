inct=[vpmrLa;pomrLa];
for nx=1:numel(inct)
n=inct(nx);
trP=trigsPL{n};
trL=DiscreteData(n).LightStart;
trLen=DiscreteData(n).LightLength/2e4;
trC=(trL-trP)/2e4;
trC=trC-0.02;
temp=nan(size(trC,1),3);
for x=1:size(trC,1)
    if all(trC(x,:)>0)
        temp(x,1)=trC(x,1);
        temp(x,2)=trLen(1);
         temp(x,3)=trC(x,2)-(trC(x,1)+trLen(1));
    else
        I=find(trC(x,:)<0,1,'last');
        temp(x,1)=trC(x,I);
        temp(x,2)=trLen(I);
        if I<size(trC,2)
        temp(x,3)=trC(x,I+1)-(trC(x,I)+trLen(I));
        end
    end
end
latLP{nx}=temp;
end

latLPall=cat(1,latLP{:});
latLPall(latLPall(:,3)>5,3)=nan;
latLPmmm=cat(1,min(latLPall),max(latLPall),median(latLPall,'omitnan'));
%%
figure
histogram(latLPall)
figure
scatter(latLPall(:,1),latLPall(:,2))
%%
cellfun(@min,latLP,'UniformOutput',0)'
%%
WL=cat(2,DiscreteData([vpmra;pomra]).WhiskingLength)';


%%

inct=[vpmra;pomra];
FillDist=nan(numel(inct),4);
for nx=1:numel(inct)
    n=inct(nx);
    CellCoord=RecDB{n,'CellAtlasCoords'};
    AN=RecDB{n,'AnimalName'}{1};
    FillCord=ExpDB{AN,'FillAtlasCoords'};
FillDist(nx,1:3)=CellCoord-FillCord;
FillDist(nx,4)=norm(CellCoord-FillCord);

end
