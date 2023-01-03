function DiscreteData=mk_nowhisk_sections(DiscreteData)




for n=1:numel(DiscreteData)
    W_start=DiscreteData(n).WhiskingStart;
    if ~isempty(W_start)
        W_end=DiscreteData(n).WhiskingStart+DiscreteData(n).WhiskingLength;
        NW_start=W_end+20;
        NW_end=W_start-20;
        if W_start(1)<=20
            NW_end(1)=[];
        else
            NW_start=[1 NW_start];
        end
        if W_end(end)>=DiscreteData(n).LengthInd-20
            NW_start(end)=[];
        else
            NW_end=[NW_end DiscreteData(n).LengthInd];
        end


        DiscreteData(n).NoWhiskingStart=NW_start;
        DiscreteData(n).NoWhiskingLength=NW_end-NW_start;
    else
        DiscreteData(n).NoWhiskingStart=[];
        DiscreteData(n).NoWhiskingLength=[];
    end

end