function [not_regular_puffs]=exclude_reg_puff(trigP)

[Cpuff]=returnBursts2(trigP,9000);%SinB(n)=mean(cellfun(@numel,Bursts));
dCPuff=cellfun(@diff,Cpuff,'UniformOutput',0);
reg=false(size(dCPuff));
for n=1:numel(dCPuff)
    if numel(dCPuff{n})>3
        if ~any(dCPuff{n}(1:3)~=4400)
            reg(n)=1;
        end
    end
    dCPuff{n}=[nan ;dCPuff{n}];
end
not_regular_puffs=cell2mat(Cpuff(~reg)');
%dCPuff=dCPuff(~reg);
% Cpuff=Cpuff(reg);
% dCPuff=dCPuff(reg);
% regular_puffs=cell2mat(Cpuff');
% not_regular_puffs=setdiff(trigP,regular_puffs);