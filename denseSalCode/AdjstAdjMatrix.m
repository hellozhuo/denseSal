function adjaM = AdjstAdjMatrix(adjcMatrix, bdIds)

% adjcMatrix = (adjcMatrix * adjcMatrix + adjcMatrix) > 0;

adjcMatrix(bdIds,bdIds) = 1;

adjaM = full(adjcMatrix);
spNum = size(adjaM,1);
adjaM(1:spNum+1:end) = 0;%set diagonal elements to 0
adjaM = uint8(adjaM);