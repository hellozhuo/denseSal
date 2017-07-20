% By Zhuo Su in 04/2017
% function [idxImgMerge,adjcMatrixMerge,pixelListMerge,meanLabColMerge]...
%     = MergeIdx(meanLabCol,colDistM,idxImg,pixelList,adjcMatrix)
function [bdIds,idxImg2,pixelList2,adjcMatrixMul, meanLabColMul]...
    = MergeIdx(meanLabCol1,colDistM,idxImg,pixelList,adjcMatrix)

spNum = size(adjcMatrix,1);
adjcMatrix(1:1+spNum:end) = 0;
pixelNum = zeros(spNum,1);
for i = 1:spNum
    pixelNum(i) = length(pixelList{i});
end

% pixelLim = max(pixelNum)+1;
info = zeros(spNum,1);
idxImg2 = zeros(size(idxImg));

spNum2 = 0;
nodes = (1:spNum);
while ~isempty(nodes)
    merged = false;
    node = nodes(1);
    adj = adjcMatrix(node,:);
    adj = find(adj~=0);
    adjDis = colDistM(node,adj);
    [~, mergeL] = sort(adjDis);
    for loc = 1:length(mergeL)
        [x,y] = ismember(adj(mergeL(loc)),nodes);
        if(x)
            spNum2 = spNum2+1;
            idxImg2([pixelList{node};pixelList{adj(mergeL(loc))}]) = spNum2;
            info(node) = spNum2;
            info(adj(mergeL(loc))) = spNum2;
            nodes([1,y]) = [];
            merged = true;
            break;
        end
    end
    if(~merged)
        spNum2 = spNum2+1;
        idxImg2(pixelList{node}) = spNum2;
        info(node) = spNum2;
        nodes(1) = [];
    end
end

bdIds = cell(1,2);
bdIds{1} = GetBndPatchIds(idxImg);
bdIds{2} = GetBndPatchIds(idxImg2);

adjcMatrixMul = zeros(spNum+spNum2,spNum+spNum2);
adjcMatrix(1:1+spNum:end) = 1;
% adjcMatrix = LinkNNAndBoundary(adjcMatrix, bdIds{1});
adjcMatrix = LinkBoundarySPs(adjcMatrix, bdIds{1});

adjcMatrix2 = GetAdjMatrix(idxImg2, spNum2);
% adjcMatrix2 = LinkNNAndBoundary(adjcMatrix2, bdIds{2});
adjcMatrix2 = LinkBoundarySPs(adjcMatrix2, bdIds{2});

adjcMatrixMul(1:spNum,1:spNum) = adjcMatrix;
adjcMatrixMul(spNum+1:end,spNum+1:end) = adjcMatrix2;
for i = 1:spNum
    adjcMatrixMul(i,info(i)+spNum) = 1;
    adjcMatrixMul(info(i)+spNum,i) = 1;
end

meanLabCol2 = zeros(spNum2, size(meanLabCol1,2));
pixelList2 = cell(spNum2,1);
for i = 1:spNum2
    adj = find(adjcMatrixMul(spNum+i,1:spNum)~=0);
    meanLabCol2(i,:) = sum(meanLabCol1(adj,:),1)/length(adj);
    sumPixel = 0;
    for j = 1:length(adj)
        sumPixel = sumPixel+length(pixelList{adj(j)});
    end
    pixelList2{i} = zeros(sumPixel,1);
    beginPos = 0;
    for j = 1:length(adj)
        pixelList2{i}(beginPos+1:beginPos+length(pixelList{adj(j)})) = pixelList{adj(j)};
        beginPos = beginPos+length(pixelList{adj(j)});
    end
end

meanLabColMul = [meanLabCol1;meanLabCol2];

% adjcMatrixList = cell(1,2);
% adjcMatrixList{1} = adjcMatrix;
% adjcMatrixList{2} = adjcMatrix2;





        
        
