function [stage2, stage1] = RankingMMR(adjcMatrixMul, idxImg, startPos,spNum,colDistMMul)

%% Original code from:
% The core function for Manifold Ranking: 
% C. Yang, L. Zhang, H. Lu, X. Ruan, and M.-H. Yang. Saliency
% detection via graph-based manifold ranking. In CVPR, 2013.

% Code Author: Wangjiang Zhu
% Email: wangjiang88119@gmail.com
% Date: 3/24/2014

%% Changed:
%By Zhuo Su in 04/2017

alpha=0.99;
theta=10;
spNumMul = size(adjcMatrixMul, 1);
scaNum = length(idxImg);

%% Construct Super-Pixel Graph
% adjcMatrix_nn = LinkNNAndBoundary2(adjcMatrix, bdIds); 
% W = getLinkB(idxImg,spNum, inputImg);
% This super-pixels linking method is from the author's code, but is 
% slightly different from that in our Saliency Optimization

W = SetSmoothnessMatrix(colDistMMul, adjcMatrixMul, theta);
W(1:1+spNumMul:end) = 0;
%%
% The smoothness setting is also different from that in Saliency
% Optimization, where exp(-d^2/(2*sigma^2)) is used
D = diag(sum(W));
optAff =(D-alpha*W)\eye(spNumMul);
optAff(1:spNumMul+1:end) = 0;  %set diagonal elements to be zero

%% Stage 1
gI = 5;
% top
Yt=zeros(spNumMul,1);
for mk = 1:scaNum
    bst = unique(idxImg{mk}(1,:))+startPos{mk};
    Yt(bst) = 1;
end
% bst=unique(idxImg(1, :));
% Yt(bst)=1;
bsalt = ones(spNumMul,1);
if(gI~=1)
    bsalt=optAff*Yt;
    bsalt=(bsalt-min(bsalt(:)))/(max(bsalt(:))-min(bsalt(:)));
    bsalt=1-bsalt;
end
% bottom
Yb=zeros(spNumMul,1);
% bsb=unique(idxImg(end, :));
% Yb(bsb)=1;
for mk = 1:scaNum
    bsb = unique(idxImg{mk}(end,:))+startPos{mk};
    Yb(bsb) = 1;
end
bsalb = ones(spNumMul,1);
if(gI~=2)
    bsalb=optAff*Yb;
    bsalb=(bsalb-min(bsalb(:)))/(max(bsalb(:))-min(bsalb(:)));
    bsalb=1-bsalb;
end
% left
Yl=zeros(spNumMul,1);
% bsl=unique(idxImg(:, 1));
% Yl(bsl)=1;
for mk = 1:scaNum
    bsl = unique(idxImg{mk}(:,1))+startPos{mk};
    Yl(bsl) = 1;
end
bsall = ones(spNumMul,1);
if(gI~=3)
    bsall=optAff*Yl;
    bsall=(bsall-min(bsall(:)))/(max(bsall(:))-min(bsall(:)));
    bsall=1-bsall;
end
% right
Yr=zeros(spNumMul,1);
% bsr=unique(idxImg(:, end));
% Yr(bsr)=1;
for mk = 1:scaNum
    bsr = unique(idxImg{mk}(:,end))+startPos{mk};
    Yr(bsr) = 1;
end
bsalr = ones(spNumMul,1);
if(gI~=4)
    bsalr=optAff*Yr;
    bsalr=(bsalr-min(bsalr(:)))/(max(bsalr(:))-min(bsalr(:)));
    bsalr=1-bsalr;
end
% combine
stage1 = cell(1,scaNum);
stage1Mul=(bsalt.*bsalb.*bsall.*bsalr);
stage1Mul=(stage1Mul-min(stage1Mul(:)))/(max(stage1Mul(:))-min(stage1Mul(:)));
% stage1Mul = zeros(spNumMul,1);
for mk = 1:scaNum
    stage1{mk} = stage1Mul(startPos{mk}+1:startPos{mk}+spNum{mk});
%     stage1Mul(startPos{mk}+1:startPos{mk}+spNum{mk}) = stage1{mk};
end

%% Stage 2
stage2 = cell(1,scaNum);
th=mean(stage1Mul);
stage2Mul=optAff*(stage1Mul >= th);
stage2Mul = (stage2Mul-min(stage2Mul(:)))/(max(stage2Mul(:))-min(stage2Mul(:)));
for mk = 1:scaNum
    stage2{mk} = stage2Mul(startPos{mk}+1:startPos{mk}+spNum{mk});
end

% gamma = 50;
% % th = graythresh(stage2);
% th = mean(stage2Mul);
% fgProb = 1./(1+exp(-gamma*(stage2Mul - th)));
% % fgProb = stage2;
% bgProb = stage2Mul < th;

function W = SetSmoothnessMatrix(colDistM, adjcMatrix_nn, theta)
allDists = colDistM(adjcMatrix_nn > 0);
maxVal = max(allDists);
minVal = min(allDists);

colDistM(adjcMatrix_nn == 0) = Inf;
colDistM = (colDistM - minVal) / (maxVal - minVal + eps);
W = exp(-colDistM * theta);

function adjcMatrix = LinkNNAndBoundary2(adjcMatrix, bdIds)
%link boundary SPs
adjcMatrix(bdIds, bdIds) = 1;

%link neighbor's neighbor
adjcMatrix = (adjcMatrix * adjcMatrix + adjcMatrix) > 0;
adjcMatrix = double(adjcMatrix);

spNum = size(adjcMatrix, 1);
adjcMatrix(1:spNum+1:end) = 0;  %diagnal elements set to be zero