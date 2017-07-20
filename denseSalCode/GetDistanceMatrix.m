function distM = GetDistanceMatrix(feature,dopca)
% Get pair-wise distance matrix between each rows in feature
% Each row of feature correspond to a sample

% Code Author: Wangjiang Zhu
% Email: wangjiang88119@gmail.com
% Date: 3/24/2014
if(nargin<2)
    dopca = false;
end

spNum = size(feature, 1);
DistM2 = zeros(spNum, spNum);

if(dopca)
    [~, fea2, vari] = pca(feature);
    fea2(:,end) = [];
%     vari(end) = [];
    for n = 1:size(fea2, 2)
        DistM2 = DistM2 + vari(n)*( repmat(fea2(:,n), [1, spNum]) - repmat(fea2(:,n)', [spNum, 1]) ).^2;
    end
else
    for n = 1:size(feature, 2)
        DistM2 = DistM2 + ( repmat(feature(:,n), [1, spNum]) - repmat(feature(:,n)', [spNum, 1]) ).^2;
    end
end
distM = sqrt(DistM2);