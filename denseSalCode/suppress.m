% By Zhuo Su in 04/2017

function fg = suppress(GMM_fg, N, wei)

% fgU = uint8(GMM_fg*255);
fg = GMM_fg;

% S_N1(boundsp) = S_N1(boundsp) - 0.6;
% neg_Ind = find(S_N1 < 0);
% if numel(neg_Ind) > 0
%    S_N1(neg_Ind) = 0.001; 
% end
while(1)
    most_sal = fg > 0.93;
    if sum(most_sal(:)) < wei*N
        fg(most_sal) = 0;
        fg = mat2gray(fg);
        wei = wei-0.004;
        fg = suppress(fg,N,wei);       
        fg(most_sal) = 1;
    %     sal_diff = setdiff(1:supNum, most_sal_sup);
    %     S_N1(sal_diff) = normalization(S_N1(sal_diff), 0);
    else
        return;
    end
end

% gamma = 10;
% th = mean2(GMM_fg);
% fg = 1./(1+exp(-gamma*(GMM_fg - th)));
% imshow(fg);


% [a, b] = imhist(fgU);
% cut = 50;
% cutA = a(cut+1:end);
% [maxA, locA] = sort(cutA, 'descend');
% b = uint8(b);
% b1 = b(locA(1)+cut);
% b2 = b(locA(2)+cut);
% b3 = b(locA(3)+cut);
% % locA = locA+cut;
% % totalA = sum(maxA);
% th = mean2(fgU);
% fgU(fgU>220) = b1;

