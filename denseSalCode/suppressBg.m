% By Zhuo Su in 04/2017

function fg = suppressBg(GMM_fg, N, wei, count_)

fg = GMM_fg;

if(nargin<3)
    wei = 0.6;
end

if(nargin<4)
    count_ = 1;
end

while(1)
    if(count_>=3)
        return;
    end
    sal = fg>0.4;
    sal2 = fg<=0.4;
    if sum(sal(:))>wei*N
        fg(sal2) = 0.4;
        fg = mat2gray(fg);
        count_ = count_+1;
%         wei = wei-0.1; % Caution: easy to cause a death loop!!!!
        fg = suppressBg(fg,N,wei,count_);       
    else
        return;
    end
end