% By Zhuo Su in 04/2017

function SaveMap(partialImg, frameRecord, imgName)

h = frameRecord(1);
w = frameRecord(2);

top = frameRecord(3);
bot = frameRecord(4);
left = frameRecord(5);
right = frameRecord(6);

partialH = bot - top + 1;
partialW = right - left + 1;

fill_value = 0;

if partialH ~= h || partialW ~= w
    feaImg = ones(h, w) * fill_value;
    feaImg(top:bot, left:right) = partialImg;
else
    feaImg = partialImg;
end

imwrite(feaImg, imgName);