%%code for
% Saliency detection via manifold ranking and dense  conditional random field 
%     by Zhuo Su; Hong Zheng; Baochang Zhang
%     Submitted to Pattern Recognition Letters

clear all
close all;clc;
addpath('SLIC');
addpath('mexGMM');

%% 1. Parameter Settings

SRC = '..\images';  
RES = '..\results';
srcSuffix = '.jpg';     %suffix for your input image

if ~exist(RES, 'dir')
    mkdir(RES);
end

%% 2. Saliency Map Calculation
files = dir(fullfile(SRC, strcat('*', srcSuffix)));

t = 0;
t_slic = 0;
t_rk = 0;
t_opt = 0;

imNum = length(files);
for kk = 1:imNum 
    disp(kk);
    srcName = files(kk).name;
    noSuffixName = srcName(1:end-length(srcSuffix));
    %% Pre-Processing: Remove Image Frames
    srcImg = imread(fullfile(SRC, srcName));
    tstart = tic;
    tic;
    if(size(srcImg,3)==1)
        srcImg = repmat(srcImg,1,1,3);
    end
    
    [noFrameImg, frameRecord] = removeframe(srcImg, 'sobel');
    [h, w, chn] = size(noFrameImg);
    
    %% Segment input rgb image into patches (SP/Grid)
    pixNumInSP = 350;                         %pixels in each superpixel
    spnumber = round( h * w / pixNumInSP );     %super-pixel number for current image
    cmpct = 20;
    
    scaNum = 2;
    idxImg = cell(1,scaNum);
    pixelList = cell(1,scaNum);

    [idxImg{1}, adjcMatrix1, pixelList{1}] = SLIC_Split(noFrameImg, spnumber, cmpct);
    
    t_slic = t_slic+toc;
    tic;
    
    meanRgbCol1 = GetMeanColor(noFrameImg, pixelList{1});
    meanLabCol1 = colorspace_ori('Lab<-', double(meanRgbCol1)/255);
    colDistM1 = GetDistanceMatrix(meanLabCol1);
    [~,idxImg{2},pixelList{2},adjcMatrixMul, meanLabColMul]...
        = MergeIdx(meanLabCol1,colDistM1,idxImg{1},pixelList{1},adjcMatrix1);

    colDistMMul = GetDistanceMatrix(meanLabColMul);
     
    spNum = cell(1,2);
    spNum{1} = length(pixelList{1});
    spNum{2} = length(pixelList{2});
    startPos = cell(1,2);
    startPos{1} = 0;
    startPos{2} = spNum{1};
     
    [stage2,stage1] = RankingMMR(adjcMatrixMul, idxImg, startPos, spNum,colDistMMul);

    stage2ave = zeros(size(idxImg{1}));
    for mk = 1:scaNum      
        smapName=fullfile(RES, strcat(noSuffixName, '_stage2_',num2str(mk),'.png')); 
        sal = SaveSaliencyMap(stage2{mk}, pixelList{mk}, frameRecord, smapName, true, false, false);
%         SaveMap(sal, frameRecord, smapName);
        stage2ave = stage2ave+sal;
    end
    stage2ave = stage2ave/scaNum;
%     smapName=fullfile(RES, strcat(noSuffixName, '_stage2_ave.png')); 
%     SaveMap(stage2ave, frameRecord, smapName);
    
    N = size(stage2ave,1)*size(stage2ave,2);     
    rate = 0.03;
    stage2ave = suppress(stage2ave,N,rate);
    rate2 = 0.6;
    stage2ave = suppressBg(stage2ave,N,rate2);
    
    t_rk = t_rk+toc;
    tic;
    
    %%  dense saliency
    GMM_fg = single(stage2ave);
    th = mean2(GMM_fg);
    GMM_bg = GMM_fg<th;
    se = strel('disk',5);   
    GMM_bg = imerode(GMM_bg,se);
    GMM_bg = single(GMM_bg);
    
    %% 71M
    iter = uint8(1);
    wei = 0.7;
    dst = GetFilteredLogFg(noFrameImg, GMM_fg, GMM_bg, GMM_fg, iter, wei);
    dst = mat2gray(dst); 
    
    radius = floor(33*sqrt(mean(dst(:))));
    dst2 = morphSmooth(dst, max(radius, 2));
    
    t_opt = t_opt+toc;
    
    t = t+toc(tstart);
    
%     smapName=fullfile(RES, strcat(noSuffixName, '_stage2_aveE.png')); 
%     SaveMap(stage2ave, frameRecord, smapName);
%     smapName=fullfile(RES, strcat(noSuffixName, '_guide.png')); 
%     SaveMap(dst, frameRecord, smapName);
    smapName=fullfile(RES, strcat(noSuffixName, '_denseSal.png')); 
    SaveMap(dst2, frameRecord, smapName);
end

avet = t/imNum;
ave_slic = t_slic/imNum;
ave_rk = t_rk/imNum;
ave_opt = t_opt/imNum;

fprintf('Total time: %0.4f\n', avet);
fprintf('SLIC: %0.4f\n', ave_slic);
fprintf('Ranking: %0.4f\n', ave_rk);
fprintf('Optimization: %0.4f\n', ave_opt);

% dataset = 'ECSSD'; 
% fid = fopen(['..\',dataset,'_record.txt'],'at');
% fprintf(fid, 'Total time: %0.4f\n', avet);
% fprintf(fid, 'SLIC: %0.4f\n', ave_slic);
% fprintf(fid, 'Ranking: %0.4f\n', ave_rk);
% fprintf(fid, 'Optimization: %0.4f\n', ave_opt);
% fclose(fid);
disp('Done!');