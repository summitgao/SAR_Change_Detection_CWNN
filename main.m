
clear;
clc;
close all;

addpath('./Utils');
addpath('./cwnn');

PatSize = 7;
SamSize = 28;

fprintf(' ... ... read image file ... ... ... ....\n');
im1   = imread('./pic/bern_1.bmp');
im2   = imread('./pic/bern_2.bmp');
im_gt = imread('./pic/bern_gt.bmp');
fprintf(' ... ... read image file finished !!! !!!\n\n');

im1 = double(im1(:,:,1));
im2 = double(im2(:,:,1));
im_gt = double(im_gt(:,:,1));
[ylen, xlen] = size(im1);
im_di = di_gen(im1, im2);
pixel_vector = reshape(im_di, ylen*xlen, 1);

fprintf('... ... clustering begin ... ...\n');
im_lab = hclustering(pixel_vector, im_di);
fprintf('@@@ @@@ clustering finished @@@@\n');

clear pixel_vector im_di;

fprintf(' ... ... ... samples initializaton begin ... ... .....\n');
fprintf(' ... ... ... Patch Size : %d pixels ... ....\n', PatSize);


pos_idx = find(im_lab == 1);
neg_idx = find(im_lab == 0);
tst_idx = find(im_lab == 0.5);


rand('seed', 2);
pos_idx = pos_idx(randperm(numel(pos_idx)));
neg_idx = neg_idx(randperm(numel(neg_idx)));

[ylen, xlen] = size(im1);

% normalization
im1 = im1/max(im1(:));
im2 = im2/max(im2(:));


mag = (PatSize-1)/2;
imTmp = zeros(ylen+PatSize-1, xlen+PatSize-1);
imTmp((mag+1):end-mag,(mag+1):end-mag) = im1; 
im1 = im2col_general(imTmp, [PatSize, PatSize]);
imTmp((mag+1):end-mag,(mag+1):end-mag) = im2; 
im2 = im2col_general(imTmp, [PatSize, PatSize]);
clear imTmp mag;


im = zeros(SamSize, SamSize, ylen*xlen);

parfor idx = 1 : size(im1, 2)
    imtmp1 = reshape(im1(:, idx), [PatSize, PatSize]);
    imtmp2 = reshape(im2(:, idx), [PatSize, PatSize]);
    
    imtmp1 = imresize(imtmp1, [SamSize/2, SamSize], 'bilinear');
    imtmp2 = imresize(imtmp2, [SamSize/2, SamSize], 'bilinear');
    
    im(:,:,idx) = [imtmp1; imtmp2];
end
clear im1 im2 idx imtmp1 imtmp2;


sam_num = 2000;
pos_num = 1000;

if length(neg_idx) < 1000
    neg_num = length(neg_idx);
else
    neg_num = sam_num - pos_num;
end

% =================================
% --- prepare the traninng data ---
trn_pat = zeros(SamSize, SamSize, 2*sam_num);
trn_lab = zeros(2, 2*sam_num);
trn_lab(1, 1:2*neg_num) = 1;
trn_lab(2, 1+2*neg_num:2*sam_num) = 1;


for i = 1:neg_num
    trn_pat(:,:,i) = im(:,:, neg_idx(i));
    
    % virtual sample
    idx1 = ceil(neg_num*rand());
    idx2 = ceil(neg_num*rand());
    ratio = rand();
    trn_pat(:,:,neg_num+i) = im(:,:,neg_idx(idx1))*ratio + ...
        im(:,:,neg_idx(idx2))*(1-ratio);
end

for i = 1 : pos_num
    trn_pat(:,:, i+2*neg_num) = im(:,:, pos_idx(i));

    % virtual sample
    idx1 = ceil(pos_num*rand());
    idx2 = ceil(pos_num*rand());
    ratio = rand();
    trn_pat(:,:, i+2*neg_num+pos_num) = im(:,:,pos_idx(idx1))*ratio + ...
        im(:,:,pos_idx(idx2))*(1-ratio);
end
clear idx1 idx2 i ratio;

% ==================================
% --- prepare the testing data -----
tst_pat = zeros(SamSize, SamSize, pos_num);
for i = 1 : numel(tst_idx)
  tst_pat(:,:,i) = im(:,:, tst_idx(i));
end


rand('seed', 2);

cwnn.layers = {
    struct('type', 'i') 
    struct('type', 'c', 'outputmaps', 6, 'kernelsize', 5) 
    struct('type', 's', 'scale', 2)
    struct('type', 'c', 'outputmaps', 12, 'kernelsize', 5) 
    struct('type', 's', 'scale', 2) 
    struct('type', 'c', 'outputmaps', 96, 'kernelsize', 4)
};



opts.alpha = 1;
opts.batchsize = 50;
opts.numepochs = 50;

fprintf(' ... ... ... CWNN SETUP ... ... ...\n');
cwnn = cnnsetup(cwnn, trn_pat, trn_lab);
fprintf(' ... ... ... CWNN TRAIN ... ... ...\n');
cwnn = cnntrain(cwnn, trn_pat, trn_lab, opts);


fprintf('\n ... ... ... CWNN TEST  ... ... ...\n');
cwnn = cnnff(cwnn, tst_pat);

% === prepare the testing result =====
res_lab = im_lab;
tst_res = zeros(numel(tst_idx), 1);
for i = 1: numel(tst_idx)
  if cwnn.o(1,i) > cwnn.o(2,i)
    tst_res(i) = 0;
  else
    tst_res(i) = 1;
  end
end
% ====================================

res_lab(tst_idx) = tst_res';

% === refine the results =============
nos_th = 30;
[res_lab,num] = bwlabel(res_lab);
for i = 1:num
    idx = find(res_lab==i);
    if numel(idx) <= nos_th
        res_lab(idx)=0;
    end
end
res_lab = res_lab>0;

[res_lab,num] = bwlabel(~res_lab);
for i = 1:num
    idx = find(res_lab==i);
    if numel(idx) <= nos_th
        res_lab(idx)=0;
    end
end
res_lab = res_lab>0;
% =====================================

% === output the results ==============
[FA,MA,OE,CA] = DAcom(im_gt, res_lab);
fid = fopen('rec.txt', 'a');
fprintf(fid, 'FALSE ALRAMS : %d \n', FA);
fprintf(fid, 'MISSED PIXEL : %d \n', MA);
fprintf(fid, 'OVERALL ERROR: %d \n', OE);
fprintf(fid, 'PCC          : %f \n\n\n', CA);
fclose(fid);

fprintf('FALSE ALRAMS : %d \n', FA);
fprintf('MISSED PIXEL : %d \n', MA);
fprintf('OVERALL ERROR: %d \n', OE);
fprintf('PCC          : %f \n\n\n', CA);

figure; imshow(res_lab);





