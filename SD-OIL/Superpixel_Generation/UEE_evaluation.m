clear all;clc;close all;

%%  SILC 超像素分割
% 读取图像
L = Load('superpixel_result.mat');
gt_label = imread('Label.png');

%% 区域一致性
% 初始化误差计数器
underseg_error = 0;

% 遍历每个真实标签区域
for k = 1:max(gt_label(:))
    % 提取当前区域
    region = (gt_label == k);
    
    % 找到该区域内的超像素标签
    sp_labels = unique(L(region));
    
    % 计算该区域的误差
    for sp = sp_labels'
        overlap = sum(sum(region & (SAM_angle == sp)));
        underseg_error = underseg_error + overlap * (sum(sum(SAM_angle == sp)) - overlap);
    end
end

% 归一化误差
underseg_error = underseg_error / numel(gt_label);
fprintf('Undersegmentation Error: %f\n', underseg_error);


