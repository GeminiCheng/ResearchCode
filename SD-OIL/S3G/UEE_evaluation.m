clear all;clc;close all;

%%  SILC �����طָ�
% ��ȡͼ��
L = Load('superpixel_result.mat');
gt_label = imread('Label.png');

%% ����һ����
% ��ʼ����������
underseg_error = 0;

% ����ÿ����ʵ��ǩ����
for k = 1:max(gt_label(:))
    % ��ȡ��ǰ����
    region = (gt_label == k);
    
    % �ҵ��������ڵĳ����ر�ǩ
    sp_labels = unique(L(region));
    
    % �������������
    for sp = sp_labels'
        overlap = sum(sum(region & (SAM_angle == sp)));
        underseg_error = underseg_error + overlap * (sum(sum(SAM_angle == sp)) - overlap);
    end
end

% ��һ�����
underseg_error = underseg_error / numel(gt_label);
fprintf('Undersegmentation Error: %f\n', underseg_error);


