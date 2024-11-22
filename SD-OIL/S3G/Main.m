close all;clear all;clc
%% ��ʼ������
load('Img.mat');
img = im2double(Img36);
[m,n,z] = size(img);
M=reshape(img,m*n,z);
%% �����طָ�
s=32;
errTh=10^-2;
wDs=0.5^2;
SAM_angle = S3G(img,s, errTh, wDs);
%% ��ʾ����
marker = zeros(m,n);
for i = 1:m
    for j = 1:n
        top = SAM_angle(max(1,i-1),j);
        bottom = SAM_angle(min(m,i+1),j);
        left = SAM_angle(i,max(1,j-1));
        right= SAM_angle(i,min(n,j+1));
        if ~(top==bottom && bottom==left && left==right)
            marker(i,j)=1;
        end
    end
end
figure,imshow(marker);
%% ���ָ�����ʾ
Result_Img = imread('Result_Img.tif');
for i=1:m
    for j=1:n
        if marker(i,j)==1
            Result_Img(i,j,:)=[255,0,0];
        end
    end
end
figure,imshow(Result_Img);
title('�����طָ���');
%% �����ؽ������
SAM = Normalize(SAM_angle);
SAM = uint8(SAM);
SAM_Img = cat(3,SAM,SAM,SAM);
for i=1:m
    for j=1:n
        if marker(i,j)==1
            SAM_Img(i,j,:)=[255,0,0];
        end
    end
end
figure,imshow(SAM_Img);
title('�����طָ���');
imwrite(SAM,'SAM.tif');