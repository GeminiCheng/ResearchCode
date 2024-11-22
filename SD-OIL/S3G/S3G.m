function Label = S3G(img,s, errTh, wDs)
% 基于SAM的超像素分割
% img为输入图像
% s*s为超像素尺寸
% errTh为控制迭代结束的联合向量残差上限
% wDs是一个参数，用于控制颜色和空间距离的权重
[m,n,dim] = size(img);

%% 计算栅格顶点与中心的坐标
h = floor(m / s); % 栅格的行数
w = floor(n / s); % 栅格的列数
rowR = floor((m - h * s) / 2); % 多余部分首尾均分的行偏移量
colR = floor((n - w * s) / 2); % 多余部分首尾均分的列偏移量
rowStart = (rowR + 1) : s : (m - s + 1); % 栅格的起始行索引
rowStart(1) = 1; 
rowEnd = rowStart + s; % 栅格的结束行索引
rowEnd(1) = rowR + s; 
rowEnd(end) = m;
colStart = (colR + 1) : s : (n - s + 1); % 栅格的起始列索引
colStart(1) = 1;
colEnd = colStart + s; % 栅格的结束列索引
colEnd(1) = colR + s;
colEnd(end) = n;
rowC = floor((rowStart + rowEnd - 1) / 2); % 栅格中心的行坐标
colC = floor((colStart + colEnd - 1) / 2); % 栅格中心的列坐标

% 显示划分结果
temp=zeros(m,n);
temp(rowStart,:)=1;
temp(:,colStart)=1;
for i=1:h
    for j=1:w
        temp(rowC(i),colC(j))=1;
    end
end
figure,imshow(temp);
%% 计算梯度
gradient_data = gradient(img);
G = sqrt(sum(gradient_data.^2, 3));
%% 选择栅格中心点3*3邻域中梯度最小点作为起始点
rowC_std = repmat(rowC', [1, w]);
colC_std = repmat(colC, [h, 1]);
rowC = rowC_std;
colC = colC_std;
for i = 1:h
    for j = 1:w
        block = G(rowC(i, j) - 1:rowC(i, j) + 1, colC(i, j) - 1:colC(i, j) + 1);
        [~, idx] = min(block(:));
        jOffset = floor((idx + 2) / 3);
        iOffset = idx - 3 * (jOffset - 1);
        rowC(i, j) = rowC(i, j) + iOffset;
        colC(i, j) = colC(i, j) + jOffset;
    end
end

%% 超像素分割
Label = zeros(m, n) - 1; % 初始化标签矩阵，用于存储每个像素的超像素标签
dis = Inf * ones(m, n); % 初始化距离矩阵，用于存储每个像素与聚类中心的距离
M = reshape(img, m * n, dim); % 将数据重排为一个像素在行中连续的矩阵

% 计算联合色域值和空域值
colorC = zeros(h, w, dim); % 初始化联合色域值矩阵
for i = 1:h
    for j = 1:w
        colorC(i, j, :) = img(rowC(i), colC(j), :); % 计算每个超像素的中心点的颜色值
    end
end
uniMat = cat(3, colorC, rowC, colC); % 联合色域值和空域值矩阵
uniMat = reshape(uniMat, h * w, dim + 2); % 将联合色域值和空域值矩阵重排为二维矩阵

iter = 1; % 初始化迭代次数
while(1)
    uniMat_old = uniMat; % 保存旧的联合色域值和空域值矩阵
    for k = 1:h * w
        c = floor((k - 1) / h) + 1; % 计算当前聚类中心的列索引
        r = k - h * (c - 1); % 计算当前聚类中心的行索引
        rowCidx = rowC(r, c); % 获取当前聚类中心的行索引
        colCidx = colC(r, c); % 获取当前聚类中心的列索引
        rowStart = max(1, rowC_std(r, c) - s); % 计算当前聚类中心的行范围起始值
        rowEnd = min(m, rowC_std(r, c) + s - 1); % 计算当前聚类中心的行范围结束值
        colStart = max(1, colC_std(r, c) - s); % 计算当前聚类中心的列范围起始值
        colEnd = min(n, colC_std(r, c) + s - 1); % 计算当前聚类中心的列范围结束值
        colorC = M((colCidx - 1) * m + rowCidx, :); % 获取当前聚类中心的颜色值
        
        % 计算像素与当前聚类中心的距离，并更新距离矩阵和标签矩阵
        for i = rowStart:rowEnd
            for j = colStart:colEnd
                colorCur = M((j - 1) * m + i, :); % 获取当前像素的颜色值
                % dc = sam_d(colorC,colorCur); % 计算颜色距离
                dc = norm(colorC-colorCur);% 计算颜色距离
                ds = norm(i - rowCidx, j - colCidx); % 计算空间距离
                d = dc^2 + wDs * (ds / s)^2; % 计算像素与当前聚类中心的距离
                if d < dis(i, j)
                    dis(i, j) = d; % 更新距离矩阵
                    Label(i, j) = k; % 更新标签矩阵
                end
            end
        end
    end
    
    % 更新聚类中心
    colorC = zeros(h, w, dim); % 初始化新的联合色域值矩阵
    for k = 1:h * w
        num = 0; % 初始化像素数
        sumColor = zeros(1, dim); % 初始化颜色和
        sumR = 0; % 初始化行和
        sumC = 0; % 初始化列和
        c = floor((k - 1) / h) + 1; % 计算当前聚类中心的列索引
        r = k - h * (c - 1); % 计算当前聚类中心的行索引
        rowCidx = rowC_std(r, c); % 获取当前聚类中心的行索引
        colCidx = colC_std(r, c); % 获取当前聚类中心的列索引
        rowStart = max(1, rowCidx - s); % 计算当前聚类中心的行范围起始值
        rowEnd = min(m, rowCidx + s - 1); % 计算当前聚类中心的行范围结束值
        colStart = max(1, colCidx - s); % 计算当前聚类中心的列范围起始值
        colEnd = min(n, colCidx + s - 1); % 计算当前聚类中心的列范围结束值
        
        % 计算每个聚类中心的新的位置和颜色值
        for row = rowStart:rowEnd
            for col = colStart:colEnd
                if Label(row, col) == k
                    num = num + 1; % 计算像素数
                    sumR = sumR + row; % 计算行和
                    sumC = sumC + col; % 计算列和
                    color = reshape(img(row, col, :), 1, dim); % 获取当前像素的颜色值
                    sumColor = sumColor + color; % 计算颜色和
                end
            end
        end
        colorC(r, c, :) = sumColor / num; % 计算新的颜色值
        rowC(r, c) = round(sumR / num); % 计算新的行坐标
        colC(r, c) = round(sumC / num); % 计算新的列坐标
    end
    uniMat = cat(3, colorC, rowC, colC); % 更新联合色域值和空域值矩阵
    uniMat = reshape(uniMat, h * w, dim + 2); % 重排联合色域值和空域值矩阵
    diff = uniMat - uniMat_old; % 计算联合色域值和空域值矩阵的差异
    diff(:, 1:2) = sqrt(wDs) * diff(:, 1:2) / s; % 对差异进行缩放
    err = norm(diff) / sqrt(h * w); % 计算差异的范数
    if err < errTh % 如果差异小于阈值，结束迭代
        break;
    end
end

%% 后处理， 按照边界接触点数最多原则分配小连通域的标签
for k = 1:h * w
    c = floor((k - 1) / h) + 1;
    r = k - h * (c - 1);
    rowCidx = rowC_std(r, c);
    colCidx = colC_std(r, c);
    rowStart = max(1, rowCidx - s);
    rowEnd = min(m, rowCidx + s - 1);
    colStart = max(1, colCidx - s);
    colEnd = min(n, colCidx + s - 1);
    block = Label(rowStart:rowEnd, colStart:colEnd);
    block(block ~= k) = 0;
    block(block == k) = 1;
    label = bwlabel(block);
    szlabel = max(label(:)); %标签个数
    bh = rowEnd - rowStart + 1;
    bw = colEnd - colStart + 1;  %block的宽高
    
    if szlabel < 2  %无伴生连通域，略过
        continue;
    end
    
    labelC = label(rowCidx - rowStart + 1, colCidx - colStart + 1); %主连通域的标记值
    top = max(1, rowStart - 1);
    bottom = min(m, rowEnd + 1);
    left = max(1, colStart - 1);
    right = min(n, colEnd + 1);
    for i = 1:szlabel %遍历连通域
        if i == labelC %主连通域不处理
            continue;
        end
        marker = zeros(bottom - top + 1, right - left + 1); %生成一个外扩一圈的marker，标记哪些点已经被统计过接触情况
        bw = label;
        bw(bw ~= i) = 0;
        bw(bw == i) = 1; %当前连通域标记图
        contourBW = bwperim(bw); %求取外轮廓
        idxArr = find(double(contourBW) == 1);
        labelArr = zeros(4 * length(idxArr), 1);  %记录轮廓点的4邻域点标记值的向量
        num = 0;
        for idx = 1:size(idxArr) %遍历轮廓点,统计其4邻域点的标记值
            bc = floor((idxArr(idx) - 1) / bh) + 1;
            br = idxArr(idx) - bh * (bc - 1); %轮廓点在block中的行列信息
            row = br + rowStart - 1;
            col = bc + colStart - 1; %轮廓点在大图中的行列信息
            rc = [row - 1, col;...
                row + 1, col;...
                row, col - 1;...
                row, col + 1];
            for p = 1:4
                row = rc(p, 1);
                col = rc(p, 2);
                if ~(row >= 1 && row <= m && col >= 1 && col <= n && Label(row, col) ~= k)
                    continue;
                end
                if marker(row - top + 1, col - left + 1) == 0 %未被统计过
                    marker(row - top + 1, col - left + 1) = 1;
                    num = num + 1;
                    labelArr(num) = Label(row, col);
                end
            end
        end
        
        labelArr(find(labelArr == 0)) = []; %去除零元素
        uniqueLabel = unique(labelArr);
        numArr = zeros(length(uniqueLabel), 1);
        for p = 1:length(uniqueLabel)
            idx = find(labelArr == uniqueLabel(p));
            numArr(p) = length(idx);
        end
        idx = find(numArr == max(numArr));
        maxnumLabel = uniqueLabel(idx(1)); %接触最多的标签
        
        for row = rowStart:rowEnd
            for col = colStart:colEnd
                if bw(row - rowStart + 1, col - colStart + 1) == 0
                    continue;
                end
                Label(row, col) = maxnumLabel;
            end
        end
    end
end

% 显示连通域处理后聚类结果
temp = mod(Label, 20) + 1;
figure;
%imagesc(label2rgb(temp - 1, 'jet', 'w', 'shuffle'));
%axis image; axis off;

imagesc(temp);
colormap gray; % 设置 colormap 为灰度
%colorbar; % 显示颜色栏
