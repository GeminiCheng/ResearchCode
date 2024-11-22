function Label = S3G(img,s, errTh, wDs)
% ����SAM�ĳ����طָ�
% imgΪ����ͼ��
% s*sΪ�����سߴ�
% errThΪ���Ƶ������������������в�����
% wDs��һ�����������ڿ�����ɫ�Ϳռ�����Ȩ��
[m,n,dim] = size(img);

%% ����դ�񶥵������ĵ�����
h = floor(m / s); % դ�������
w = floor(n / s); % դ�������
rowR = floor((m - h * s) / 2); % ���ಿ����β���ֵ���ƫ����
colR = floor((n - w * s) / 2); % ���ಿ����β���ֵ���ƫ����
rowStart = (rowR + 1) : s : (m - s + 1); % դ�����ʼ������
rowStart(1) = 1; 
rowEnd = rowStart + s; % դ��Ľ���������
rowEnd(1) = rowR + s; 
rowEnd(end) = m;
colStart = (colR + 1) : s : (n - s + 1); % դ�����ʼ������
colStart(1) = 1;
colEnd = colStart + s; % դ��Ľ���������
colEnd(1) = colR + s;
colEnd(end) = n;
rowC = floor((rowStart + rowEnd - 1) / 2); % դ�����ĵ�������
colC = floor((colStart + colEnd - 1) / 2); % դ�����ĵ�������

% ��ʾ���ֽ��
temp=zeros(m,n);
temp(rowStart,:)=1;
temp(:,colStart)=1;
for i=1:h
    for j=1:w
        temp(rowC(i),colC(j))=1;
    end
end
figure,imshow(temp);
%% �����ݶ�
gradient_data = gradient(img);
G = sqrt(sum(gradient_data.^2, 3));
%% ѡ��դ�����ĵ�3*3�������ݶ���С����Ϊ��ʼ��
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

%% �����طָ�
Label = zeros(m, n) - 1; % ��ʼ����ǩ�������ڴ洢ÿ�����صĳ����ر�ǩ
dis = Inf * ones(m, n); % ��ʼ������������ڴ洢ÿ��������������ĵľ���
M = reshape(img, m * n, dim); % ����������Ϊһ�����������������ľ���

% ��������ɫ��ֵ�Ϳ���ֵ
colorC = zeros(h, w, dim); % ��ʼ������ɫ��ֵ����
for i = 1:h
    for j = 1:w
        colorC(i, j, :) = img(rowC(i), colC(j), :); % ����ÿ�������ص����ĵ����ɫֵ
    end
end
uniMat = cat(3, colorC, rowC, colC); % ����ɫ��ֵ�Ϳ���ֵ����
uniMat = reshape(uniMat, h * w, dim + 2); % ������ɫ��ֵ�Ϳ���ֵ��������Ϊ��ά����

iter = 1; % ��ʼ����������
while(1)
    uniMat_old = uniMat; % ����ɵ�����ɫ��ֵ�Ϳ���ֵ����
    for k = 1:h * w
        c = floor((k - 1) / h) + 1; % ���㵱ǰ�������ĵ�������
        r = k - h * (c - 1); % ���㵱ǰ�������ĵ�������
        rowCidx = rowC(r, c); % ��ȡ��ǰ�������ĵ�������
        colCidx = colC(r, c); % ��ȡ��ǰ�������ĵ�������
        rowStart = max(1, rowC_std(r, c) - s); % ���㵱ǰ�������ĵ��з�Χ��ʼֵ
        rowEnd = min(m, rowC_std(r, c) + s - 1); % ���㵱ǰ�������ĵ��з�Χ����ֵ
        colStart = max(1, colC_std(r, c) - s); % ���㵱ǰ�������ĵ��з�Χ��ʼֵ
        colEnd = min(n, colC_std(r, c) + s - 1); % ���㵱ǰ�������ĵ��з�Χ����ֵ
        colorC = M((colCidx - 1) * m + rowCidx, :); % ��ȡ��ǰ�������ĵ���ɫֵ
        
        % ���������뵱ǰ�������ĵľ��룬�����¾������ͱ�ǩ����
        for i = rowStart:rowEnd
            for j = colStart:colEnd
                colorCur = M((j - 1) * m + i, :); % ��ȡ��ǰ���ص���ɫֵ
                % dc = sam_d(colorC,colorCur); % ������ɫ����
                dc = norm(colorC-colorCur);% ������ɫ����
                ds = norm(i - rowCidx, j - colCidx); % ����ռ����
                d = dc^2 + wDs * (ds / s)^2; % ���������뵱ǰ�������ĵľ���
                if d < dis(i, j)
                    dis(i, j) = d; % ���¾������
                    Label(i, j) = k; % ���±�ǩ����
                end
            end
        end
    end
    
    % ���¾�������
    colorC = zeros(h, w, dim); % ��ʼ���µ�����ɫ��ֵ����
    for k = 1:h * w
        num = 0; % ��ʼ��������
        sumColor = zeros(1, dim); % ��ʼ����ɫ��
        sumR = 0; % ��ʼ���к�
        sumC = 0; % ��ʼ���к�
        c = floor((k - 1) / h) + 1; % ���㵱ǰ�������ĵ�������
        r = k - h * (c - 1); % ���㵱ǰ�������ĵ�������
        rowCidx = rowC_std(r, c); % ��ȡ��ǰ�������ĵ�������
        colCidx = colC_std(r, c); % ��ȡ��ǰ�������ĵ�������
        rowStart = max(1, rowCidx - s); % ���㵱ǰ�������ĵ��з�Χ��ʼֵ
        rowEnd = min(m, rowCidx + s - 1); % ���㵱ǰ�������ĵ��з�Χ����ֵ
        colStart = max(1, colCidx - s); % ���㵱ǰ�������ĵ��з�Χ��ʼֵ
        colEnd = min(n, colCidx + s - 1); % ���㵱ǰ�������ĵ��з�Χ����ֵ
        
        % ����ÿ���������ĵ��µ�λ�ú���ɫֵ
        for row = rowStart:rowEnd
            for col = colStart:colEnd
                if Label(row, col) == k
                    num = num + 1; % ����������
                    sumR = sumR + row; % �����к�
                    sumC = sumC + col; % �����к�
                    color = reshape(img(row, col, :), 1, dim); % ��ȡ��ǰ���ص���ɫֵ
                    sumColor = sumColor + color; % ������ɫ��
                end
            end
        end
        colorC(r, c, :) = sumColor / num; % �����µ���ɫֵ
        rowC(r, c) = round(sumR / num); % �����µ�������
        colC(r, c) = round(sumC / num); % �����µ�������
    end
    uniMat = cat(3, colorC, rowC, colC); % ��������ɫ��ֵ�Ϳ���ֵ����
    uniMat = reshape(uniMat, h * w, dim + 2); % ��������ɫ��ֵ�Ϳ���ֵ����
    diff = uniMat - uniMat_old; % ��������ɫ��ֵ�Ϳ���ֵ����Ĳ���
    diff(:, 1:2) = sqrt(wDs) * diff(:, 1:2) / s; % �Բ����������
    err = norm(diff) / sqrt(h * w); % �������ķ���
    if err < errTh % �������С����ֵ����������
        break;
    end
end

%% ���� ���ձ߽�Ӵ��������ԭ�����С��ͨ��ı�ǩ
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
    szlabel = max(label(:)); %��ǩ����
    bh = rowEnd - rowStart + 1;
    bw = colEnd - colStart + 1;  %block�Ŀ��
    
    if szlabel < 2  %�ް�����ͨ���Թ�
        continue;
    end
    
    labelC = label(rowCidx - rowStart + 1, colCidx - colStart + 1); %����ͨ��ı��ֵ
    top = max(1, rowStart - 1);
    bottom = min(m, rowEnd + 1);
    left = max(1, colStart - 1);
    right = min(n, colEnd + 1);
    for i = 1:szlabel %������ͨ��
        if i == labelC %����ͨ�򲻴���
            continue;
        end
        marker = zeros(bottom - top + 1, right - left + 1); %����һ������һȦ��marker�������Щ���Ѿ���ͳ�ƹ��Ӵ����
        bw = label;
        bw(bw ~= i) = 0;
        bw(bw == i) = 1; %��ǰ��ͨ����ͼ
        contourBW = bwperim(bw); %��ȡ������
        idxArr = find(double(contourBW) == 1);
        labelArr = zeros(4 * length(idxArr), 1);  %��¼�������4�������ֵ������
        num = 0;
        for idx = 1:size(idxArr) %����������,ͳ����4�����ı��ֵ
            bc = floor((idxArr(idx) - 1) / bh) + 1;
            br = idxArr(idx) - bh * (bc - 1); %��������block�е�������Ϣ
            row = br + rowStart - 1;
            col = bc + colStart - 1; %�������ڴ�ͼ�е�������Ϣ
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
                if marker(row - top + 1, col - left + 1) == 0 %δ��ͳ�ƹ�
                    marker(row - top + 1, col - left + 1) = 1;
                    num = num + 1;
                    labelArr(num) = Label(row, col);
                end
            end
        end
        
        labelArr(find(labelArr == 0)) = []; %ȥ����Ԫ��
        uniqueLabel = unique(labelArr);
        numArr = zeros(length(uniqueLabel), 1);
        for p = 1:length(uniqueLabel)
            idx = find(labelArr == uniqueLabel(p));
            numArr(p) = length(idx);
        end
        idx = find(numArr == max(numArr));
        maxnumLabel = uniqueLabel(idx(1)); %�Ӵ����ı�ǩ
        
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

% ��ʾ��ͨ����������
temp = mod(Label, 20) + 1;
figure;
%imagesc(label2rgb(temp - 1, 'jet', 'w', 'shuffle'));
%axis image; axis off;

imagesc(temp);
colormap gray; % ���� colormap Ϊ�Ҷ�
%colorbar; % ��ʾ��ɫ��
