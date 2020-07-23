## Task 6 3D Measurement SfM
**Name:** Tianyi Yu  
**Degree:** BA  
**ID:** 20550150

### Problem 2
---
```matlab
prob2I1 = imread(fullfile(imgFolder,cell2mat(ImgStruct.imgFileName(1))));
prob2I2 = imread(fullfile(imgFolder,cell2mat(ImgStruct.imgFileName(2))));
prob2index = [];

selectedimg = [1 2];

index1=[];
for i = 1:PointInfo.nPt3D
    a = cell2mat(PointInfo.matchImgIdx(i));
    if ismember(selectedimg,a)==[1 1]
        index1 = [index1,i];
    end
end

mp1=[];
for i = 1:PointInfo.nPt3D
    if ismember(i,index1)==1
        c = cell2mat(PointInfo.matchImgPos(i));
        mp1=[mp1;c(1,:)];
    end
end

mp2=[];
for i = 1:PointInfo.nPt3D
    if ismember(i,index1)==1
        d = cell2mat(PointInfo.matchImgPos(i));
        mp2=[mp2;d(2,:)];
    end
end

figure; ax = axes;
showMatchedFeatures(prob2I1,prob2I2,mp1,mp2,'montage','Parent',ax);
title(ax, 'Candidate point matches');
legend(ax, 'Matched points 1','Matched points 2');
```
![](Prob2.png)

### Problem 3
---
```matlab
selectedpt=1000;
prob3idx = cell2mat(PointInfo.matchImgIdx(1000))

prob3I1 = imread(fullfile(imgFolder,cell2mat(ImgStruct.imgFileName(prob3idx(1)))));
prob3I2 = imread(fullfile(imgFolder,cell2mat(ImgStruct.imgFileName(prob3idx(2)))));
prob3I3 = imread(fullfile(imgFolder,cell2mat(ImgStruct.imgFileName(prob3idx(3)))));

mp = cell2mat(PointInfo.matchImgPos(selectedpt))
mp3 = mp(1,:);
mp4 = mp(2,:);
mp5 = mp(3,:);

figure;
subplot(1,3,1); imshow(prob3I1); hold on; plot(mp3(1),mp3(2),'ro');
subplot(1,3,2); imshow(prob3I2); hold on; plot(mp4(1),mp4(2), 'ro');
subplot(1,3,3); imshow(prob3I3); hold on; plot(mp5(1),mp5(2), 'ro');
```


### Problem 4
---
```matlab
px = 3888/2;
py = 2592/2;


[f1 f2 f3 f4 f5 f6] = ImgStruct.fMat{1:6};
K1 = [f1 0 px; 0 f1 py; 0 0 1];
K2 = [f2 0 px; 0 f2 py; 0 0 1];
K3 = [f3 0 px; 0 f3 py; 0 0 1];
K4 = [f4 0 px; 0 f4 py; 0 0 1];
K5 = [f5 0 px; 0 f5 py; 0 0 1];
K6 = [f6 0 px; 0 f6 py; 0 0 1];

[R1 R2 R3 R4 R5 R6] = ImgStruct.RMat{1:6};

[C1 C2 C3 C4 C5 C6] = ImgStruct.CMat{1:6};

I = eye(3);

P1 = K1*R1*[I -C1'];
P2 = K2*R2*[I -C2'];
P3 = K3*R3*[I -C3'];
P4 = K4*R4*[I -C4'];
P5 = K5*R5*[I -C5'];
P6 = K6*R6*[I -C6'];

ImgStruct.PMat = {P1,P2,P3,P4,P5,P6};
```


### Problem 5
---
```matlab
selectedpt = 1000;
p5pt3d = PointInfo.Pt3D{selectedpt};
p5pt3d = [p5pt3d 1]';
ptidx = PointInfo.matchImgIdx{selectedpt};

prob5I1 = imread(fullfile(imgFolder,cell2mat(ImgStruct.imgFileName(ptidx(1)))));
prob5I2 = imread(fullfile(imgFolder,cell2mat(ImgStruct.imgFileName(ptidx(2)))));
prob5I3 = imread(fullfile(imgFolder,cell2mat(ImgStruct.imgFileName(ptidx(3)))));

p5pt1 = ImgStruct.PMat{ptidx(1)}*p5pt3d;
p5pt1 = [p5pt1(1)/p5pt1(3) p5pt1(2)/p5pt1(3)];
p5pt2 = ImgStruct.PMat{ptidx(2)}*p5pt3d;
p5pt2 = [p5pt2(1)/p5pt2(3) p5pt2(2)/p5pt2(3)];
p5pt3 = ImgStruct.PMat{ptidx(3)}*p5pt3d;
p5pt3 = [p5pt3(1)/p5pt3(3) p5pt3(2)/p5pt3(3)];

figure;
subplot(1,3,1); imshow(prob5I1); hold on; plot(p5pt1(1),p5pt1(2),'ro');
subplot(1,3,2); imshow(prob5I2); hold on; plot(p5pt2(1),p5pt2(2),'ro');
subplot(1,3,3); imshow(prob5I3); hold on; plot(p5pt3(1),p5pt3(2),'ro');
```
![](Prob5.png)

### Problem 6
---
```matlab
selectedimg = [1 2];

index3=[];
p6pt3d = [];
for i = 1:PointInfo.nPt3D
    a = PointInfo.matchImgIdx{i};
    if ismember(selectedimg,a)==[1 1]
        index3 = [index3,i];
        p6pt3d = [p6pt3d [PointInfo.Pt3D{i} 1]'];
    end
end

p6pta = [ImgStruct.PMat{1}*p6pt3d]';
pta(:,1) = p6pta(:,1)./p6pta(:,3);
pta(:,2) = p6pta(:,2)./p6pta(:,3);

p6ptb = [ImgStruct.PMat{2}*p6pt3d]';
ptb(:,1) = p6ptb(:,1)./p6ptb(:,3);
ptb(:,2) = p6ptb(:,2)./p6ptb(:,3);

nPt = 421;

P1 = ImgStruct.PMat{1};
P2 = ImgStruct.PMat{2};

x1 = pta';
x2 = ptb';

funRow = @(u,v,up,vp) [u*up v*up up u*vp v*vp vp u v 1];

% normalize 8 point algorithm
T1 = [2/max(x1(1,:)) 0 -1; 0 2/max(x1(2,:)) -1; 0 0 1]; 
T2 = [2/max(x2(1,:)) 0 -1; 0 2/max(x2(2,:)) -1; 0 0 1];

x1T = T1*[x1;ones(1,nPt)]; x1T = bsxfun(@rdivide, x1T(1:2,:), x1T(3,:));
x2T = T2*[x2;ones(1,nPt)]; x2T = bsxfun(@rdivide, x2T(1:2,:), x2T(3,:));
id = randperm(nPt,8); 
A = zeros(8,9);
for ii=1:8
   A(ii,:) = funRow(x1T(1,id(ii)), x1T(2,id(ii)), x2T(1,id(ii)), x2T(2,id(ii))); 
end
[~, ~, V] = svd(A);
FT = reshape(V(:,9), 3, 3)';
[U, D, V] = svd(FT);
FT = U*diag([D(1,1) D(2,2) 0])*V';

F = T2'*FT*T1

% test transpose(x')*F*x = 0 
xFx_norm = zeros(421,1);
for ii=1:421
    l = [x2(:,ii);1]'*F;
    dist = abs(l(1)*x1(1,ii) + l(2)*x1(2,ii) + l(3))/norm(l(1:2));
    xFx_norm(ii) = dist;
end

mean(xFx_norm)

% F2=[e]xP'P+
e2 = ImgStruct.PMat{2}*[ImgStruct.CMat{1} 1]';
e2x = [0 -e2(3) e2(2); e2(3) 0 -e2(1); -e2(2) e2(1) 0];
F2 = e2x*ImgStruct.PMat{2}*pinv(ImgStruct.PMat{1});

xFx_norm = zeros(421,1);
for ii=1:421
    l = [x2(:,ii);1]'*F2;
    dist = abs(l(1)*x1(1,ii) + l(2)*x1(2,ii) + l(3))/norm(l(1:2));
    xFx_norm(ii) = dist;
end

mean(xFx_norm)
```

### Problem 7
---
```matlab
selectedimg = [1 2];

index3=[];
p6pt3d = [];
for i = 1:PointInfo.nPt3D
    a = PointInfo.matchImgIdx{i};
    if ismember(selectedimg,a)==[1 1]
        index3 = [index3,i];
        p6pt3d = [p6pt3d [PointInfo.Pt3D{i} 1]'];
    end
end

p6pta = [ImgStruct.PMat{1}*p6pt3d]';
pta(:,1) = p6pta(:,1)./p6pta(:,3);
pta(:,2) = p6pta(:,2)./p6pta(:,3);

p6ptb = [ImgStruct.PMat{2}*p6pt3d]';
ptb(:,1) = p6ptb(:,1)./p6ptb(:,3);
ptb(:,2) = p6ptb(:,2)./p6ptb(:,3);

nPt = 421;

P1 = ImgStruct.PMat{1};
P2 = ImgStruct.PMat{2};

x1 = pta';
x2 = ptb';

funRow = @(u,v,up,vp) [u*up v*up up u*vp v*vp vp u v 1];


% normalize 8 point algorithm
T1 = [2/max(x1(1,:)) 0 -1; 0 2/max(x1(2,:)) -1; 0 0 1]; 
T2 = [2/max(x2(1,:)) 0 -1; 0 2/max(x2(2,:)) -1; 0 0 1];

x1T = T1*[x1;ones(1,nPt)]; x1T = bsxfun(@rdivide, x1T(1:2,:), x1T(3,:));
x2T = T2*[x2;ones(1,nPt)]; x2T = bsxfun(@rdivide, x2T(1:2,:), x2T(3,:));

A = zeros(8,9);
for ii=1:8
   A(ii,:) = funRow(x1T(1,id(ii)), x1T(2,id(ii)), x2T(1,id(ii)), x2T(2,id(ii))); 
end
[~, ~, V] = svd(A);
FT = reshape(V(:,9), 3, 3)';
[U, D, V] = svd(FT);
FT = U*diag([D(1,1) D(2,2) 0])*V';

F = T2'*FT*T1

% test transpose(x')*F*x = 0 
xFx_norm = zeros(421,1);
for ii=1:421
    l = [x2(:,ii);1]'*F;
    dist = abs(l(1)*x1(1,ii) + l(2)*x1(2,ii) + l(3))/norm(l(1:2));
    xFx_norm(ii) = dist;
end
mean(xFx)
mean(xFx_norm)

% F2=[e]xP'P+
e2 = ImgStruct.PMat{2}*[ImgStruct.CMat{1} 1]';
e2x = [0 -e2(3) e2(2); e2(3) 0 -e2(1); -e2(2) e2(1) 0];
F2 = e2x*ImgStruct.PMat{2}*pinv(ImgStruct.PMat{1});

xFx_norm = zeros(421,1);
for ii=1:421
    l = [x2(:,ii);1]'*F2;
    dist = abs(l(1)*x1(1,ii) + l(2)*x1(2,ii) + l(3))/norm(l(1:2));
    xFx_norm(ii) = dist;
end
mean(xFx)
mean(xFx_norm)

pta1 = [pta(1,:) 1]';
pta2 = [pta(2,:) 1]';

ep1 = F*pta1
ep2 = F*pta2

x = 0:0.5:3888;

p7y1 = (-ep1(1)*x-ep1(3))/ep1(2);
p7y2 = (-ep2(1)*x-ep2(3))/ep2(2);

p7pt1 = ImgStruct.PMat{2}*p6pt3d(:,1);
p7pt2 = ImgStruct.PMat{2}*p6pt3d(:,2);

pt1 = [p7pt1(1)/p7pt1(3) p7pt1(2)/p7pt1(3)];
pt2 = [p7pt2(1)/p7pt2(3) p7pt2(2)/p7pt2(3)];

%prob7I1 = imread(fullfile(uimgFolder,ImgStruct.imgFileName{1}));
prob7I2 = imread(fullfile(uimgFolder,ImgStruct.imgFileName{2}));

figure;
subplot(1,2,1); imshow(prob7I2); hold on; plot(x,p7y1,'-r'); plot(x,p7y2,'-b'); plot(pt1(1),pt1(2),'ro'); plot(pt2(1),pt2(2),'bo');
subplot(1,2,2); imshow(prob7I2); hold on; plot(pt1(1),pt1(2),'ro'); plot(pt2(1),pt2(2),'bo');
```
![](Prob7.png)

### Problem 8
---
```matlab
prob8I1 = imread(fullfile(uimgFolder,ImgStruct.imgFileName{1}));
prob8I2 = imread(fullfile(uimgFolder,ImgStruct.imgFileName{2}));
% h = cpselect(prob8I1,prob8I2);

% First 3D point
x1 = movingPoints(1,1);
x2 = movingPoints(1,2);
x3 = 1;
y1 = fixedPoints(1,1);
y2 = fixedPoints(1,2);
y3 = 1;

movingP = ImgStruct.PMat{1};
fixedP = ImgStruct.PMat{2};

coeff1 = [x3*movingP(1,:)-x1*movingP(3,:);x3*movingP(2,:)-x2*movingP(3,:);...
    y3*fixedP(1,:)-y1*fixedP(3,:);y3*fixedP(2,:)-y2*fixedP(3,:)];

[~,~,pt3d1] = svd(coeff1);
pt3d1 = pt3d1(:,end);
p8pt1 = [pt3d1(1)/pt3d1(4) pt3d1(2)/pt3d1(4) pt3d1(3)/pt3d1(4)];

% Second 3D point
xx1 = movingPoints(2,1);
xx2 = movingPoints(2,2);
xx3 = 1;
yy1 = fixedPoints(2,1);
yy2 = fixedPoints(2,2);
yy3 = 1;

coeff2 = [xx3*movingP(1,:)-xx1*movingP(3,:);xx3*movingP(2,:)-xx2*movingP(3,:);...
    yy3*fixedP(1,:)-yy1*fixedP(3,:);yy3*fixedP(2,:)-yy2*fixedP(3,:)];

[~,~,pt3d2] = svd(coeff2);
pt3d2 = pt3d2(:,end);
p8pt2 = [pt3d2(1)/pt3d1(4) pt3d2(2)/pt3d1(4) pt3d2(3)/pt3d1(4)];

% Calculate distance
d1 = norm(p8pt1-p8pt2)
scalingfactor = 2/d1

worldPoints = triangulate(movingPoints,fixedPoints,movingP',fixedP');
wpt1 = worldPoints(1,:);
wpt2 = worldPoints(2,:);
d2 = norm(wpt1-wpt2)
```
![](Prob8.png)

### Problem 9
---
```matlab
prob9I1 = imread(fullfile(uimgFolder,ImgStruct.imgFileName{1}));
prob9I2 = imread(fullfile(uimgFolder,ImgStruct.imgFileName{2}));
h  = cpselect(prob9I1,prob9I2);

% First 3D point
x1 = movingPoints1(1,1);
x2 = movingPoints1(1,2);
x3 = 1;
y1 = fixedPoints1(1,1);
y2 = fixedPoints1(1,2);
y3 = 1;

movingP = ImgStruct.PMat{1};
fixedP = ImgStruct.PMat{2};

coeff1 = [x3*movingP(1,:)-x1*movingP(3,:);x3*movingP(2,:)-x2*movingP(3,:);...
    y3*fixedP(1,:)-y1*fixedP(3,:);y3*fixedP(2,:)-y2*fixedP(3,:)];

[~,~,pt3d1] = svd(coeff1);
pt3d1 = pt3d1(:,end);
p8pt1 = [pt3d1(1)/pt3d1(4) pt3d1(2)/pt3d1(4) pt3d1(3)/pt3d1(4)];

% Second 3D point
xx1 = movingPoints1(2,1);
xx2 = movingPoints1(2,2);
xx3 = 1;
yy1 = fixedPoints1(2,1);
yy2 = fixedPoints1(2,2);
yy3 = 1;

coeff2 = [xx3*movingP(1,:)-xx1*movingP(3,:);xx3*movingP(2,:)-xx2*movingP(3,:);...
    yy3*fixedP(1,:)-yy1*fixedP(3,:);yy3*fixedP(2,:)-yy2*fixedP(3,:)];

[~,~,pt3d2] = svd(coeff2);
pt3d2 = pt3d2(:,end);
p8pt2 = [pt3d2(1)/pt3d1(4) pt3d2(2)/pt3d1(4) pt3d2(3)/pt3d1(4)];

% Calculate distance
d3 = norm(p8pt1-p8pt2);
d3_true = d3*scalingfactor
```






