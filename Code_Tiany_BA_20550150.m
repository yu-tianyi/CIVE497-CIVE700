%% Problem 1

clc; clear all; close all;

% (a)

denom = (0-3)*(-2-8)-(4-0)*(8-10);
numerx = (0*0-4*3)*(8-10)-(0-3)*(8*8-(-2)*10);
numery = (0*0-4*3)*((-2)-8)-(4-0)*(8*8-(-2)*10);

px = numerx/denom; 
py = numery/denom;
fprintf('(px,py) = (%1.2f, %1.2f) \n', px, py);


% (b)
l1 = cross([0,4,1]',[3,0,1]');
l2 = cross([8,-2,1]',[10,8,1]');
x = cross(l1,l2);
px = x(1)/x(3);
py = x(2)/x(3);
fprintf('(px,py) = (%1.2f, %1.2f) \n', px, py);
%% Problem 2
%%
% (a)

clc; clear all; close all;

a = 1;
b = 0;
c = 1;
d = 4;
e = 2;
f = -29;

conic = [a b/2 d/2; b/2 c e/2; d/2 e/2 f];
m = 1;
n = 4;
lmn = conic*[m;n;1];

% alpha*x + beta*y + 1 = 0
alpha = lmn(1)/lmn(3);
beta = lmn(2)/lmn(3);
fprintf('line euqation: %fx+%fy+1=0 \n', alpha, beta);

% y-axis: x=0
lya  = [1 0 0];

x = cross(lmn, lya);
px = x(1)/x(3);
py = x(2)/x(3);
fprintf('(px,py) = (%1.2f, %1.2f) \n', px, py);

% (b)

m2 = 3;
n2 = -4;
lmn2 = conic*[m2;n2;1];

x2 = cross(lmn, lmn2);
px = x2(1)/x2(3);
py = x2(2)/x2(3);
fprintf('(px,py) = (%1.2f, %1.2f) \n', px, py);

%% Problem 3
%%
clc; clear all; close all;

I = imread('20190206_155248.jpg');
I=rgb2gray(I);

points = detectHarrisFeatures(I);

strongest = points.selectStrongest(10);
imshow(I); hold on;
plot(strongest); hold off;

ss = strongest.Location

x = [2030.2 1542.7 2237.7 2794.7];
y = [638.4 1264.8 1760.2 1000.8];

% First vanishing point
l1 = cross([2030.2,638.4,1]',[2794.7,1000.8,1]');
l2 = cross([1542.7,1264.8,1]',[2237.7,1760.2,1]');
x1 = cross(l1,l2);
x1coordinate = [x1(1)/x1(3) x1(2)/x1(3)]

% Second vanishing point
l3 = cross([2030.2,638.4,1]',[1542.7,1264.8,1]');
l4 = cross([2794.7,1000.8,1]',[2237.7,1760.2,1]');
x2 = cross(l3,l4);
x2coordinate = [x2(1)/x2(3) x2(2)/x2(3)]


% Vanishing line
l5 = cross(x1,x2)

%% Problem 4
%%
clc; clear all; close all;

x = [3 7 8 7];
y = [1 4 2 -1];

% First vanishing point
l1 = cross([3,1,1]',[7,4,1]');
l2 = cross([8,2,1]',[7,-1,1]');
x1 = cross(l1,l2);

% Second vanishing point
l3 = cross([3,1,1]',[7,-1,1]');
l4 = cross([8,2,1]',[7,4,1]');
x2 = cross(l3,l4);

% Vanishing line
l5 = cross(x1,x2);

h = [1,0,0; 0,1,0; l5(1),l5(2),l5(3)]

% Affine transformation
v_old = [x;y;1 1 1 1];
v_new = h*v_old;
newpoint1 = [v_new(1,1)/v_new(3,1) v_new(2,1)/v_new(3,1)]
newpoint2 = [v_new(1,2)/v_new(3,2) v_new(2,2)/v_new(3,2)]
newpoint3 = [v_new(1,3)/v_new(3,3) v_new(2,3)/v_new(3,3)]
newpoint4 = [v_new(1,4)/v_new(3,4) v_new(2,4)/v_new(3,4)]
%% Problem 5
%%
% (b)
clc; clear all; close all;

col1 = [1;2;1;1];
col2 = [2;2;-1;2];
col3 = [0;6;1;4];
col4 = [1;4;0;3];
col5 = [1;-1;1;2];

S = [col1 col2 col3 col4 col5];
Sr = rref(S);
basisofsubspace = [col1 col2 col3 col5]

% (c)

clc; clear all; close all;
% In linear algebra, the column space (also called the range or image) of a matrix A is the span (set of all possible linear combinations) of its column vectors. The column space of a matrix is the image or range of the corresponding matrix transformation.


col1 = [-1;2;-2];
col2 = [3;0;4];
col3 = [3;6;2];
col4 = [2;1;4];

A = sym([col1 col2 col3 col4]);
Acolumnspace = rref(A)
Arowspace = transpose(Acolumnspace)
Anullspace = null(A)


% (d)
clc; clear all; close all;
col1 = [1;-2;0;2];
col2 = [0;1;5;10];
col3 = [0;-3;-14;-28];
col4 = [0;-2;-9;-18];
col5 = [2;-4;0;4];

A = sym([col1 col2 col3 col4 col5]);
Anullspace = null(A)
nullityofA = 2

% (e)

clc; clear all; close all;
col1 = [1;-2];
col2 = [0;1];
col3 = [0;-3];
col4 = [0;-2];
col5 = [2;-4];

A = sym([col1 col2 col3 col4 col5]);
Anullspace = null(A)
nullityofA = 3


%% Problem 6
%%
clear; close all; clc; format shortG;
 
%% Parameter
imgBoardFile = 'IMG_0069.JPG';
imgPicFile = 'adrian-trinkaus-681192-unsplash.jpg';
 
info = imfinfo(imgPicFile);
sizePic = [info.Width info.Height]
 
%% Step1: Pick four corners of your white board in (a)
imgBoard = imread(imgBoardFile);
figure(1); imshow(imgBoard);
p = drawpolygon('LineWidth',7,'Color','black');
corner = p.Position
 
%% Step2: Compute H (Your Section)
a1 = corner(1,1);
a2 = corner(2,1);
a3 = corner(3,1);
a4 = corner(4,1);
b1 = corner(1,2);
b2 = corner(2,2);
b3 = corner(3,2);
b4 = corner(4,2);
 
x1 = 0;
x2 = 0;
x3 = sizePic(1);
x4 = sizePic(1);
y1 = 0;
y2 = sizePic(2);
y3 = sizePic(2);
y4 = 0;
 
 
A = sym([x1 y1 1 0 0 0 -x1*a1 -y1*a1 -a1;...
    0 0 0 x1 y1 1 -x1*b1 -y1*b1 -b1;...
    x2 y2 1 0 0 0 -x2*a2 -y2*a2 -a2;...
    0 0 0 x2 y2 1 -x2*b2 -y2*b2 -b2;...
    x3 y3 1 0 0 0 -x3*a3 -y3*a3 -a3;...
    0 0 0 x3 y3 1 -x3*b3 -y3*b3 -b3;...
    x4 y4 1 0 0 0 -x4*a4 -y4*a4 -a4;...
    0 0 0 x4 y4 1 -x4*b4 -y4*b4 -b4]);
 
HH = null(A);
HH = double(HH);
H = [HH(1) HH(2) HH(3);HH(4) HH(5) HH(6);HH(7) HH(8) HH(9)]'
 
 
%% Step3: Overlay your picture (I think there is a better way to do this)
imgPic = imread(imgPicFile);
[imgPicTran, RB] = imwarp(imgPic, projective2d(H));
BWPic = roipoly(imgPicTran, corner(:,1)-RB.XWorldLimits(1), corner(:,2)-RB.YWorldLimits(1));
 
BWBoard = ~roipoly(imgBoard, corner(:,1), corner(:,2));
RA = imref2d(size(BWBoard));
 
imgBoardMask = bsxfun(@times, imgBoard, cast(BWBoard, 'like', imgBoard));
imgPicTranMask = bsxfun(@times, imgPicTran, cast(BWPic, 'like', imgPicTran));
 
imgFinal(:,:,1) = imfuse(imgBoardMask(:,:,1),RA, imgPicTranMask(:,:,1),RB,'diff');
imgFinal(:,:,2) = imfuse(imgBoardMask(:,:,2),RA, imgPicTranMask(:,:,2),RB,'diff');
imgFinal(:,:,3) = imfuse(imgBoardMask(:,:,3),RA, imgPicTranMask(:,:,3),RB,'diff');
 
imshow(imgFinal); imwrite(imgFinal, 'result.jpg');
%% Problem 7
%%
clear; close all; clc; format shortG;
 
%% Parameter
imgBoardFile = 'IMG_3331.JPG';
 
 
%% Step1: Pick four corners of your white board in (a)
imgBoard = imread(imgBoardFile);
figure(1); imshow(imgBoard);
p = drawpolygon('LineWidth',7,'Color','black');
corner = p.Position;
 
%% Step2: Compute H (Your Section)
a1 = corner(1,1);
a2 = corner(2,1);
a3 = corner(3,1);
a4 = corner(4,1);
b1 = corner(1,2);
b2 = corner(2,2);
b3 = corner(3,2);
b4 = corner(4,2);
 
x1 = 0;
x2 = 0;
x3 = 11;
x4 = 11;
y1 = 0;
y2 = 8.5;
y3 = 8.5;
y4 = 0;
 
 
A = sym([x1 y1 1 0 0 0 -x1*a1 -y1*a1 -a1;...
    0 0 0 x1 y1 1 -x1*b1 -y1*b1 -b1;...
    x2 y2 1 0 0 0 -x2*a2 -y2*a2 -a2;...
    0 0 0 x2 y2 1 -x2*b2 -y2*b2 -b2;...
    x3 y3 1 0 0 0 -x3*a3 -y3*a3 -a3;...
    0 0 0 x3 y3 1 -x3*b3 -y3*b3 -b3;...
    x4 y4 1 0 0 0 -x4*a4 -y4*a4 -a4;...
    0 0 0 x4 y4 1 -x4*b4 -y4*b4 -b4]);
 
HH = null(A);
HH = double(HH);
H = [HH(1) HH(2) HH(3);HH(4) HH(5) HH(6);HH(7) HH(8) HH(9)]
 
figure(2); imshow(imgBoard);
q = drawpolyline('LineWidth',7,'Color','black');
 
corner2 = q.Position;
 
corner2 = [corner2(1,1) corner2(2,1); corner2(1,2) corner2(2,2); 1 1];
 
corner2 = inv(H)*corner2;
 
point1 = [corner2(1,1)/corner2(3,1);corner2(2,1)/corner2(3,1)];
 
point2 = [corner2(1,2)/corner2(3,2);corner2(2,2)/corner2(3,2)];
 
measurement = norm(point1 - point2)