### Problem 1
##### Exercise 001

```matlab
for i = 1:10
    for j = 1:10
    M3(i,j) = 1;
    end
end

M3;

for k = 1:5
    for l = 1:5
        if k == l
            M5(k,l)=1;
        else
            M5(k,l)=0;
        end
    end
end

M5;
```
##### Exercise 002
```matlab
V = [1, 1, 1, 1, 1];
E2 = diag(V,1);
```
##### Exercise 003
```matlab
M1 = reshape(1:64, 8, 8);

M2 = transpose(M1);

M3 = reshape(1:64, 8, 8);
M3(1:4,1:4) = 1;

M4 = reshape(1:64, 8, 8);
M4(4:5, 4:5) = 0;

M5 = reshape(1:64, 8, 8);
M_5 = repelem(eye(4, 4), 2, 2);
M5(M_5 == 1) = 1;

M6 = reshape(1:64, 8, 8);
M_6 = rot90(M_5);
M6(M_6 == 1) = 1;

M7 = reshape(1:64, 8, 8);
M_7 = repelem(eye(4, 4), 2, 2);
M7(M_7 == 1) = 0;

M8 = reshape(1:64, 8, 8);
M_8 = repelem(eye(4, 4), 2, 2);
M8(M_8 == 0) = 100

M9 = reshape(1:64, 8, 8);
M_9 = repelem(eye(4, 4), 2, 2);
M9(M_9 == 0) = 0
```
### Problem 2
##### Exercise 001
```matlab
clc, clear, close all
img = imread('coloredChips.png');

[w,h,c] = size(img); % w: width, h: height, c: channel
fprintf('width: %d, height: %d, and channel(depth): %d\n ', w, h, c)
px1 = [150 0 0];
px2 = [255 100 100];

loc1 = [];
count = 1;
for ii=1:size(img,1)
    for jj=1:size(img,2)
        if img(ii,jj,1) >= px1(1) && img(ii,jj,1) <= px2(1) ...
           && img(ii,jj,2) >= px1(2) && img(ii,jj,2) <= px2(2) ...
           && img(ii,jj,3) >= px1(3) && img(ii,jj,3) <= px2(3)
                    loc1(count,:)=[ii,jj];
                    count = count + 1;
                    img(ii,jj,:) = 0;
        end
    end
end

imshow(img)
```
### Problem 3
##### Exercise 001
```matlab
clc, clear all, close all

img = 255*ones(200, 200, 'uint8');

img(1:100, 101:200, :) = zeros(100, 100, 'uint8');

img(101:200, 1:100, :) = img(1:100, 101:200, :);

im = repmat(img, 4);

imshow(im);

rectangle('position', [1 1 799 799])
```
### Problem 4
##### Exercise 002
```matlab
clc, clear all, close all

sampImgDir = 'sample_images'; % directory of sample images
filenameImg = fullfile(sampImgDir, 'color_patch.png'); % image filename
img = imread(filenameImg); % read an image

subplot(611); imshow(img)

imSize = [size(img,1) size(img,2)];

% Retrive the 4 squares in the first column

img1 = img(1:50, 1:50, :);
img2 = img(51:100, 1:50, :);
img3 = img(101:150, 1:50, :);
img4 = img(151:200, 1:50, :);
hold on
subplot(612); imshow(img1); 
subplot(613); imshow(img2); 
subplot(614); imshow(img3); 
subplot(615); imshow(img4);
hold off

rect = [101, 101, 50, 50];

% Get pixel values for blue color

p = impixel(img3, 1, 1);

% Replace the rectangle by blue color
for i = rect(1):rect(1)+rect(3)
    for j = rect(2):rect(2)+rect(4)
        img(i, j, 1) = p(1);
        img(i, j, 2) = p(2);
        img(i, j, 3) = p(3);
    end
end

subplot(616); imshow(img)
```
### Problem 5
##### Exercise 003
```matlab
clc, clear all, close all

sampImgDir = 'sample_images'; % directory of sample images

% Read images

filenameImg = fullfile(sampImgDir, 'cameraman.tif'); 
img = imread(filenameImg); 

filenameImg2 = fullfile(sampImgDir, 'llama.jpg'); 
imgg = imread(filenameImg2); 

% Rotate and flip images

img90 = imrotate(img, 90);
img180 = imrotate(img, 180);
img270 = imrotate(img, 270);
imglr = fliplr(img);
imgud = flipud(img);

imgg90 = imrotate(imgg, 90);
imgg180 = imrotate(imgg, 180);
imgg270 = imrotate(imgg, 270);
imgglr = fliplr(imgg);
imggud = flipud(imgg);

% Plot images

imshow(img)
imshow(img90)
imshow(img180)
imshow(img270)
imshow(imglr)
imshow(imgud)

imshow(imgg)
imshow(imgg90)
imshow(imgg180)
imshow(imgg270)
imshow(imgglr)
imshow(imggud)
```
### Problem 6
**Q.003:** Please explain why c is changed. 

Because the value of c indicates how many color channels there are in the image. Red and blue are added to the grey scale image, so the number of color channels changed. That's why c is changed.

**Q.004:** Please explain the difference between above two different methods for inserting circles. 

The first method with plot function places the circles on top of the true color images while the second mathod with insertShape function overwrites the pixel values.

**Q.006:** Please compare these results. Why does the third one produce "False"?

Because 0 is black and 255 is white, 255 is the maximum on a greyscale. grey300 would still be a matrix with 255, so grey300(10, 10) is not equal to 300.

**Q.007:** What is byte or bit? 

A 'bit' is the smallest unit of storage and can be 0 or 1. A byte is a collection of 8 bits and thus can represent 2^8 = 256 values.

**Q.010:** Are high resolution images always high quality? Do high-quality photos have high resolution? 

High resolution images do not always have high quality. If the image is grainy and pixelated, higher dpi doesn not help. High quality photos usually have high resolutions for them ot be clear, vivid and detailed.

**Q.011:** What did you learn from this example? 

Images are composed of different colors. Each color channel represents a color and the image can be split into different color channels. They can be merged to obtain the original image as well.

**Q.012:** What is loseless image compresssion? 

It is a way of reducing the size of an image while maintaining the quality.

**Q.013:** With this result, please estimate the filesize of a 8,688 x 5,792 image stored as tif (lossless compression)?

8688 x 5792 = 50320896 bit = 6290112 byte = 6142.69 KB = 6.00 MB

**Q.016:** What is the difference between im2double and double? 

The function double() only converts the specified array to a variable of type double, keeping the same values. The function im2double() also rescales the image to [0, 1], which is convenient when working other image processing functions.


