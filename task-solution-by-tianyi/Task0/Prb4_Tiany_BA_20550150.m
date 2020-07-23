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