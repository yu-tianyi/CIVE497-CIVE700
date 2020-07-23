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