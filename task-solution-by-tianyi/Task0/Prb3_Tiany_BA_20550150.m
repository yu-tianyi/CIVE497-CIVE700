clc, clear all, close all

img = 255*ones(200, 200, 'uint8');

img(1:100, 101:200, :) = zeros(100, 100, 'uint8');

img(101:200, 1:100, :) = img(1:100, 101:200, :);

im = repmat(img, 4);

imshow(im);

rectangle('position', [1 1 799 799])