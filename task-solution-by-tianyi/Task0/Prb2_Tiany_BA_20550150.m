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