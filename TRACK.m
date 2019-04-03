clear
clc
close all

BW = imread('FINCAN3.png');
I = imbinarize(BW);
I2 = rgb2gray(BW);
I2 = imbinarize(I2);
imshow(I2)
impixelinfo

r1 = 282;
c1 = 115;

contour = bwtraceboundary(I2,[c1 r1],'S');

ANSS = regionprops(I2,'Centroid');
centroids = cat(1, ANSS.Centroid);


hold on
plot(contour(:,2),contour(:,1),'g','LineWidth',2)
hold on
plot(centroids(:,1), centroids(:,2), 'b*')

