%TAHAP 1
%Normalisasi gambar original
%resize gambar kedalam dimensi pixel 400x400
Image_rgb = imread('*******');
Image_rgb = imresize(Image_rgb, [400 400]);
Image_rgb = double(Image_rgb);

%partisi gambar dalam ruang warna rgb
Image_red = Image_rgb(:,:,1);
Image_green = Image_rgb(:,:,2);
Image_blue = Image_rgb(:,:,3);

%Normalisasi in action
[row,col] = size(Image_rgb(:,:,1));
for y = 1:row %-->numberof rows in image
for x = 1:col %-->number of columns in the image
Red = Image_red(y,x);
Green = Image_green(y,x);
Blue = Image_blue(y,x);
NormalizedRed = Red/sqrt(Red^2 + Green^2 + Blue^2);
NormalizedGreen = Green/sqrt(Red^2 + Green^2 + Blue^2);
NormalizedBlue = Blue/sqrt(Red^2 + Green^2 + Blue^2);

Image_red(y,x) = NormalizedRed;
Image_green(y,x) = NormalizedGreen;
Image_blue(y,x) = NormalizedBlue;
end
end

Image_rgb(:,:,1) = Image_red;
Image_rgb(:,:,2) = Image_green;
Image_rgb(:,:,3) = Image_blue;

Image_rgb = Image_rgb .* Image_rgb;
Image_rgb = Image_rgb .* Image_rgb;

%TAHAP 2
%cek hasil gambar normalisasi
figure; imshow(Image_rgb);

%TAHAP 3
%Segmentasi menggunakan clustering Kmean
%%imread_rgb adalah gambar hasil normalisasi
he = Image_rgb;
imshow(he), title('udang');
text(size(he,2),size(he,1)+15,...
'udang windu aceh', ...
'FontSize',7,'HorizontalAlignment','right');
cform = makecform('srgb2lab');
lab_he = applycform(he,cform);
ab = double(lab_he(:,:,2:3));
nrows = size(ab,1);
ncols = size(ab,2);
ab = reshape(ab,nrows*ncols,2);
nColors = 3;

% Pengulangan clustering
[cluster_idx cluster_center] = kmeans(ab,nColors,'distance','sqEuclidean', ...
'Replicates',3);
pixel_labels = reshape(cluster_idx,nrows,ncols);
imshow(pixel_labels,[]), title('image labeled by cluster index');
segmented_images = cell(1,3);
rgb_label = repmat(pixel_labels,[1 1 3]);
for k = 1:nColors
color = he;
color(rgb_label ~= k) = 0;
segmented_images{k} = color;
end

%TAHAP4
%cek gambar ketiga segmen satu satu yaaa!
imshow(segmented_images{1}), title('objects in cluster 1');

imshow(segmented_images{2}), title('objects in cluster 2');

imshow(segmented_images{3}), title('objects in cluster 3');

%Tahap5
%Kalkulasi intensitas RGB pada masing-masing cluster kemudian direratakan.
redChannel1 = segmented_images {1} (:,:, 1);
greenChannel1 = segmented_images {1} (:,:, 2);
blueChannel1 = segmented_images {1} (:,:, 3);
redMean1 = mean(mean(redChannel1));
greenMean1 = mean(mean(greenChannel1));
blueMean1 = mean(mean(blueChannel1));
redChannel2 = segmented_images {2} (:,:, 1);
greenChannel2 = segmented_images {2} (:,:, 2);
blueChannel2 = segmented_images {2} (:,:, 3);
redMean2 = mean(mean(redChannel2));
greenMean2 = mean(mean(greenChannel2));
blueMean2 = mean(mean(blueChannel2));
redChannel3 = segmented_images {3} (:,:, 1);
greenChannel3 = segmented_images {3} (:,:, 2);
blueChannel3 = segmented_images {3} (:,:, 3);
redMean3 = mean(mean(redChannel3));
greenMean3 = mean(mean(greenChannel3));
blueMean3 = mean(mean(blueChannel3));
Allredmean = (redMean1 + redMean2 + redMean3)/3
Allgreenmean = (greenMean1 + greenMean2 + greenMean3)/3
Allbluemean = (blueMean1 + blueMean2 + blueMean3)/3

%catat nilai Allredmean Allgreenmean dan Allbluemean

%Tahap6
%kalkulasi intensitas RGB pada gambar original ternormalisasi
redChannel = Image_rgb(:,:, 1);
greenChannel = Image_rgb(:,:, 2);
blueChannel = Image_rgb(:,:, 3);
redMean = mean(mean(redChannel))
greenMean = mean(mean(greenChannel))
blueMean = mean(mean(blueChannel))

%kalkulasi intensitas RGB gambar original yang belum ternormalisasi
Image_ori = imread('*******');
Image_ori = imresize(Image_ori, [400 400]);
I = Image_ori;
redChannelori = I (:,:, 1);
greenChannelori = I (:,:, 2);
blueChannelori = I (:,:, 3);
redMeanori = mean(mean(redChannelori))
greenMeanori = mean(mean(greenChannelori))
blueMeanori = mean(mean(blueChannelori))

%Tahap 7
% kalkulasi luasan area warna hasil segmentasi Kmean
% untuk mempermudah ubah nama objek segmen 1, 2, dan 3
SG1 = segmented_images {1};
SG2 = segmented_images {2};
SG3 = segmented_images {3};

% mengubah segmen 1 kedalam images black and white
% segmen 1
thresholdsg1 = graythresh (SG1);
BWSG1 = im2bw (SG1, thresholdsg1);
imshow (BWSG1)
% simpan gambar BWSG1
% hitung luas pixel BWSG1
bwarea (BWSG1)

%segmen 2
thresholdsg2 = graythresh (SG2);
BWSG2 = im2bw (SG2, thresholdsg2);
imshow (BWSG2)
% simpan gambar BWSG2
% hitung luas pixel BWSG2
bwarea (BWSG2)

%segmen 3
thresholdsg3 = graythresh (SG3);
BWSG3 = im2bw (SG3, thresholdsg3);
imshow (BWSG3)
% simpan gambar BWSG3
% hitung luas pixel BWSG3
bwarea (BWSG3)

% luas area full body
% gunakan gambar original yang di resize 400 x 400
Image_ori = imread('******');
Image_ori = imresize(Image_ori, [400 400]);
I = Image_ori;
thresholdori = graythresh (I);
BWI = im2bw (I, thresholdori);
imshow (BWI)
% simpan gambar BWI
% hitung luas pixel BWI
% luas area yang dihitung adalah luas pada warna putih = 1 karena yang hitam = 0
bwarea (BWI)

-------------------------
%penting
%jika ada kegagalan string konversi ke binary image maka perlu dilakukan manual konversi dengan treshold 0.1
BW = im2bw (RGB2, 0.1);
imshow (BW)
bwarea (BW)

%============================================Section 2=====

%TAHAP 1
%Normalisasi gambar original
%resize gambar kedalam dimensi pixel 400x400
Image_rgb = imread('***********');
Image_rgb = imresize(Image_rgb, [400 400]);
Image_rgb = double(Image_rgb);

%partisi gambar dalam ruang warna rgb
Image_red = Image_rgb(:,:,1);
Image_green = Image_rgb(:,:,2);
Image_blue = Image_rgb(:,:,3);

%Normalisasi in action
[row,col] = size(Image_rgb(:,:,1));
for y = 1:row %-->numberof rows in image
for x = 1:col %-->number of columns in the image
Red = Image_red(y,x);
Green = Image_green(y,x);
Blue = Image_blue(y,x);
NormalizedRed = Red/sqrt(Red^2 + Green^2 + Blue^2);
NormalizedGreen = Green/sqrt(Red^2 + Green^2 + Blue^2);
NormalizedBlue = Blue/sqrt(Red^2 + Green^2 + Blue^2);

Image_red(y,x) = NormalizedRed;
Image_green(y,x) = NormalizedGreen;
Image_blue(y,x) = NormalizedBlue;
end
end

Image_rgb(:,:,1) = Image_red;
Image_rgb(:,:,2) = Image_green;
Image_rgb(:,:,3) = Image_blue;

Image_rgb = Image_rgb .* Image_rgb;
Image_rgb = Image_rgb .* Image_rgb;

%Tahap 2
% analisa LAB gambar ternormalisasi
% masukan data
fabric = Image_rgb;
figure(1), imshow(fabric), title('udang');
%menghitung sampel warna pada objek
load regioncoordinates;
nColors = 6;
sample_regions = false([size(fabric,1) size(fabric,2) nColors]);
for count = 1:nColors
sample_regions(:,:,count) = roipoly(fabric,region_coordinates(:,1,count),...
region_coordinates(:,2,count));
end
imshow(sample_regions(:,:,2)),title('sample region for red');
%konversi gambar RGB ke Lab
cform = makecform('srgb2lab');
lab_fabric = applycform(fabric,cform);
%mengitung rerata nilai a dan b
a = lab_fabric(:,:,2);
b = lab_fabric(:,:,3);
color_markers = repmat(0, [nColors, 2]);
for count = 1:nColors
color_markers(count,1) = mean2(a(sample_regions(:,:,count)));
color_markers(count,2) = mean2(b(sample_regions(:,:,count)));
end
disp(sprintf('[%0.3f,%0.3f]',color_markers(2,1),color_markers(2,2)));
%klasifikasi masing-masing pixel dengan metode nearest neighbour
color_labels = 0:nColors-1;
a = double(a);
b = double(b);
distance = repmat(0,[size(a), nColors]);
for count = 1:nColors
distance(:,:,count) = ( (a - color_markers(count,1)).^2 + ...
(b - color_markers(count,2)).^2 ).^0.5;
end
[value, label] = min(distance,[],3);
label = color_labels(label);
clear value distance;
rgb_label = repmat(label,[1 1 3]);
segmented_images = repmat(uint8(0),[size(fabric), nColors]);
for count = 1:nColors
color = fabric;
color(rgb_label ~= color_labels(count)) = 0;
segmented_images(:,:,:,count) = color;
end
imshow(segmented_images(:,:,:,2)), title('red objects');
imshow(segmented_images(:,:,:,3)), title('green objects');
imshow(segmented_images(:,:,:,4)), title('purple objects');
imshow(segmented_images(:,:,:,5)), title('magenta objects');
imshow(segmented_images(:,:,:,6)), title('yellow objects');
purple = [119/255 73/255 152/255];
plot_labels = {'k', 'r', 'g', purple, 'm', 'y'};
figure
for count = 1:nColors
plot(a(label==count-1),b(label==count-1),'.','MarkerEdgeColor', ...
plot_labels{count}, 'MarkerFaceColor', plot_labels{count});
hold on;
end
title('Scatterplot of the segmented pixels in ''a*b*'' space');
xlabel('''a*'' values');
ylabel('''b*'' values');
mata = a(:);
matb = b(:);

%tahap 3
%analisa LAB gambar original
fabric = imread('*******');
figure(1), imshow(fabric), title('udang');
%menghitung sampel warna pada objek
load regioncoordinates;
nColors = 6;
sample_regions = false([size(fabric,1) size(fabric,2) nColors]);
for count = 1:nColors
sample_regions(:,:,count) = roipoly(fabric,region_coordinates(:,1,count),...
region_coordinates(:,2,count));
end
imshow(sample_regions(:,:,2)),title('sample region for red');
%konversi gambar RGB ke Lab
cform = makecform('srgb2lab');
lab_fabric = applycform(fabric,cform);
%mengitung rerata nilai a dan b
a = lab_fabric(:,:,2);
b = lab_fabric(:,:,3);
color_markers = repmat(0, [nColors, 2]);
for count = 1:nColors
color_markers(count,1) = mean2(a(sample_regions(:,:,count)));
color_markers(count,2) = mean2(b(sample_regions(:,:,count)));
end
disp(sprintf('[%0.3f,%0.3f]',color_markers(2,1),color_markers(2,2)));
%klasifikasi masing-masing pixel dengan metode nearest neighbour
color_labels = 0:nColors-1;
a = double(a);
b = double(b);
distance = repmat(0,[size(a), nColors]);
for count = 1:nColors
distance(:,:,count) = ( (a - color_markers(count,1)).^2 + ...
(b - color_markers(count,2)).^2 ).^0.5;
end
[value, label] = min(distance,[],3);
label = color_labels(label);
clear value distance;
rgb_label = repmat(label,[1 1 3]);
segmented_images = repmat(uint8(0),[size(fabric), nColors]);
for count = 1:nColors
color = fabric;
color(rgb_label ~= color_labels(count)) = 0;
segmented_images(:,:,:,count) = color;
end
imshow(segmented_images(:,:,:,2)), title('red objects');
imshow(segmented_images(:,:,:,3)), title('green objects');
imshow(segmented_images(:,:,:,4)), title('purple objects');
imshow(segmented_images(:,:,:,5)), title('magenta objects');
imshow(segmented_images(:,:,:,6)), title('yellow objects');
purple = [119/255 73/255 152/255];
plot_labels = {'k', 'r', 'g', purple, 'm', 'y'};
figure
for count = 1:nColors
plot(a(label==count-1),b(label==count-1),'.','MarkerEdgeColor', ...
plot_labels{count}, 'MarkerFaceColor', plot_labels{count});
hold on;
end
title('Scatterplot of the segmented pixels in ''a*b*'' space');
xlabel('''a*'' values');
ylabel('''b*'' values');




