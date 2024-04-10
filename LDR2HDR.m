clear
close all
clc

%% Importdata
resolution = '_low';
FilePath = 'Data\';
for i = 1:4
    I(:,:,:,i) = double(imread([FilePath,num2str(i),resolution,'.png']));
end
width = size(I,1);
height = size(I,2);
figure
subplot(2,4,1)
imshow(uint8(I(:,:,:,1)),[])
subplot(2,4,2)
imshow(uint8(I(:,:,:,2)),[])
subplot(2,4,3)
imshow(uint8(I(:,:,:,3)),[])
subplot(2,4,4)
imshow(uint8(I(:,:,:,4)),[])
subplot(2,4,5)
histogram(I(:,:,:,1))
xlim([0,255])
subplot(2,4,6)
histogram(I(:,:,:,2))
xlim([0,255])
subplot(2,4,7)
histogram(I(:,:,:,3))
xlim([0,255])
subplot(2,4,8)
histogram(I(:,:,:,4))
xlim([0,255])

data = importdata('E:\Onedrive\PhDCourses\EC520\Project\Data\ps_hdr.tiff');


%% Splitting the image into chunks
xsize = size(I,1);
ysize = size(I,2);
num = size(I,4);
xchunk = 10;
ychunk = 15;
if mod(xsize,xchunk) ~=0 || mod(ysize,ychunk) ~=0
    error('Size of chunk is not an integer');
end
x_part_size = xsize/xchunk;
y_part_size = ysize/ychunk;
I_part = zeros(x_part_size, y_part_size, 3, num, xchunk, ychunk);
I_part_g = zeros(x_part_size, y_part_size, 3, num, xchunk, ychunk);

%% gradient
G1 = [-3 0 3; -10 0 10; -3 0 3];
G2 = [-3 -10 -3; 0 0 0; 3 10 3];
G1 = [-1 0 1; -2 0 2; -1 0 1];
G2 = [-1 -2 -1; 0 0 0; 1 2 1];
I_g1 = zeros(size(I));
I_g2 = zeros(size(I));
for m = 1:num
    for n = 1:3
        I_g1(:,:,n,m) = abs(conv2(I(:,:,n,m),G1,'same'));
        I_g2(:,:,n,m) = abs(conv2(I(:,:,n,m),G2,'same'));
    end
end
I_g = sqrt(I_g1.^2+I_g2.^2);
figure
subplot(2,2,1)
imshow(sum(I_g(:,:,:,1),3),[])
subplot(2,2,2)
imshow(sum(I_g(:,:,:,2),3),[])
subplot(2,2,3)
imshow(sum(I_g(:,:,:,3),3),[])
subplot(2,2,4)
imshow(sum(I_g(:,:,:,4),3),[])
% figure
% subplot(2,2,1)
% imshow(sum(I_g1(:,:,:,1),3),[])
% subplot(2,2,2)
% imshow(sum(I_g1(:,:,:,2),3),[])
% subplot(2,2,3)
% imshow(sum(I_g1(:,:,:,3),3),[])
% subplot(2,2,4)
% imshow(sum(I_g1(:,:,:,4),3),[])
% figure
% subplot(2,2,1)
% imshow(sum(I_g2(:,:,:,1),3),[])
% subplot(2,2,2)
% imshow(sum(I_g2(:,:,:,2),3),[])
% subplot(2,2,3)
% imshow(sum(I_g2(:,:,:,3),3),[])
% subplot(2,2,4)
% imshow(sum(I_g2(:,:,:,4),3),[])


exposureTimes = [1/400, 1/800, 1/2000, 1/5000];

% shape of I_part: x, y, color, number of image, number of x, number of y
for m = 1:xchunk
    for n = 1:ychunk
        I_part(:, :, :, :, m, n) = I((m-1)*x_part_size+1:m*x_part_size,(n-1)*y_part_size+1:n*y_part_size, :, :);
        I_part_g(:, :, :, :, m, n) = I_g((m-1)*x_part_size+1:m*x_part_size,(n-1)*y_part_size+1:n*y_part_size, :, :);
    end
end
I_part_g_sum = sum(I_part_g, [1,2,3]);
I_part_g_reshape = reshape(I_part_g_sum,[num,xchunk,ychunk]);



for m = 1:xchunk
    for n = 1:ychunk
        [maxg(m,n), pos(m,n)] = max(I_part_g_sum(:,:,:,:,m,n));
        I_pick((m-1)*x_part_size+1:m*x_part_size,(n-1)*y_part_size+1:n*y_part_size, :, :) = I_part(:, :, :, pos(m,n), m, n);
        I_pick_g((m-1)*x_part_size+1:m*x_part_size,(n-1)*y_part_size+1:n*y_part_size, :, :) = I_part_g(:, :, :, pos(m,n), m, n);
    end
end
I_pick = I_pick/max(I_pick(:))*255;
figure
imshow(uint8(I_pick))

figure
imshow(sum(I_pick_g,3),[])




midexposure = mean(exposureTimes);
weight = exposureTimes;
HDR = weight(1)*I(:,:,:,1)+weight(2)*I(:,:,:,2)+weight(3)*I(:,:,:,3)+weight(4)*I(:,:,:,4);

HDR = HDR/max(HDR(:))*255;
figure
imshow(uint8(HDR))
























