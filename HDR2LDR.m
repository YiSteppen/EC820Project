close all
clear 
clc

%% 
HDR = double(imread('Data\NY_hdr.tif'));

L = 1:65535;
L = L';
%%
LDR_linear = round(HDR/257);


%%


alpha = 0.00005;
L_nlr = L./(1+L);
L_alpha = alpha*L./(1+alpha*L);
L_alpha = L_alpha/max(L_alpha(:));
L_prime = alpha*HDR;
LDR_nlr = L_prime./(1+L_prime);
LDR_nlr = LDR_nlr/max(LDR_nlr(:))*255;




%%
L_white = 5;
LDR_low = alpha*HDR;
L_d = alpha*L.*(1+alpha*L./L_white^2)./(1+alpha*L);
L_d = L_d/max(L_d(:));
LDR_d = LDR_low.*(1+LDR_low./L_white^2)./(1+LDR_low);
LDR_d = LDR_d/max(LDR_d(:))*255;

figure
subplot(2,3,1)
imshow(uint8(LDR_linear),[])
subplot(2,3,2)
imshow(uint8(LDR_nlr),[])
subplot(2,3,3)
imshow(uint8(LDR_d),[])
subplot(2,3,4)
plot(L, L/65535)
subplot(2,3,5)
plot(L,L_alpha)
subplot(2,3,6)
plot(L, L_d)


L_d = zeros([65535,5]);
LDR_d = zeros([size(LDR_nlr),5]);
figure
L_white_range = [0.3,1,3,6,18];
for i = 1:5
    L_white = L_white_range(i);
    L_d(:,i) = alpha*L.*(1+alpha*L./L_white^2)./(1+alpha*L);
    L_d(:,i) = L_d(:,i)/max(L_d(:,i),[],'all');
    LDR_d(:,:,:,i) = LDR_low.*(1+LDR_low./L_white^2)./(1+LDR_low);
    LDR_d(:,:,:,i) = LDR_d(:,:,:,i)/max(LDR_d(:,:,:,i),[],'all')*255;
    subplot(2,5,i)
    imshow(uint8(LDR_d(:,:,:,i)),[])
    subplot(2,5,i+5)
    plot(L, L_d(:,i))
end

