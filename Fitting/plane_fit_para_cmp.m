clear;clc;
load('plane_fit_para_m0.mat')
par0 = plane_fit_para;
load('plane_fit_para_m1.mat')
par1 = plane_fit_para;
test0 = squeeze(par0(2,:,:));
test1 = squeeze(par1(2,:,:));
a0 = mean(par0(1,:,:),'all');
a1 = mean(par1(1,:,:),'all');
b0 = mean(par0(2,:,:),'all');
b1 = mean(par1(2,:,:),'all');
c0 = mean(par0(3,:,:),'all');
c1 = mean(par1(3,:,:),'all');
figure
hold on
subplot(1,3,1)
imagesc((test0-test1)')
title('m0-m1')
subplot(1,3,2)
imagesc(test0')
title('m0')
subplot(1,3,3)
imagesc(test1')
title('m1')
