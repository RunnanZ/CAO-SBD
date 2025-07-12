function [k]=k_cal(Blurred,kx,ky)
dirname='save'%
sig_noise=0.01;

%load blurred/sizeL
bmp_outname=sprintf('%s/diagfe_filt_sps_im%d_ker%d',dirname,i,j);
% 调用去卷积算法函数，传入加载的模糊图像y，核大小，真实图像x，噪声水平，和输出图像的文件名
[k]=k_calculate(Blurred,kx,ky,sig_noise,bmp_outname,1);%dispOn：0不显示处理中图像
end