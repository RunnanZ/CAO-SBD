%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main file for experimental data for CAO-SBD.
% The Code is created based on the method described in the following paper
% [1]. Zhang R, Du H, Zhou N, et al. Computational Adaptive Optics for Fluorescence Microscopy via Sparse Blind Deconvolution[J].
%      Laser & Photonics Reviews, 2025: 2500032. https://doi.org/10.1002/lpor.202500032
%
% [2]. Levin A, Weiss Y, Durand F, et al. Understanding and evaluating blind deconvolution algorithms[C]//
%      2009 IEEE conference on computer vision and pattern recognition. IEEE, 2009: 1964-1971.&#x20;

%    The Code is modified and extended from Levin's code
%    Contact: Runnan Zhang (runnanzhang@njust.edu.cn) or Chao Zuo (zuochao@njust.edu.cn)
%    Date  : 2025/7/4

clc
clear
close all

load('BPAE_0001_C002T001.mat'); % load confocal image
figure,imshow(I,[])


%% Setting up system parameters
[Ny, Nx] = size(I);
NA = 1.45;                                                                 % NA of objective lens
Lambda = 0.461;                                                            % Wavelength
k = 2*pi/Lambda;                                                           % Wave vector
Pixelsize = 0.031;                                                         % Pixel size
Max_frequency = NA/Lambda;                                                 % coherent resolution limit
delta_x = 1/(Pixelsize*Nx);                                                % frequency sampling X.
delta_y = 1/(Pixelsize*Ny);                                                % frequency sampling Y.
fx = (-fix(Nx/2):1:fix((Nx-1)/2))*delta_x;                                 % frequency coordinate fx.
fy = (-fix(Ny/2):1:fix((Ny-1)/2))*delta_y;                                 % frequency coordinate fy.
xx = (-fix(Nx/2):1:fix((Nx-1)/2))/delta_x;                                 % spatial coordinate X.
yy = (-fix(Ny/2):1:fix((Ny-1)/2))/delta_y;                                 % spatial coordinate Y.

[fx2D, fy2D] = meshgrid(fx,fy);
[x2D, y2D] = meshgrid(xx,yy);


%% Calculate transfer functions

Pupil = zeros(Ny,Nx);
Pupil((fx2D.^2+fy2D.^2)<=Max_frequency.^2) = 1;

OTF = Correlation(Pupil,Pupil); % incoherence OTF
figure;
imshow(abs(OTF),[]);title('OTF Function')
% PSF = abs(ifftshift(ifft2(ifftshift(OTF)))); % wide field PSF


OTF = Correlation(OTF,OTF);% confocal OTF
figure;
imshow(abs(OTF),[]);title('confococal OTF Function')
PSF = abs(ifftshift(ifft2(ifftshift(OTF)))); % confocal PSF

figure;
imshow(PSF,[]);title('PSF calculate')
Mask = zeros(size(OTF));
Mask(abs(OTF)>0.0001)=1;
Mask = fftshift(Mask);

Blurred = mat2gray(I);

% Sets the estimate kernel size for blind PSF estimation
kx = 15;
ky = 15;

% calculate k
k = k_cal(Blurred, kx, ky);



%% Zernike fitting
H_k = psf2otf(k,sizeI);
PSF1 = otf2psf(H_k);
PSF_reshape = reshape(PSF1,[Ny*Nx,1])/10; %  PSF_reshape   k
dim = [Nx, Ny];
zernike_poly        = genZernikePoly(fx2D, fy2D, NA, Lambda, 21); % zernike poly
% zernike_coeff_k     =[0;1;0;0;0;...
%     0;0;0;0;0;...
%     0;0;0;0;0;...
%     0;0;0];
% Q0 = [1;1;1];
% Q0 =[1;1;1;0;0;...
%      0;0];
n = 7;
for i = 1:n
    X(:,i) = zernike_poly(:,i);
end
Q0 = zeros(n, 1);
for i = 1:3
    Q0(i,1) = 1;
end
Pupil_reshape = Pupil(:);
% PSF2 = fftshift(PSF2);
% PSF_reshape = PSF2(:);
% AAAAA = reshape(PSF_reshape,[512 512]);
% figure,imshow(AAAAA,[])
Fun = @(Q,X)...
    reshape(abs(ifftshift(ifft2(reshape(Pupil_reshape.*exp(1i.*( X * Q)),dim)))).^2,[Ny*Nx,1]);


PSF_Q0 = Fun(Q0,X);
PSF_correction = reshape(PSF_Q0,dim);
figure,
subplot(1,2,1),imshow(PSF_correction,[]);
subplot(1,2,2),imshow(k,[])

[a] = nlinfit(X,PSF_reshape,Fun,Q0);
PSF_Q0 = Fun(a,X);
PSF_correction = reshape(PSF_Q0,dim);
% figure,
% subplot(1,2,1),imshow(PSF_correction,[]);
% subplot(1,2,2),imshow(PSF,[]);

%% blind deconvolution

niter = 30;

ImageEstimate = Blurred;
psfEstimate = PSF_correction; % Convolution kernel selection: PSF, k, PSF_correction
% figure,imshow(psfEstimate,[]);
sizeh = size(psfEstimate);
ImageEstimateFlip = flipPSF(ImageEstimate);
psfFlip = flipPSF(psfEstimate);

smallValue = 0;
figure,
for i = 1:niter
    disp(i)

    HI = psf2otf(ImageEstimate);
    Hpsf = psf2otf(psfEstimate,sizeI);
    Conv = otf2psf(HI.*Hpsf);
    DV = Blurred./Conv;
    DV_Conv = otf2psf(psf2otf(DV).*psf2otf(ImageEstimateFlip),sizeh);
    psfEstimate = DV_Conv.*psfEstimate;

    H_psfEstimate = psf2otf(psfEstimate,sizeI).*Mask;
    psfEstimate = otf2psf(H_psfEstimate,sizeh);

    psfEstimate = max(psfEstimate,smallValue);     % non-negative constraint
    sumPSF = sum(psfEstimate(:));
    psfEstimate = psfEstimate./(sumPSF + (sumPSF==0)*eps);

    psfFlip = flipPSF(psfEstimate);
    H_psfFlip = psf2otf(psfFlip,sizeI);


    HI = psf2otf(ImageEstimate);
    Hpsf = psf2otf(psfEstimate,sizeI);
    Conv = otf2psf(HI.*Hpsf);
    DV = Blurred./Conv;

    DV_Conv = otf2psf(psf2otf(DV).*H_psfFlip,sizeI);
    ImageEstimate = DV_Conv.*ImageEstimate;

    ImageEstimate = max(ImageEstimate,smallValue);     % non-negative constraint
    sumImage = sum(ImageEstimate(:));
    ImageEstimate = ImageEstimate./(sumImage + (sumImage==0)*eps);
    ImageEstimateFlip = flipPSF(ImageEstimate);


    % subplot(1,2,1),imshow(I,[]);
    % subplot(1,2,2),imshow(mat2gray(ImageEstimate),[]);

    
    figure(99),
    subplot(1,2,1),imshow(Blurred,[]);title('blurred')
    subplot(1,2,2),imshow(ImageEstimate,[]);title('CAO-SBD')% Reconstruction result
    colormap(hot)
    pause(0.01)
end
ImageEstimate = mat2gray(ImageEstimate);




%% traditional RL deconvblind
Blurred=mat2gray(I);
J = deconvlucy(Blurred,PSF,30);

figure,
subplot(1,3,1),imshow(Blurred,[]);title('blurred')
subplot(1,3,2),imshow(J,[]);title('RL deconvolution')
subplot(1,3,3),imshow(ImageEstimate,[]);title('CAO-SBD')
colormap(hot)