function Result = Correlation(input1,input2)
%此处显示有关此函数的摘要
%两个输入做自相关运算
%input1---照明相干系数
%input2---孔径函数
FT2 = @(x) fftshift(fft2(x));
iFT2 = @(x) ifftshift(ifft2(ifftshift(x)));
F_input1 = FT2(input1);
F_conj_input2 = FT2(conj(input2));
Result1 = iFT2(F_input1.*F_conj_input2);

Result=Result1/max(max(Result1));
end
% 
% 
% qian=F_input1.*F_conj_input2;% X=fft(Aexp(iϕ))⋅fft(Aexp(−iϕ))   
% mul1=abs(F_input1).*abs(F_input1);% X=∣fft(Aexp(iϕ))∣ ^2
% % 求F_input1或者F_conj_input2​
% test1=sqrt(F_input1.*F_conj_input2);%sqrt(X)=∣fft(Aexp(iϕ))∣
% test2=
%  %逆推
%  qianfft=fftshift(fft2(fftshift(Result1)));
%     Result1_imag0 = complex( real(Result1), zeros(size(Result1)));
%     F_input1x2_imag0=fftshift(fft2(fftshift(Result1_imag0)));%对应F_input1.*F_conj_input2
%     F_input1x2 = complex( real(F_input1x2_imag0), zeros(size(F_input1x2_imag0)));
%  
% %F_input1.*F_conj_input2 虚部应该为0
% 
% F_input1
% 
% 
% IMAG1=imag(F_input1.*F_conj_input2);%Result1
% IMAG1_1=max(IMAG1,1e-6);%max=1.9e-6
% max_value = max(IMAG1_1(:));
% 
% % [max_value, max_index] = max(IMAG1_1(:));
% % [row, col] = ind2sub(size(IMAG1_1), max_index);