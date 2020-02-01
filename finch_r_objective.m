% FINCH reconstruction
% all units are in mm
theta1 = 2*pi/3;
z1 = 3.000;
z2 = 3.005;
y2 = 0.005;
%Np1 = 2500; Np2 = 1000;
I1 = finchH_objective(0,0,z1,theta1) + finchH_objective(0.001,0,z1,theta1) + finchH_objective(-0.001,0,z1,theta1);
theta2 = 2*pi/3;
I2 = finchH_objective(0,0,z1,theta2) + finchH_objective(0.001,0,z1,theta2) + finchH_objective(-0.001,0,z1,theta2);
theta3 = 4*pi/3;
I3 = finchH_objective(0,0,z1,theta3) + finchH_objective(0.001,0,z1,theta3) + finchH_objective(-0.001,0,z1,theta3);
%{
figure();
imagesc(I3);
figure();
imagesc(I3_moved);
I_complete = I3+I3_moved;
figure();
imagesc(I_complete);
%}


% add noise

Nph = 10000;

I1 = (Nph/sum(sum(I1)))*I1;
I2 = (Nph/sum(sum(I2)))*I2;
I3 = (Nph/sum(sum(I3)))*I3;

I1 = poissrnd(I1);
I2 = poissrnd(I2);
I3 = poissrnd(I3);

% assemble

IF = 1/(3*sqrt(3))*(I1*(exp(1j*theta3)-exp(1j*theta2)) + ...
    I2*(exp(1j*theta1)-exp(1j*theta3)) + ...
    I3*(exp(1j*theta2)-exp(1j*theta1)));

figure();
imagesc(abs(IF));
figure();
imagesc(angle(IF));

wave = 0.0005;
%z = -111.3419;
%z= -306.0482;
z = 300.0;
%z = ones(512,512);
%z  = -400*z;


N = 2*256;
x = linspace(-3,3,N);
[xx, yy] = meshgrid(x);
%xx = xx*(400/3); yy = yy*(400/3);
%{
Huygens Convolution
g_const = 1/(1j*wave);
g_den = sqrt(xx.^2 + yy.^2 + z.^2);
k = (2*pi)/wave;
g_num = exp(1j*k*g_den);
g = g_const*(g_num./g_den);
ft_g = fft2(fftshift(g));
ft_ccd = fft2(fftshift(IF));
ft_ccd2 = fft2(ft_ccd);
E_huygens = (fftshift(ifft2(ft_ccd2.*ft_g)));
%}
%Frensel Transform
s = fftshift(ifft2(fft2(fftshift(IF)).*fft2(fftshift(exp((1j*pi/wave/z)*(xx.^2 + yy.^2))))));
%Guide star reconstruction
%s = fftshift(ifft2(fft2(fftshift(IF)).*fft2(fftshift(IF))));

figure();
%imagesc(abs(E_huygens));
imagesc(abs(s));
colorbar;
