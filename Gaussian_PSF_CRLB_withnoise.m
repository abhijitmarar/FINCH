%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CRLB calculation for a Gaussian PSF with noise
% 02/07/2020
% Abhijit Marar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clear all

f_o = 3;                                % Focal length of objective (mm)    
NA = 1.42;                              % Numerical aperture of objective
wave = 670e-6;                          % Wavelength of light (mm)
n = 1.515;                              % Refractive index
FWHM = (0.61*wave)/NA;                  % FWHM
sigma_o = FWHM/(2*sqrt(2*log(2)));      % standard deviation of Gaussian
z = 2.990:50e-6:3.010;                  % Distance between sample and objective
z_s = (z-3.000);
N = 6000;                               % No. of photons                                
%d = (2*wave)/NA^2;                     % Depth of focus (mm)
d = wave/(n*(1-sqrt(1-(NA^2/n^2))));
bg = 1000;                              % No. of bg photons/area    


sigma_g = zeros(1,length(z_s));
fisher_x = zeros(1,length(z_s));
fisher_z = zeros(1,length(z_s));
sigma_xy = zeros(1,length(z_s));
sigma_z = zeros(1,length(z_s));
dzg_dz = zeros(1,length(z_s));

for i = 1:length(z_s)
    sigma_g(i) = sigma_o*sqrt(1+(z_s(i)/d)^2);
    q_gauss = @(x,y) (1./(2*pi*sigma_g(i).^2)).*exp(-(x.^2+y.^2)/(2*sigma_g(i).^2));
    dq_dx = @(x,y) (-x./(2.*pi.*sigma_g(i).^4)).*exp(-(x.^2+y.^2)/(2.*sigma_g(i).^2));
    dq_dzg = @(x,y) q_gauss(x,y).*((x.^2+y.^2-2.*sigma_g(i).^2)./(sigma_g(i).^3));
    dzg_dz(i) = ((sigma_o*z_s(i))/d^2)*(1+(z_s(i)/d)^2)^(-1/2); 
    F_xx =  @(x,y) ((N./((bg/N)+q_gauss(x,y))).*dq_dx(x,y).*dq_dx(x,y));
    F_zz = @(x,y) ((N./((bg/N)+q_gauss(x,y))).*dq_dzg(x,y).*dq_dzg(x,y).*dzg_dz(i).*dzg_dz(i));
    fisher_x(i) = integral2(F_xx,-6*sigma_g(i),6*sigma_g(i),-6*sigma_g(i),6*sigma_g(i));
    fisher_z(i) = integral2(F_zz,-6*sigma_g(i),6*sigma_g(i),-6*sigma_g(i),6*sigma_g(i));
    sigma_xy(i) = 1/sqrt(fisher_x(i));
    sigma_z(i) = 1/sqrt(fisher_z(i));
end


defocus = z_s*1.e+3;
%{
figure
subplot(1,2,1)
plot(defocus,CRLB_x.*1e+6,'LineWidth',3)
axis([-1 1 0 30])
title('CRLB_{xy}')
xlabel('Distance between sample and objective (\mum)')
ylabel('\sigma_{x}, \sigma_{y} (nm)')

subplot(1,2,2)
plot(defocus,CRLB_z*1e+6,'LineWidth',3)
axis([-1 1 0 100])
title('CRLB_{z}')
xlabel('Distance between sample and objective (\mum)')
ylabel('\sigma_{z}(nm)')

%}
