%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CRLB calculation for Astigmatic PSF without noise
% 02/07/2020
% Abhijit Marar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all

f_o = 3;                                % Focal length of objective (mm)    
NA = 1.42;                              % Numerical aperture of objective
wave = 670e-6;                          % Wavelength of light (mm)
FWHM = (0.61*wave)/NA;                  % FWHM
sigma_o = FWHM/(2*sqrt(2*log(2)));      % standard deviation of Gaussian
z = 2.999:20e-6:3.001;                  % Distance between sample and objective
z_s = (z-3.000);
N = 1000;                               % No. of photons                                
d = 800e-6;                             % Depth of focus (mm)
gamma = 400e-6;                         % Amount of astigmatism (mm)

sigma_x = zeros(1,length(z_s));
sigma_y = zeros(1,length(z_s));
sigma_z = zeros(1,length(z_s));
CRLB_x = zeros(1,length(z_s));
CRLB_y = zeros(1,length(z_s));

for i = 1:length(z_s)
    sigma_x(i) = sigma_o*sqrt(1+((z_s(i)+gamma)/d)^2);
    sigma_y(i) = sigma_o*sqrt(1+((z_s(i)-gamma)/d)^2);
    CRLB_x(i) = sigma_x(i)/sqrt(N);
    CRLB_y(i) = sigma_y(i)/sqrt(N);
end

defocus = z_s*1.e+3;

figure
subplot(1,2,1)
plot(defocus,CRLB_x*1e+6,'LineWidth',3)
hold on
plot(defocus,CRLB_y*1e+6,'LineWidth',3)

axis([-1 1 0 30])
title('CRLB_{xy}')
xlabel('Distance between sample and objective (\mum)')
ylabel('\sigma_{x}, \sigma_{y} (nm)')
    