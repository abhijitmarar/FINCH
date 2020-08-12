%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CRLB calculation for a Gaussian PSF
% 01/20/2020
% Abhijit Marar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f_o = 3;                                % Focal length of objective (mm)    
NA = 1.42;                              % Numerical aperture of objective
wave = 670e-6;                          % Wavelength of light (mm)
n = 1.515;                              % Refractive Index
FWHM = (0.61*wave)/NA;                  % FWHM
sigma_o = FWHM/(2*sqrt(2*log(2)));      % standard deviation of Gaussian
z = 2.990:50e-6:3.010;                  % Distance between sample and objective
z_s = (z-3.000);
N = 6000;                               % No. of photons                                
%d = (2*wave)/NA^2;                      % Depth of focus (mm)
d = wave/(n*(1-(1-(NA/n)^2)^(1/2)));

sigma_g = zeros(1,length(z_s));
sigma_xy = zeros(1,length(z_s));
sigma_z = zeros(1,length(z_s));
for i = 1:length(z_s)
    sigma_g(i) = sigma_o*sqrt(1+(z_s(i)/d)^2);
    sigma_xy(i) = sigma_g(i)/sqrt(N);
    sigma_z(i) = (1/sqrt(N))*((d^2/(2*abs(z_s(i))))+(abs(z_s(i))/2));
end

defocus = z_s*1.e+3;
%{
figure
subplot(1,2,1)
plot(defocus,sigma_xy*1e+6,'LineWidth',3)
axis([-1 1 0 30])
title('CRLB_{xy}')
xlabel('Distance between sample and objective (\mum)')
ylabel('\sigma_{x}, \sigma_{y} (nm)')


subplot(1,2,2)
plot(defocus,sigma_z*1e+6,'LineWidth',3)
axis([-1 1 0 100])
title('CRLB_{z}')
xlabel('Distance between sample and objective (\mum)')
ylabel('\sigma_{z}(nm)')
%}

