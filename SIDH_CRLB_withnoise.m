%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CRLB calculation for SIDH with Poisson noise
% 02/12/2020
% Abhijit Marar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DEFINE CONSTANTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f_o = 3;                            % Focal length of objective (mm)    
NA = 1.42;                          % Numerical aperture of objective
D_bpp = (2*f_o*NA);                 % Diameter of back pupil plane
wave = 515e-6;                      % Wavelength of light (mm)
k = 2*pi/wave;                      % Wavenumber
z_s = 2.990:200e-6:3.010;           % Distance between sample and objective
d_slm = 3;                          % Distance between objective and SLM
f_slm = 300;                        % Focal length of diffractive lens
z_h = 500;
%z_h = [50,75,100,125,150];          % Distance between SLM and camera before focus
%z_h = [450,475,500,525,550];       % Distance between SLM and camera after focus
r_h = radius_hologram;
r_h(r_h == 0) = eps;
z_r(z_r == 0) = eps;
N = 1000;                          % No. of photons in hologram
bg = 3;                            % Bg photons            
defocus = (z_s*1e+3-3e+3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = zeros(length(z_h),length(z_s));
fisher_x = zeros(length(z_h),length(z_s));
fisher_z = zeros(length(z_h),length(z_s));
CRLB_x = zeros(length(z_h),length(z_s));
CRLB_z = zeros(length(z_h),length(z_s));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CRLB_x
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
p_noise(1) = subplot(1,2,1);
for j = 1:length(z_h)
    for i = 1:length(z_s)
        A(j,i) = 1/(pi*(r_h(j,i)^2+2*(z_r(j,i)/k)*sin((k/(2*z_r(j,i)))*r_h(j,i)^2)));
        q_sidh = @(a,b) A(j,i).*(1+cos((k/(2*z_r(j,i))).*(a.^2+b.^2)));
        dq_dxs = @(a,b) 2.*A(j,i).*trans_mag(j,i).*(k/(2*z_r(j,i))).*sin((k/(2*z_r(j,i))).*(a.^2+b.^2)).*a;
        F_xx = @(a,b) ((N./(q_sidh(a,b)+bg)).*dq_dxs(a,b).*dq_dxs(a,b));
        polar_Fxx = @(theta,r) F_xx(r.*cos(theta),r.*sin(theta)).*r;
        fisher_x(j,i) = integral2(polar_Fxx,0,2*pi,0,r_h(j,i),'method','iterated');
        CRLB_x(j,i) = 1.e+6/sqrt(fisher_x(j,i));
    end
    plot(defocus,CRLB_x(j,:),'LineWidth',3);
    hold on
end
axis([-10 10 0 50]);
title('CRLB_{xy}');
legend(strcat('z_h =',num2str(z_h'),' mm'),'Location','northwest','FontWeight','bold');
xlabel('Distance between sample and objective (\mum)');
ylabel('\sigma_{x}, \sigma_{y} (nm)');
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CRLB_z
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p_noise(2) = subplot(1,2,2);
for j = 1: length(z_h)
    for i = 1:length(z_s)
        A(j,i) = 1/(pi*(r_h(j,i)^2+2*(z_r(j,i)/k)*sin((k/(2*z_r(j,i)))*r_h(j,i)^2)));
        q_sidh = @(a,b) A(j,i).*(1+cos((k/(2*z_r(j,i))).*(a.^2+b.^2)));
        dq_dzr = @(a,b) ((A(j,i).*(a.^2+b.^2).*k)./(2.*z_r(j,i).^2)).*sin((k/(2*z_r(j,i))).*(a.^2+b.^2));
        F_zz = @(a,b) ((N./(q_sidh(a,b)+bg)).*dq_dzr(a,b).*dq_dzr(a,b).*df_zr(j,i).*df_zr(j,i));
        polar_Fzz = @(theta,r) F_zz(r.*cos(theta),r.*sin(theta)).*r;
        fisher_z(j,i) = integral2(polar_Fzz,0,2*pi,0,r_h(j,i),'method','iterated');
        CRLB_z(j,i) = 1.e+6/sqrt(fisher_z(j,i));
    end
    plot(defocus,CRLB_z(j,:),'LineWidth',3)
    hold on
end
axis([-10 10 0 50]);
title('CRLB_{z}');
legend(strcat('z_h =',num2str(z_h'),' mm'),'Location','northwest','FontWeight','bold');
xlabel('Distance between sample and objective (\mum)');
ylabel('\sigma_{z}(nm)');
hold off






