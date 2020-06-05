%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SIDH characterization
% 02/1/2020
% Abhijit Marar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SIDH with 2 spherical waves, Simple configuration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DEFINE CONSTANTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f_o = 3;                                    % Focal length of objective (mm)    
NA = 1.42;                                  % Numerical aperture of objective
D_bpp = (2*f_o*NA);                         % Diameter of back pupil plane
R_o = f_o*NA;                               % Radius of beam at interferometer
wave = 670e-6;                              % Wavelength of light (mm)
delta_c = 16e-3;                            % Pixel size of camera (mm)
z_s = 2.990:5e-6:3.010;                    % Distance between sample and objective
d_slm = 3;                                  % Distance between objective and SLM
f_slm1 = 200;                               % Focal length of shorter focal length mirror
f_slm2 = 400;                               % Focal length of longer focal length mirror
s_fac = (f_slm2-f_slm1)/(f_slm2+f_slm1);    % s-factor
%z_h = (2*f_slm1*f_slm2)/(f_slm1+f_slm2);   % Distance between interferometer and camera (maximum overlap)
%z_h = [500,525,550,575,600];
z_h = 150;
zh_min = (4*R_o*delta_c)/wave;              % Min distance between interferometer and camera
defocus = (z_s*1e+3-3e+3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%RADIUS OF SHORTER FOCAL LENGTH SPHERICAL WAVE AT CAMERA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms d_SLM f_O z_S z_H f_SLM1 f_SLM2
mat_shorter_wave = simplify([1 z_H;0 1]*[1 0; -1/f_SLM1 1]*[1 d_SLM;0 1]*[1 0; -1/f_O 1]*[1 z_S;0 1]);
rad_shorter_wave = zeros(length(z_h),length(z_s));
figure;
for j = 1:length(z_h)
    for i = 1:length(z_s)
        mat_shorter_wave_const = subs(mat_shorter_wave,[d_SLM,f_O,z_S,z_H,f_SLM1],[d_slm,f_o,z_s(i),z_h(j),f_slm1]);
        rad_shorter_wave(j,i) = abs(mat_shorter_wave_const(1,2)*NA);
    end
    plot(defocus,rad_shorter_wave(j,:),'LineWidth',6);
    hold on    
end
title('Radius of f_{d1}');
xlabel('Distance between sample and objective (\mum)');
ylabel('Radius(mm)');
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%RADIUS OF LONGER FOCAL LENGTH SPHERICAL WAVE AT CAMERA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mat_longer_wave = simplify([1 z_H;0 1]*[1 0; -1/f_SLM2 1]*[1 d_SLM;0 1]*[1 0; -1/f_O 1]*[1 z_S;0 1]);
rad_longer_wave = zeros(length(z_h),length(z_s));
figure;
for j = 1:length(z_h)
    for i = 1:length(z_s)
        mat_longer_wave_const = subs(mat_longer_wave,[d_SLM,f_O,z_S,z_H,f_SLM2],[d_slm,f_o,z_s(i),z_h(j),f_slm2]);
        rad_longer_wave(j,i) = abs(mat_longer_wave_const(1,2)*NA);
    end
    plot(defocus,rad_longer_wave(j,:),'LineWidth',6);
    hold on   
end
title('Radius of f_{d2}');
xlabel('Distance between sample and objective (\mum)');
ylabel('Radius(mm)');
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%RADIUS OF HOLOGRAM AT CAMERA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
radius_hologram = zeros(length(z_h),length(z_s));
figure;
for j = 1:length(z_h)
    for i = 1:length(z_s)
        radius_hologram(j,i) = min(rad_shorter_wave(j,i),rad_longer_wave(j,i));
    end
    plot(defocus,radius_hologram(j,:),'LineWidth',6);
    hold on
end
title('Radius of Hologram');
xlabel('Distance between sample and objective (\mum)');
ylabel('Radius(mm)');
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CALCULATION OF RECONSTRUCTION DISTANCE (Z_R)
%Siegel, Nisan, Joseph Rosen, and Gary Brooker. Optics express 20.18 (2012): 19822-19835.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms z_d
zr_focus = simplify(abs(((z_H-f_SLM1)*(z_H-f_SLM2))/(f_SLM1-f_SLM2)));
z_D = simplify((z_S*(f_O-d_SLM)+f_O*d_SLM)/(f_O-z_S));
z_f1 = simplify((z_H*z_D)-f_SLM1*(z_D+z_H));
z_f2 = simplify((z_H*z_D)-f_SLM2*(z_D+z_H));
zr_out_of_focus = simplify(abs((z_f1*z_f2)/(z_d^2*(f_SLM1-f_SLM2))));
z_r = zeros(length(z_h),length(z_s));
figure;
for j = 1:length(z_h)
    for i = 1:length(z_s)
        if (z_s(i) == f_o)
            z_r(j,i) = double(subs(zr_focus, [z_H,f_SLM1,f_SLM2],[z_h(j),f_slm1,f_slm2]));
        else
            z_D_const = subs(z_D,[z_S,f_O,d_SLM],[z_s(i),f_o,d_slm]);
            z_f1_const = subs(z_f1,[z_H,z_D,f_SLM1],[z_h(j),z_D_const,f_slm1]);
            z_f2_const = subs(z_f2,[z_H,z_D,f_SLM2],[z_h(j),z_D_const,f_slm2]);
            z_r(j,i) = double(subs(zr_out_of_focus,[z_f1,z_f2,z_d,f_SLM1,f_SLM2],[z_f1_const,z_f2_const,z_D_const,f_slm1,f_slm2]));
        end    
    end
    plot(defocus,z_r(j,:),'LineWidth',6);
    hold on
end
title('Reconstruction Distance');
xlabel('Distance between sample and objective (\mum)');
ylabel('Radius(mm)');
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CALCULATION OF TRANSVERSE MAGNIFICATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
f_e = simplify((z_S*f_O)/(f_O-z_S));
trans_mag = zeros(length(z_h),length(z_s));
mag_out_of_focus_num = (f_e*z_H);
mag_out_of_focus_den = (f_e+d_SLM);
figure;
for j = 1:length(z_h)
    for i = 1:length(z_s)
        if (z_s(i) == f_o)
            trans_mag(j,i) = z_h(j)/f_o;
        else
            fe_const = subs(f_e, [f_O,z_S],[f_o,z_s(i)]);
            mag_num_cnst = double(subs(mag_out_of_focus_num,[f_e,z_H],[fe_const,z_h(j)]));
            mag_den_cnst = z_s(i)*double(subs(mag_out_of_focus_den,[f_e,d_SLM],[fe_const,d_slm]));
            trans_mag(j,i) = mag_num_cnst/mag_den_cnst;
        end
    end
    plot(defocus,trans_mag(j,:),'LineWidth',6);
    hold on
end
title('Transverse Magnification');
legend(strcat('z_h =',num2str(z_h'),' mm'),'Location','best','FontWeight','bold');
xlabel('Distance between sample and objective (\mum)');
ylabel('Radius(mm)');
hold off
