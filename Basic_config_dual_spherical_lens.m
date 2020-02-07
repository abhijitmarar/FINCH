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
wave = 650e-6;                              % Wavelength of light (mm)
delta_c = 16e-3;                            % Pixel size of camera (mm)
z_s = 2.990:200e-6:3.010;                   % Distance between sample and objective
d_slm = 3;                                  % Distance between objective and SLM
f_slm1 = 470;                               % Focal length of shorter focal length mirror
f_slm2 = 523;                               % Focal length of longer focal length mirror
s_fac = (f_slm2-f_slm1)/(f_slm2+f_slm1);    % s-factor
%z_h = (2*f_slm1*f_slm2)/(f_slm1+f_slm2);   % Distance between interferometer and camera (maximum overlap)
z_h = 995;
zh_min = (4*R_o*delta_c)/wave;              % Min distance between interferometer and camera
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%RADIUS OF SHORTER FOCAL LENGTH SPHERICAL WAVE AT CAMERA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms d_SLM f_O z_S z_H f_SLM1 f_SLM2
mat_shorter_wave = simplify([1 z_H;0 1]*[1 0; -1/f_SLM1 1]*[1 d_SLM;0 1]*[1 0; -1/f_O 1]*[1 z_S;0 1]);
rad_shorter_wave = zeros(1,length(z_s));
for i = 1:length(z_s)
    mat_shorter_wave_const = subs(mat_shorter_wave,[d_SLM,f_O,z_S,z_H,f_SLM1],[d_slm,f_o,z_s(i),z_h,f_slm1]);
    rad_shorter_wave(i) = abs(mat_shorter_wave_const(1,2)*NA);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%RADIUS OF LONGER FOCAL LENGTH SPHERICAL WAVE AT CAMERA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mat_longer_wave = simplify([1 z_H;0 1]*[1 0; -1/f_SLM2 1]*[1 d_SLM;0 1]*[1 0; -1/f_O 1]*[1 z_S;0 1]);
rad_longer_wave = zeros(1,length(z_s));
for i = 1:length(z_s)
    mat_longer_wave_const = subs(mat_longer_wave,[d_SLM,f_O,z_S,z_H,f_SLM2],[d_slm,f_o,z_s(i),z_h,f_slm2]);
    rad_longer_wave(i) = abs(mat_longer_wave_const(1,2)*NA);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%RADIUS OF HOLOGRAM AT CAMERA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
radius_hologram = zeros(1,length(z_s));
for i = 1:length(z_s)
    radius_hologram(i) = min(rad_shorter_wave(i),rad_longer_wave(i));
end
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
z_r = zeros(1,length(z_s));
for i = 1:length(z_s)
    if (z_s(i) == f_o)
        z_r(i) = double(subs(zr_focus, [z_H,f_SLM1,f_SLM2],[z_h,f_slm1,f_slm2]));
    else
        z_D_const = subs(z_D,[z_S,f_O,d_SLM],[z_s(i),f_o,d_slm]);
        z_f1_const = subs(z_f1,[z_H,z_D,f_SLM1],[z_h,z_D_const,f_slm1]);
        z_f2_const = subs(z_f2,[z_H,z_D,f_SLM2],[z_h,z_D_const,f_slm2]);
        z_r(i) = double(subs(zr_out_of_focus,[z_f1,z_f2,z_d,f_SLM1,f_SLM2],[z_f1_const,z_f2_const,z_D_const,f_slm1,f_slm2]));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CALCULATION OF TRANSVERSE MAGNIFICATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
f_e = simplify((z_S*f_O)/(f_O-z_S));
trans_mag = zeros(1,length(z_s));
mag_out_of_focus_num = (f_e*z_H);
mag_out_of_focus_den = (f_e+d_SLM);
for i = 1:length(z_s)
    if (z_s(i) == f_o)
        trans_mag(i) = z_h/f_o;
    else
        fe_const = subs(f_e, [f_O,z_S],[f_o,z_s(i)]);
        mag_num_cnst = double(subs(mag_out_of_focus_num,[f_e,z_H],[fe_const,z_h]));
        mag_den_cnst = z_s(i)*double(subs(mag_out_of_focus_den,[f_e,d_SLM],[fe_const,d_slm]));
        trans_mag(i) = mag_num_cnst/mag_den_cnst;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PLOTTING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
defocus = (z_s*1e+3-3e+3);

h(1) = subplot(3,2,1);
plot(defocus,rad_shorter_wave)
title('Radius of spherical wave with shorter focal length')
xlabel('Distance between sample and objective (\mum)')
ylabel('Radius(mm)')


h(2) = subplot(3,2,2);
plot(defocus,rad_longer_wave)
title('Radius of spherical wave with longer focal length')
xlabel('Distance between sample and objective (\mum)')
ylabel('Radius(mm)')

h(3) = subplot(3,2,3);
plot(defocus,radius_hologram)
title('Radius of hologram')
xlabel('Distance between sample and objective (\mum)')
ylabel('Radius(mm)')

h(4) = subplot(3,2,4);
plot(defocus,trans_mag)
title('Transverse magnification')
xlabel('Distance between sample and objective (\mum)')
ylabel('M_{T}')

h(5) = subplot(3,2,5);
plot(defocus,z_r)
title('Reconstruction Distance')
xlabel('Distance between sample and objective (\mum)')
ylabel('z_{r}(mm)')
pos = get(h,'Position');
new = mean(cellfun(@(v)v(1),pos(1:2)));
set(h(5),'Position',[new,pos{end}(2:end)])
