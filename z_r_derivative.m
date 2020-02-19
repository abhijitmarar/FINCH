%Derivative of z_r w.r.t z_s

f_o = 3;
dist = 2.990:200e-6:3.010;
d_slm = 3;                          % Distance between objective and SLM
f_slm = 300;                        % Focal length of diffractive lens
camera_dist = [50,75,100,125,150];
%camera_dist = 150;
defocus = (dist*1e+3-3e+3);

syms fo z fslm dslm z_h z_s

fe = (fo*z_s)/(fo-z_s);
f1 = (fslm*(fe+dslm))/(fslm-fe-dslm);
zr = ((f1+z_h)*(fe+dslm+z_h))/(f1-fe-dslm);

fe_prime = simplify(diff(fe,z_s));
f1_prime = simplify(diff(f1,z_s));
zr_prime = simplify(diff(zr,z_s));

df_zr = zeros(length(z_h),length(z_s));
figure;
for j = 1: length(camera_dist)
    for i = 1:length(dist)
        zr_prime_cnst = simplify(subs(zr_prime,[fo,z_h,dslm,fslm,z_s],[f_o,camera_dist(j),d_slm,f_slm,dist(i)]));
        df_zr(j,i) = double(abs(zr_prime_cnst));
    end
    plot(defocus,df_zr(j,:),'LineWidth',3);
    hold on
end
hold off


