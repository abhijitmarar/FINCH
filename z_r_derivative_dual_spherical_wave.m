%Derivative of z_r w.r.t z_s for SIDH with dual spherical waves

f_o = 3;
dist = 2.990:200e-6:3.010;

d_slm = 3;                                  % Distance between objective and SLM
f_slm1 = 100;                               % Focal length of shorter focal length mirror
f_slm2 = 300;                               % Focal length of longer focal length mirror

camera_dist = 200;

syms fo z fslm1 fslm2 dslm z_h z_s

z_d = simplify((z_s*(fo-dslm)+fo*dslm)/(fo-z_s));
z_f1 = simplify((z_h*z_d)-fslm1*(z_d+z_h));
z_f2 = simplify((z_h*z_d)-fslm2*(z_d+z_h));

zr_out_of_focus = simplify((z_f1*z_f2)/(z_d^2*(fslm1-fslm2)));

zd_prime = simplify(diff(z_d,z_s));
zf1_prime = simplify(diff(z_f1,z_s));
zf2_prime = simplify(diff(z_f2,z_s));
zr_prime = simplify(diff(zr_out_of_focus,z_s));

df_zr = zeros(1,length(dist));

for i = 1:length(dist)
    if (dist(i) == f_o)
        df_zr(i) = Inf;
    else
        zr_prime_cnst = simplify(subs(zr_prime,[fo,z_h,dslm,fslm1,fslm2,z_s],[f_o,camera_dist,d_slm,f_slm1,f_slm2,dist(i)]));
        df_zr(i) = double(abs(zr_prime_cnst));
    end
end
   
defocus = (dist*1e+3-3e+3);

figure()
plot(defocus,df_zr);