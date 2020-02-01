% Some FINCH analysis
% Optics Letters 32(8) 912 (2007)
% Optics Express 19(6) p. 5047 (2011)
% units mm

fo = 3; % objective focal length
d1 = 100; % distance from objective to SLM
zh = 400; % distance from SLM to CCD
fd = 200; % focal length on SLM

dz = 1.e-3*linspace(-5,5,64); % defocus
zs = fo + dz; % distance to object

fe = zs*fo./(fo-zs);
f1 = fd*(fe+d1)./(fd-(fe+d1));

zr = (f1+zh).*(fe+d1+zh)./(f1-fe-d1);

% magnificaton and NA
Mt = zh.*fe./(zs.*(fe+d1));
Dslm = 1;
a = 0; % for a point object
Dh = Dslm*abs(1-(1-a)*zh/fd);

figure();
plot(dz,zr,'-o');
title('reconstruction position');

figure();
plot(dz,Mt);
title('Transverse Magnification');

figure();
out1 = (fo/Dslm)*ones(size(dz));
out2 = abs(zr)./(Mt*Dh);
plot(dz,out1,dz,out2);
title('NA');
legend('input resolution','output resolution');