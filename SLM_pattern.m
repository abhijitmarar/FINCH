clear all
N = 2*256;
x = linspace(-3,3,N);
[xx, yy] = meshgrid(x);

orig = (1:1:512)';

wave = 0.515;
a = 400;
theta1 = 0;
theta2 = 2*pi/3;
theta3 = 4*pi/3;

randx_spher = randperm(N,256)';
randy_spher = randperm(N,256)';
randx_plane = setdiff(orig,randx_spher);
randy_plane = setdiff(orig,randy_spher);

xx_spher = xx(randx_spher,randy_spher);
yy_spher = yy(randx_spher,randy_spher);
R_phase1 = ((pi/wave/a)*(xx_spher.^2+yy_spher.^2)-theta2);

xx_plane = xx(randx_plane,randy_plane);
yy_plane = yy(randx_plane,randy_plane);
R_plane = (xx_spher.^2+yy_spher.^2)*0;
%R_phase2 = ((pi/wave/a)*(xx.^2+yy.^2)-theta2);
%R_phase3 = ((pi/wave/a)*(xx.^2+yy.^2)-theta3);
R = 0.5+0.5*exp((-1j*pi/wave/a)*(xx.^2+yy.^2)-1j*theta3);

imagesc(R)
figure()
imagesc(R_plane)
%figure()
%imagesc(R_phase3)



