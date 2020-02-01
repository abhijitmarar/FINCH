% FINCH reconstruction

theta1 = 0;
z1 = 0;
z2 = -6;
y2 = .3;
I1 = finchH(0,0,z1,theta1) + finchH(0,y2,z2,theta1); % + finchH(1,0.0,-72,theta1);
theta2 = 2*pi/3;
I2 = finchH(0,0,z1,theta2) + finchH(0,y2,z2,theta2); % + finchH(1,0.0,-72,theta2);
theta3 = 4*pi/3;
I3 = finchH(0,0,z1,theta3) + finchH(0,y2,z2,theta3); % + finchH(1,0.0,-72,theta3);

% add noise
Nph = 2500;
I1 = (Nph/sum(sum(I1)))*I1;
I2 = (Nph/sum(sum(I2)))*I2;
I3 = (Nph/sum(sum(I3)))*I3;
I1 = poissrnd(I1);
I2 = poissrnd(I1);
I3 = poissrnd(I3);

% assemble

IF = I1*(exp(1j*theta3)-exp(1j*theta2)) + ...
    I2*(exp(1j*theta1)-exp(1j*theta3)) + ...
    I3*(exp(1j*theta2)-exp(1j*theta1));

figure();
imagesc(abs(IF));
figure();
imagesc(angle(IF));

wave = 0.0005;
%z = -111.3419;
%z= -170;
z = -153.1367;

N = 2*128;
x = linspace(-1,1,N);
[xx, yy] = meshgrid(x);
s = fftshift(ifft2(fft2(fftshift(IF)).*fft2(fftshift(exp((1j*pi/wave/z)*(xx.^2 + yy.^2))))));

figure();
imagesc(abs(s));
colorbar;