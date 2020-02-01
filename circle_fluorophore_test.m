% FINCH reconstruction

theta1 = 0;
z1 = 3.000;
z2 = 2.995;
y2 = .005;
%radius = .005;
%phi = linspace(0,2*pi,20);
%x_coord = radius*cos(phi);
%y_coord = radius*sin(phi);

I1 = finchH_objective(.005*cosd(270),.005*sind(270),z1,theta1) + finchH_objective(.005*cosd(60),.005*sind(60),z2,theta1) +...
    + finchH_objective(.005*cosd(120),.005*sind(120),z1,theta1);
theta2 = 2*pi/3;
I2 = finchH_objective(.005*cosd(270),.005*sind(270),z1,theta2) + finchH_objective(.005*cosd(60),.005*sind(60),z2,theta2) +...
    + finchH_objective(.005*cosd(120),.005*sind(120),z1,theta2);
theta3 = 4*pi/3;
I3 = finchH_objective(.005*cosd(270),.005*sind(270),z1,theta3) + finchH_objective(.005*cosd(60),.005*sind(60),z2,theta3) +...
    + finchH_objective(.005*cosd(120),.005*sind(120),z1,theta3);


figure();
imagesc(I3);

%{
figure();
imagesc(I3_moved);
I_complete = I3+I3_moved;
figure();
imagesc(I_complete);
%}

%{
% add noise
Nph = 5000;
I1 = (Nph/sum(sum(I1)))*I1;
I2 = (Nph/sum(sum(I2)))*I2;
I3 = (Nph/sum(sum(I3)))*I3;
I1 = poissrnd(I1);
I2 = poissrnd(I2);
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
%z= -306.0482;
z = -260.8631;
%z = 281.2500;

N = 2*256;
x = linspace(-1,1,N);
[xx, yy] = meshgrid(x);
%xx = xx*(400/3); yy = yy*(400/3);
s = fftshift(ifft2(fft2(fftshift(IF)).*fft2(fftshift(exp((1j*pi/wave/z)*(xx.^2 + yy.^2))))));

figure();
imagesc(abs(s));
colorbar;
%}