clc;
clear all;
close all;

N = 1000;
t = linspace(0,20*pi,100*N);
r = 10;
x = r*cos(t);
y = r*sin(t);
figure
plot3(y,x,t)