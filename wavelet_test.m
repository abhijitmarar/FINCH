%wavelet test
clc
%clear all
close all

V0 = In;
%V0 = rand(256,256) + 20;
k1 = [1/16,1/4,3/8,1/4,1/16];
k2 = [1/16,0,1/4,0,3/8,0,1/4,0,1/16];
N1 = 5;
N2 = 9;
%col = zeros(256,256);
%row = zeros(256,256);

%V1 = conv2(k1,k1,V0);
%V2 = conv2(k2,k2,V1);


V1 = conv2(conv2(V0,k1','same'),k1,'same');
V2 = conv2(conv2(V1,k2','same'),k2,'same');
W = V1-V2;