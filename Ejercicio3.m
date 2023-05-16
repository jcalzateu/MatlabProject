clc,close all,clear all;
%% Juan Camilo Alzate Urrego Cod. 316505
%% Taller 2 Informatica II
%Paraboloide eliptico
c=1;
a=5;
b=7;
x=linspace(-100,100,50);
y=linspace(-100,100,50);
[x,y]=meshgrid(x,y);
Z=c*((x.^2)/(a.^2)+(y.^2/(b.^2)));
subplot(1,2,1)
mesh(x,y,Z)
title('Mesh')
xlabel('eje x')
ylabel('eje y')
axis off
grid off
subplot(1,2,2)
surf(x,y,Z)
title('surf')
xlabel('eje x')
ylabel('eje y')
zlabel('eje z')
axis off
grid off