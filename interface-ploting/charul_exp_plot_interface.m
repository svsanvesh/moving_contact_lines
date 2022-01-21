% This is a matlab code to plot INterface shapes extracted Basilisk
% simulation. 
% Author- Anvesh 
% Date- 21st jan 2022
% The setWH function is used to generate high quality plots

clear all
clc



A=load('interface_data.mat');

plot(A.curve(1,:) ,A.curve(2,:), 'b','LineWidth',2, 'MarkerSize',10 )
axis equal