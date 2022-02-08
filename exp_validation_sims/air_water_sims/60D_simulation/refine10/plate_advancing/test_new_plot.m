% This is a matlab code to plot INterface shapes extracted Basilisk
% simulation. 
% Author- Anvesh 
% Date- 7nd jan 2022
% The setWH function is used to generate high quality plots

%% Staic meniscus plot for 22.5mm domain. 

clear all
clc

Domain = 0.0225;  % The domain size is taken from the Basilisk code 
theta0 = 30*pi/180;

x = linspace(0 ,Domain,100);
y = 0.0015570 + 0.0027./(tan(theta0).*exp(x./0.0027));


plot(x , y,'k','LineWidth',2, 'MarkerSize',10 )
axis equal 