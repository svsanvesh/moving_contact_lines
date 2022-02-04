% This is a matlab code to plot INterface shapes extracted Basilisk
% simulation. 
% Author- Anvesh 
% Date- 2nd jan 2022
% The setWH function is used to generate high quality plots
%%
clear all
clc

Domain = 0.015;  % The domain size is taken from the Basilisk code i.e. 15mm 
theta0 = 30*pi/180;



% The shape of the static initial meniscus is plotted for reference . 

x = linspace(-Domain/2,Domain/2,100);
yinit = -0.0013 + 0.0027./(tan(theta0).*exp((x+ 0.0075)/0.0027));

figure(1)
plot(x +Domain/2,yinit,'k','LineWidth',2, 'MarkerSize',10 )
ylabel('Downward moving plate','FontWeight','Bold','FontSize',20);
xlabel('15mm domain ','FontWeight','Bold','FontSize',20);
set(gca,'Fontsize',12,'Fontname','Times New Roman')
set(gca,'XMinorTick','on','YMinorTick','on','TickDir','in')
axis equal 
hold on 
setWH(gcf,9,6)
%% Time step - 4.6
% The file ' interface_profile_t4.600000.dat ' is loaded in as an array and is
% used to plot the interafce shape at the corresponding timestep. 

A=load('interface_profile_t0.000000.dat');
A= sortrows(A);

plot(A(:,1) +Domain/2,A(:,2), 'b','LineWidth',2, 'MarkerSize',10 )

ylabel('Downward moving plate','FontWeight','Bold','FontSize',20);
xlabel('Particle diameter (\mum)','FontWeight','Bold','FontSize',20);
set(gca,'Fontsize',12,'Fontname','Times New Roman')
set(gca,'XMinorTick','on','YMinorTick','on','TickDir','in')
axis equal 
setWH(gcf,9,6)
hold on
%% Time step - 4.7
% The file ' interface_profile_t4.600000.dat ' is loaded in as an array and is
% used to plot the interafce shape at the corresponding timestep. 

B=load('interface_profile_t4.700000.dat');
B= sortrows(B);

plot(B(:,1) +Domain/2,B(:,2), 'm','LineWidth',2, 'MarkerSize',10 )

ylabel('Downward moving plate','FontWeight','Bold','FontSize',20);
xlabel('Particle diameter (\mum)','FontWeight','Bold','FontSize',20);
set(gca,'Fontsize',12,'Fontname','Times New Roman')
set(gca,'XMinorTick','on','YMinorTick','on','TickDir','in')
axis equal 
setWH(gcf,9,6)
hold on




