% This is a matlab code to plot INterface shapes extracted Basilisk
% simulation. 
% Author- Anvesh 
% Date- 7nd jan 2022
% The setWH function is used to generate high quality plots

% %% Staic meniscus plot for 22.5mm domain. 
% 
% clear all
% clc
% 
% Domain = 0.0225;  % The domain size is taken from the Basilisk code 
% theta0 = 30*pi/180;
% 
% x = linspace(0 ,Domain,100);
% y = 0.0015570 + 0.0027./(tan(theta0).*exp(x./0.0027));
% 
% 
% plot(x , y,'k','LineWidth',2, 'MarkerSize',10 )
% axis equal 
% hold on
%% Time step - 4.6
% The file ' interface_profile_t4.600000.dat ' is loaded in as an array and is
% used to plot the interafce shape at the corresponding timestep. 

S= importdata('interface_profile_t1.200000.dat');

S= sortrows(S);

plot(S(:,1) ,S(:,2), 'k','LineWidth',1.5, 'MarkerSize',10 )

ylabel('Downward moving plate','FontWeight','Bold','FontSize',20);
xlabel('Particle diameter (\mum)','FontWeight','Bold','FontSize',20);
set(gca,'Fontsize',12,'Fontname','Times New Roman')
set(gca,'XMinorTick','on','YMinorTick','on','TickDir','in')
axis equal 
setWH(gcf,9,6)
hold on




%% Time step - 4.6



% The file ' interface_profile_t4.600000.dat ' is loaded in as an array and is
% used to plot the interafce shape at the corresponding timestep. 

A= importdata('interface_profile_t15.200000.dat');

A= sortrows(A);

plot(A(:,1) ,A(:,2), 'm','LineWidth',1.5, 'MarkerSize',10 )

ylabel('Downward moving plate','FontWeight','Bold','FontSize',20);
xlabel('Particle diameter (\mum)','FontWeight','Bold','FontSize',20);
set(gca,'Fontsize',12,'Fontname','Times New Roman')
set(gca,'XMinorTick','on','YMinorTick','on','TickDir','in')
axis equal 
setWH(gcf,9,6)
hold on
%% Time step - 4.6
% The file ' interface_profile_t4.600000.dat ' is loaded in as an array and is
% used to plot the interafce shape at the corresponding timestep. 

B= importdata('interface_profile_t40.100000.dat');

B= sortrows(B);

plot(B(:,1) ,B(:,2), 'b','LineWidth',1.5, 'MarkerSize',10 )

ylabel('Downward moving plate','FontWeight','Bold','FontSize',20);
xlabel('Particle diameter (\mum)','FontWeight','Bold','FontSize',20);
set(gca,'Fontsize',12,'Fontname','Times New Roman')
set(gca,'XMinorTick','on','YMinorTick','on','TickDir','in')
axis equal 
setWH(gcf,9,6)
hold on
