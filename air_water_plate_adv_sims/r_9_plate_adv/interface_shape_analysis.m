clear all
clc


Domain = 0.015;
theta0 = 30*pi/180;

x = linspace(-Domain/2,Domain/2,20);
yinit = -0.0013 + 0.0027./(tan(theta0).*exp((x+ 0.0075)/0.0027));

X = linspace(0,Domain,100);
Y = -0.0013 + 0.0027./(tan(theta0).*exp(X/0.0027));

figure
plot(x,yinit,'ko');
hold on

A=load('interface_profile_t0.000000.dat');
%A=sort(A,1);

plot(A(:,1),A(:,2),'r.');


B=load('interface_profile_t5.000000.dat');
%A=sort(A,1);

plot(B(:,1),B(:,2),'b*'); 
hold off

% B(:,1) = B(:,1) + Domain/2 ;
% 
 final_x_dat= B(:,1);
 final_y_dat= B(:,2);

plot(final_x_dat, final_y_dat,'r *');

figure(1)
hold on
plot(X,Y,'k--')
hold off

figure(1)
% axis equal




%% Test the fit. 
p1 = steady_state_30d_init_meniscus.p1;
p2 = steady_state_30d_init_meniscus.p2;
p3 = steady_state_30d_init_meniscus.p3;
p4 = steady_state_30d_init_meniscus.p4;
p5 = steady_state_30d_init_meniscus.p5;
p6 = steady_state_30d_init_meniscus.p6;
p7 = steady_state_30d_init_meniscus.p7;
p8 = steady_state_30d_init_meniscus.p8;
p9 = steady_state_30d_init_meniscus.p9;

xfit = B(:,1);
yfit = p1*xfit.^8 + p2*xfit.^7 + p3*xfit.^6 + p4*xfit.^5 + p5*xfit.^4 + p6*xfit.^3 + p7*xfit.^2 + p8*xfit + p9;
figure(1)
hold on
plot(xfit,yfit,'r*')
hold off


