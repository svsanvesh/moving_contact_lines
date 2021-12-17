clear all
clc


Domain = 0.015;
theta0 = 30*pi/180;

x = linspace(-Domain/2,Domain/2,20);
yinit = -0.0013 + 0.0027./(tan(theta0).*exp((x+ 0.0075)/0.0027));

X = linspace(0,Domain,100);
Y = -0.0013 + 0.0027./(tan(theta0).*exp(X/0.0027));

figure
plot(x+Domain/2,yinit,'ko');
hold on

A=load('interface_profile_t0.000000.dat');
%A=sort(A,1);

plot(A(:,1)+Domain/2,A(:,2),'r.');


B=load('interface_profile_t5.000000.dat');
%A=sort(A,1);

plot(B(:,1)+Domain/2,B(:,2),'b.'); 
hold off

figure(1)
hold on
plot(X,Y,'k--')
hold off

figure(1)
axis equal



figure(1)
hold on
B1 = sort(B,1,'ascend');
plot(B1(:,1)+Domain/2,B1(:,2),'k')

%% Test the fit. 
p1 = Interface_Shape_Steady.p1;
p2 = Interface_Shape_Steady.p2;
p3 = Interface_Shape_Steady.p3;
p4 = Interface_Shape_Steady.p4;
p5 = Interface_Shape_Steady.p5;
p6 = Interface_Shape_Steady.p6;
p7 = Interface_Shape_Steady.p7;
p8 = Interface_Shape_Steady.p8;
p9 = Interface_Shape_Steady.p9;

xfit = B(:,1);
yfit = p1*xfit.^8 + p2*xfit.^7 + p3*xfit.^6 + p4*xfit.^5 + p5*xfit.^4 + p6*xfit.^3 + p7*xfit.^2 + p8*xfit + p9;
figure(1)
hold on
plot(xfit+Domain/2,yfit,'k*')
hold off


