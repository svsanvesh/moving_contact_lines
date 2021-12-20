clear all;
clc



x = linspace(-Domain/2,Domain/2,100);
yinit = -0.0013 + 0.0027./(tan(theta0).*exp((x+ 0.0075)/0.0027));

ms=10;
lw=1
plot(x,y,'-k*','LineWidth',lw, 'MarkerSize',ms)
ylabel('amp','FontWeight','Bold','FontSize',20);
xlabel('Particle diameter (\mum)','FontWeight','Bold','FontSize',20);
set(gca,'Fontsize',12,'Fontname','Times New Roman')
set(gca,'XMinorTick','on','YMinorTick','on','TickDir','in')
% set(gca,'Color','white');                
 xlim([0.5,10]);
%ylim([5e-4,1e3]);
setWH(gcf,4,3)
%%
print(gcf,'-dpng','-r300','test.png')
% print(gcf,'-deps','-r300','test')
