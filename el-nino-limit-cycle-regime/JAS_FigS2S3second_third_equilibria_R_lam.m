%% 第二个平衡态
clear
clc
s=1/3;
del=1/2;
R=0.01:0.01:3;
lam=0.01:0.01:3;

for i=1:300
    for j=1:300
       m=lam(i)*R(j)-s-1;
       n=4*s*(R(j)-1);
       if m^2+n<0
           qeqn2(i,j)=NaN;
           continue
       end
       x=m+1+sqrt(m^2+n);
       T1eqn2(i,j)=-(x-1)^2/4/s/R(j);
       T2eqn2(i,j)=(1-x^2)/4/s/R(j);
       qeqn2(i,j)=(x-1)/2/s*R(j)/R(j);
       h1eqn2(i,j)=lam(i)*(x-1)/2/s*R(j)/R(j);
       h2eqn2(i,j)=h1eqn2(i,j)-2*R(j)*lam(i)*(T1eqn2(i,j)-T2eqn2(i,j));
       Tseqn2(i,j)=h2eqn2(i,j)-1;
    end
end
[R1,lam1]=meshgrid(R,lam);
% B=R1.*lam1;
blueOrangeRed=[0.00,0.00,0.45;...
               0.00,0.18,0.70;...
               0.08,0.36,0.86;...
               0.25,0.52,0.96;...
               0.48,0.73,0.95;...
               0.72,0.90,0.96;...
               1.00,0.93,0.70;...
               0.99,0.72,0.42;...
               0.96,0.46,0.18;...
               0.80,0.20,0.02;...
               0.58,0.04,0.00;...
               0.35,0.00,0.00];
colorbarWidth=0.018;
colorbarGap=0.025;
figure
% gca2=subplot(2,2,3);
% [c2,h2]=contour(R1,lam1,T2eqn,[-10:-2,-1,-0.75,-0.5,-0.25,0],'k','ShowText','on','linewidth',1.5,'Labelspacing',400);
% clabel(c2,h2,'fontsize',18);
% axis([1,3,0,2]);
% subtitle(gca2,'(c) T2* eq2');
% set(gca2,'fontsize',16);
% xlabel('R');
% ylabel('χ');
% gca3=subplot(2,2,2);
% [c3,h3]=contour(R1,B,T1eqn,[0,-0.1,-0.25,-0.5,-0.75,-7:-1],'k','ShowText','on','linewidth',1.5,'Labelspacing',400);
% clabel(c3,h3,'fontsize',18);
% axis([1,3,0,2]);
% subtitle(gca3,'(b) T1* eq2');
% set(gca3,'fontsize',16);
% xlabel('R');
% ylabel('χ');
gca5=subplot(1,1,1);
q2clim=[-22.5,22.5];
q2ticks=-20:5:20;
q2levels=linspace(q2clim(1),q2clim(2),size(blueOrangeRed,1)+1);
q2plot=qeqn2;
q2plot(q2plot<q2clim(1))=q2clim(1);
q2plot(q2plot>q2clim(2))=q2clim(2);
[~,h5]=contourf(R1,lam1,q2plot,q2levels,'LineStyle','none');
hold on
axis([0,3,0,3]);
caxis(gca5,q2clim);
colormap(gca5,blueOrangeRed);
cb5=colorbar(gca5);
set(cb5,'Ticks',q2ticks,'TickLabels',arrayfun(@num2str,q2ticks,'UniformOutput',false));
plot(ones(1,301),0:0.01:3,'k-','LineWidth',2);
hold on
plot(0.01:0.01:1,(1+s)./[0.01:0.01:1],'k--','LineWidth',2);
set(cb5,'fontsize',16,'FontName','Times New Roman');
title(gca5,'q_2^*','Interpreter','tex','FontName','Times New Roman');
set(gca5,'fontsize',16,'FontName','Times New Roman');
set(gca5,'Color',[0.88,0.88,0.88]);
xlabel('R','FontName','Times New Roman');
ylabel('Λ','FontName','Times New Roman');
ax5pos=get(gca5,'Position');
set(cb5,'Position',[ax5pos(1)+ax5pos(3)+colorbarGap,ax5pos(2),colorbarWidth,ax5pos(4)]);
% gca1=subplot(2,2,4);
% [c1,h1]=contour(R1,B,h1eqn,[0,0.25,0.5,0.75,1:9],'k','ShowText','on','linewidth',1.5,'Labelspacing',230);
% clabel(c1,h1,'fontsize',18);
% axis([1,3,0,2]);
% subtitle(gca1,'(d) η1* eq2');
% set(gca1,'fontsize',16);
% xlabel('R');
% ylabel('χ');

% figure
% gca2=subplot(2,2,1);
% [c2,h2]=contour(R1,B,h2eqn,[0,-0.25,-0.5,-0.75,-1:-1:-9],'k','ShowText','on','linewidth',1.5,'Labelspacing',300);
% clabel(c2,h2,'fontsize',18);
% axis([1,3,0,2]);
% subtitle(gca2,'(e) η2* eq2');
% set(gca2,'fontsize',16);
% xlabel('R');
% ylabel('χ');
% gca4=subplot(2,2,2);
% [c4,h4]=contour(R1,B,Tseqn,[-1,-1.25,-1.5,-1.75,-2:-1:-10],'k','ShowText','on','linewidth',1.5,'Labelspacing',400);
% clabel(c4,h4,'fontsize',18);
% axis([1,3,0,2]);
% subtitle(gca4,'(f) Ts* eq2');
% set(gca4,'fontsize',16);
% xlabel('R');
% ylabel('χ');

% figure
% gca=subplot(1,1,1);
% plot([0.001:0.001:1],s+1+sqrt(4*s*(1-[0.001:0.001:1])),'k','linewidth',1.5)
% hold on
% plot([0.001:0.001:1],s+1-sqrt(4*s*(1-[0.001:0.001:1])),'k--','linewidth',1.5)
% hold on
% plot(ones(1,1333),0.001:0.001:4/3,'k','linewidth',1.5)
% axis([0,3,0,2]);
% xlabel('R');
% ylabel('χ');
% set(gca,'fontsize',28);
% str = {'q*≥0'};
% text(1.5,1,str,'fontsize',14);
% str = {'q*＜0'};
% text(1.5,1,str,'fontsize',14);
% str = {'Complex Root'};
% text(1.5,1,str,'fontsize',14);
% title('The second equilibrium');

%% 第三个平衡态
for i=1:300
    for j=1:300
       m=lam(i)*R(j)-s-1;
       n=4*s*(R(j)-1);
       if m^2+n<0
           qeqn3(i,j)=NaN;
           continue
       end
       x=m+1-sqrt(m^2+n);
       % T1eqn(i,j)=-(x-1)^2/4/s/R(j);
       % T2eqn(i,j)=(1-x^2)/4/s/R(j);
       qeqn3(i,j)=(x-1)/2/s;
       % h1eqn(i,j)=lam(i)*(x-1)/2/s*R(j)/R(j);
       % h2eqn(i,j)=h1eqn(i,j)-2*R(j)*lam(i)*(T1eqn(i,j)-T2eqn(i,j));
       % Tseqn(i,j)=h2eqn(i,j)-1;
    end
end

figure
gca6=subplot(1,1,1);
q3clim=[-5,5];
q3ticks=-5:1:5;
q3levels=linspace(q3clim(1),q3clim(2),size(blueOrangeRed,1)+1);
q3plot=qeqn3;
q3plot(q3plot<q3clim(1))=q3clim(1);
q3plot(q3plot>q3clim(2))=q3clim(2);
[~,h6]=contourf(R1,lam1,q3plot,q3levels,'LineStyle','none');
hold on;
plot(ones(1,301),0:0.01:3,'k-','LineWidth',2);
hold on
plot(0.01:0.01:1,(1+s)./[0.01:0.01:1],'k--','LineWidth',2);
axis([0,3,0,3]);
caxis(gca6,q3clim);
colormap(gca6,blueOrangeRed);
cb6=colorbar(gca6);
set(cb6,'Ticks',q3ticks,'TickLabels',arrayfun(@num2str,q3ticks,'UniformOutput',false));
set(cb6,'fontsize',16,'FontName','Times New Roman');
title(gca6,'q_3^*','Interpreter','tex','FontName','Times New Roman');
set(gca6,'fontsize',16,'FontName','Times New Roman');
set(gca6,'Color',[0.88,0.88,0.88]);
xlabel('R','FontName','Times New Roman');
ylabel('Λ','FontName','Times New Roman');
ax6pos=get(gca6,'Position');
set(cb6,'Position',[ax6pos(1)+ax6pos(3)+colorbarGap,ax6pos(2),colorbarWidth,ax6pos(4)]);

