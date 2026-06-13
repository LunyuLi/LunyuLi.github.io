clear
clc
%
sigma=0.001:0.01:10;
chi=0.001:0.002:3;
s=1/3;
R=1.5;
%
load('stop_oscillation_line_sigma.mat');
load('stop_oscillation_line_chi.mat');
load('Hopf_bifurcation_line_sigma.mat');
load('Hopf_bifurcation_line_chi.mat');
% 平衡态
for i=1:1500
    for j=1:1000
       m=chi(i)-s-1;
       n=4*s*(R-1);
       x=m+1+sqrt(m^2+n);
       T1eqnasy(i,j)=-(x-1)^2/4/s/R;
       T2eqnasy(i,j)=(1-x^2)/4/s/R;
    end
end
SSTdiff=T1eqnasy-T2eqnasy;
% jacobian矩阵特征值
for i=1:1500
    for j=1:1000
       m=chi(i)-s-1;
       n=4*s*(R-1);
       x=m+1+sqrt(m^2+n);
       %盛金公式
       a=1;
       b=-((x-1)*(6*chi(i)-x-3)/(4*s)+R-x-sigma(j)-1);
       c=-((x-1)*(4*sigma(j)*chi(i)+6*chi(i)-sigma(j)*(x+3)-(3*x+1))/(4*s)+(R-x)*(sigma(j)+1)-sigma(j));
       d=-((x-1)*(4*sigma(j)*chi(i)-sigma(j)*(3*x+1))/(4*s)+sigma(j)*(R-x));
       A=b*b-3*a*c;     
       B=b*c-9*a*d;     
       C=c*c-3*b*d;     
       DET = B*B - 4*A*C; 
       if (A == 0) && (B == 0)
            X1(i,j)= -c/b; 
            X2(i,j)=X1(i,j);
            X3(i,j)=X1(i,j);
            X4(i,j)=0;
            X5(i,j)=0;
       end
       if DET >0
          Y1 = A*b + 1.5*a*(-B + sqrt(DET));
          Y2 = A*b + 1.5*a*(-B - sqrt(DET));
          y1 = nthroot(Y1,3);  y2 = nthroot(Y2,3); %这里取Y1和Y2的立方根
          X1(i,j) = (-b-y1-y2)/(3*a);
          X2(i,j) = (-b + 0.5*(y1 + y2))/(3*a);  %X2和X3是第2和第3个特征值的实部
          X3(i,j) = (-b + 0.5*(y1 + y2))/(3*a);
          X4(i,j) = 0.5*sqrt(3)*(y1 - y2)/(3*a); %X4是特征值的虚部，X5是X4的相反数
          X5(i,j) = -0.5*sqrt(3)*(y1 - y2)/(3*a);
       end
       if DET==0 && (A ~= 0) 
           K = (b*c-9*a*d)/(b*b - 3*a*c); 
           K = round(K,14);
           X1(i,j) = -b/a + K;  
           X2(i,j) = -0.5*K; 
           X3(i,j) = X2(i,j);
           X4(i,j) =0;
           X5(i,j)=0;
       end
       if DET < 0
           sqA = sqrt(A);
           T = (A*b - 1.5*a*B)/(A*sqA);
           theta = acos(T);
           csth  = cos(theta/3);
           sn3th = sqrt(3)*sin(theta/3);
           X1(i,j) = (-b - 2*sqA*csth)/(3*a);
           X2(i,j) = (-b + sqA*(csth + sn3th))/(3*a);
           X3(i,j) = (-b + sqA*(csth - sn3th))/(3*a);
           X4(i,j) =-1;%这里是为了画出0等值线
           X5(i,j)=1;
           clear sqA T theta csth
       end
        X1symmetry(i,j)=-1;
        X2symmetry(i,j)=R-1;
        X3symmetry(i,j)=-sigma(j);
    end
end
[sigma1,chi1]=meshgrid(sigma,chi);
% 绘图代码
figure
gca1=subplot(2,3,1);
fill(s*[sigma2,10*ones(1,589),fliplr(sigma5),zeros(1,33)],[chi2,2.265:-0.001:1.677,fliplr(chi5),0.915:0.001:0.947],[1,0.5,0.5],'linestyle', 'none');
title(gca1,'(a) Stability Diagram','fontname','Times New Roman');
axis([0,2.5,0,3]);
xticks([0:0.5:2.5]);
yticks([0:0.5:3]);
set(gca1,'fontsize',18,'fontname','Times New Roman');
set(gca1,'TickDir','both');
xlabel('sσ');
ylabel('χ');

gca2=subplot(2,3,2);
[c2,h2]=contour(s*sigma1,chi1,X2,[-0.7,-0.5,0,0.5,1,10],'k-','ShowText','on','linewidth',1,'Labelspacing',200);
clabel(c2,h2,'fontsize',14,'fontname','Times New Roman');
hold on
[c2,h2]=contour(s*sigma1,chi1,X2,[0,100000000000],'r-','ShowText','on','linewidth',2,'Labelspacing',200);
clabel(c2,h2,'fontsize',14,'fontname','Times New Roman');
xlabel('sσ');
ylabel('χ');
title('(c) Re(eigenvalue 2)','fontname','Times New Roman');
set(gca2,'fontsize',18,'fontname','Times New Roman');
axis([0,2.5,0,3]);
xticks([0:0.5:2.5]);
yticks([0:0.5:3]);
set(gca2,'TickDir','both');

gca3=subplot(2,3,5);
[c3,h3]=contour(s*sigma1,chi1,X3,[-0.7,-0.5,0,0.5,1],'k-','ShowText','on','linewidth',1,'Labelspacing',200);
clabel(c3,h3,'fontsize',14,'fontname','Times New Roman');
hold on
[c3,h3]=contour(s*sigma1,chi1,X3,[0,100000000000],'r-','ShowText','on','linewidth',2,'Labelspacing',200);
clabel(c3,h3,'fontsize',14,'fontname','Times New Roman');
xlabel('sσ');
ylabel('χ');
title('(d) Re(eigenvalue 3)','fontname','Times New Roman');
set(gca3,'fontsize',18,'fontname','Times New Roman');
axis([0,2.5,0,3]);
xticks([0:0.5:2.5]);
yticks([0:0.5:3]);
set(gca3,'TickDir','both');

sigma2=[];
chi2=[];
for j=1:1000
    for i=451:1499
        if (X4(i,j)>0)&&(X4(i+1,j)==-1)
            sigma2=cat(2,sigma2,sigma1(i,j));
            chi2=cat(2,chi2,chi1(i,j));
        end
    end
end

sigma3=[];
chi3=[];
for j=1:1000
    for i=1:450
        if (X4(i,j)==-1)&&(X4(i+1,j)>0)
            sigma3=cat(2,sigma3,sigma1(i,j));
            chi3=cat(2,chi3,chi1(i,j));
        end
    end
end

sigma4=[];
chi4=[];
for j=1:1000
    for i=1:450
        if (X4(i,j)>0)&&(X4(i+1,j)==-1)
            sigma4=cat(2,sigma4,sigma1(i,j));
            chi4=cat(2,chi4,chi1(i,j));
        end
    end
end

sigma5=[];
chi5=[];
for j=1:1000
    for i=1:1499
        if (X2(i,j)<0)&&(X2(i+1,j)>0)
            sigma5=cat(2,sigma5,sigma1(i,j));
            chi5=cat(2,chi5,chi1(i,j));
        end
    end
end

gca4=subplot(2,3,3);
fill(s*[sigma2,10*ones(1,741),10:-0.001:0,zeros(1,2038)],[chi2,2.26:0.001:3,3*ones(1,10001),3:-0.001:0.963],[0.9,0.9,0.9],'linestyle','none');
hold on
fill(s*[sigma3,fliplr(sigma4),zeros(1,515)],[chi3,fliplr(chi4),0.367:0.001:0.881],[0.9,0.9,0.9],'linestyle','none');
hold on
[c4,h4]=contour(s*sigma1,chi1,X4,[0:0.5:2],'k','ShowText','on','linewidth',1,'Labelspacing',1000);
clabel(c4,h4,'fontsize',14,'fontname','Times New Roman');
xlabel('sσ');
ylabel('χ');
title('(e) Im(eigenvalue 2)','fontname','Times New Roman');
set(gca4,'fontsize',18,'fontname','Times New Roman');
axis([0,2.5,0,3]);
xticks([0:0.5:2.5]);
yticks([0:0.5:3]);
set(gca4,'TickDir','both');

gca5=subplot(2,3,6);
fill(s*[sigma2,10*ones(1,741),10:-0.001:0,zeros(1,2038)],[chi2,2.26:0.001:3,3*ones(1,10001),3:-0.001:0.963],[0.9,0.9,0.9],'linestyle','none');
hold on
fill(s*[sigma3,fliplr(sigma4),zeros(1,515)],[chi3,fliplr(chi4),0.367:0.001:0.881],[0.9,0.9,0.9],'linestyle','none');
hold on
[c5,h5]=contour(s*sigma1,chi1,X5,-[0:0.5:2],'k-','ShowText','on','linewidth',1,'Labelspacing',600);
clabel(c5,h5,'fontsize',14,'fontname','Times New Roman');
xlabel('sσ');
ylabel('χ');
title('(f) Im(eigenvalue 3)','fontname','Times New Roman');
set(gca5,'fontsize',18,'fontname','Times New Roman');
axis([0,2.5,0,3]);
xticks([0:0.5:2.5]);
yticks([0:0.5:3]);
set(gca5,'TickDir','both');

gca4=subplot(2,3,4);
SSTdiff_max=max(SSTdiff,[],'all');
levels=linspace(0,SSTdiff_max,11);
[~,h4]=contourf(gca4,s*sigma1,chi1,SSTdiff,levels,'LineStyle','none');
caxis(gca4,[0,SSTdiff_max]);
map=[
    0.00 0.00 0.35
    0.00 0.10 0.65
    0.00 0.28 0.90
    0.20 0.50 1.00
    0.45 0.72 1.00
    0.75 0.90 1.00
    1.00 0.92 0.75
    1.00 0.72 0.45
    0.95 0.45 0.20
    0.75 0.15 0.08
    0.45 0.00 0.00
];
colormap(gca4,map);
cb4=colorbar(gca4);
set(cb4,'Ticks',linspace(0,SSTdiff_max,6),'fontsize',18,'fontname','Times New Roman');
title(gca4,'(b) T1*-T2* Equilibrium','fontname','Times New Roman');
axis(gca4,[0,2.5,0,3]);
yticks(gca4,0:0.5:3);
set(gca4,'fontsize',18,'fontname','Times New Roman');
set(gca4,'TickDir','both');
xlabel(gca4,'s\sigma');
ylabel(gca4,'\chi');
