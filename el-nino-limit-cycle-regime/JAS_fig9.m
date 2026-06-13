%% Eigenvalues of Jacobian matrix at the equilibrium in s-sigma coordinates
clear
clc

R=1.5;
lam=0.8;
s=0.05:0.001:0.55;
sig=0:0.002:5;

[s1,sig1]=meshgrid(s,sig);
m=lam*R-s1-1;
n=4*s1*(R-1);
x=m+1+sqrt(m.^2+n);

a=1;
b=-((x-1).*(6*lam*R-x-3)./(4*s1)+R-x-sig1-1);
c=-((x-1).*(4*sig1*R*lam+6*lam*R-sig1.*(x+3)-(3*x+1))./(4*s1)+(R-x).*(sig1+1)-sig1);
d=-((x-1).*(4*sig1*R*lam-sig1.*(3*x+1))./(4*s1)+sig1.*(R-x));

A=b.^2-3*a*c;
B=b.*c-9*a*d;
C=c.^2-3*b.*d;
DET=B.^2-4*A.*C;

X1=nan(size(s1));
X2=nan(size(s1));
X3=nan(size(s1));
X4=nan(size(s1));
X5=nan(size(s1));

complexRoot=DET>0;
Y1=A(complexRoot).*b(complexRoot)+1.5*a*(-B(complexRoot)+sqrt(DET(complexRoot)));
Y2=A(complexRoot).*b(complexRoot)+1.5*a*(-B(complexRoot)-sqrt(DET(complexRoot)));
y1=nthroot(Y1,3);
y2=nthroot(Y2,3);
X1(complexRoot)=(-b(complexRoot)-y1-y2)/(3*a);
X2(complexRoot)=(-b(complexRoot)+0.5*(y1+y2))/(3*a);
X3(complexRoot)=X2(complexRoot);
X4(complexRoot)=0.5*sqrt(3)*(y1-y2)/(3*a);
X5(complexRoot)=-X4(complexRoot);

threeReal=DET<0;
sqA=sqrt(A(threeReal));
T=(A(threeReal).*b(threeReal)-1.5*a*B(threeReal))./(A(threeReal).*sqA);
T=max(min(T,1),-1);
theta=acos(T);
csth=cos(theta/3);
sn3th=sqrt(3)*sin(theta/3);
X1(threeReal)=(-b(threeReal)-2*sqA.*csth)/(3*a);
X2(threeReal)=(-b(threeReal)+sqA.*(csth+sn3th))/(3*a);
X3(threeReal)=(-b(threeReal)+sqA.*(csth-sn3th))/(3*a);

multipleRoot=abs(DET)<=1e-10 & A~=0;
K=(b(multipleRoot).*c(multipleRoot)-9*a*d(multipleRoot))./(b(multipleRoot).^2-3*a*c(multipleRoot));
X1(multipleRoot)=-b(multipleRoot)/a+K;
X2(multipleRoot)=-0.5*K;
X3(multipleRoot)=X2(multipleRoot);
X4(multipleRoot)=0;
X5(multipleRoot)=0;

bottomSig=nan(1,numel(s));
for j=1:numel(s)
    idx=find(complexRoot(:,j),1,'first');
    if ~isempty(idx)
        bottomSig(j)=sig(idx);
    end
end
bottomSig=fillmissing(bottomSig,'linear','EndValues','nearest');

limitCycle=complexRoot & X2>0;

fig=figure('Color','w','Units','pixels','Position',[100,100,1020,620]);
t=tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');

axA=nexttile(t,1);
drawStabilityPanel(axA,s1,sig1,limitCycle,bottomSig,s);

axC=nexttile(t,2);
drawRealPanel(axC,s1,sig1,X2,[-0.65,-0.55,-0.4,0.5,1,2],...
    '(c) Re(eigenvalue 2)',false,false);

axE=nexttile(t,3);
drawImagPanel(axE,s1,sig1,X4,bottomSig,s,'(e) Im(eigenvalue 2)',false,false,1);

axB=nexttile(t,4);
drawContourOnlyPanel(axB,s1,sig1,X1,[-4,-3,-2,-1.2,-0.8],...
    '(b) eigenvalue 1 (Real)',true,true);

axD=nexttile(t,5);
drawRealPanel(axD,s1,sig1,X3,[-0.65,-0.55,-0.4,0.5,1,2],...
    '(d) Re(eigenvalue 3)',true,false);

axF=nexttile(t,6);
drawImagPanel(axF,s1,sig1,X5,bottomSig,s,'(f) Im(eigenvalue 3)',true,false,-1);

function drawStabilityPanel(ax,s1,sig1,limitCycle,bottomSig,s)
    contourf(ax,s1,sig1,double(limitCycle),[0.5,1.5],'LineStyle','none');
    hold(ax,'on');
    colormap(ax,[1,0.45,0.48]);
    axis(ax,[0.05,0.55,0,5]);
    setupAxes(ax,'(a) Stability Diagram',false,true);
    text(ax,0.15,1.35,'Limit Cycle',...
        'FontName','Times New Roman','FontWeight','bold',...
        'FontSize',10,'Rotation',55,'HorizontalAlignment','center');
    text(ax,0.18,0.95,'Regime',...
        'FontName','Times New Roman','FontWeight','bold',...
        'FontSize',10,'Rotation',55,'HorizontalAlignment','center');
    text(ax,0.26,2.05,'Hopf Bifurcation',...
        'FontName','Times New Roman','FontWeight','bold',...
        'FontSize',10,'Rotation',-50,'HorizontalAlignment','center');
    quiver(ax,0.29,2.0,-0.08,-0.25,0,...
        'Color','k','LineWidth',0.9,'MaxHeadSize',1.4);
    text(ax,0.075,0.62,'Oscillation',...
        'FontName','Times New Roman','FontSize',8,'Rotation',-65);
    text(ax,0.067,0.46,'Stop',...
        'FontName','Times New Roman','FontSize',8,'Rotation',-65);
end

function drawRealPanel(ax,s1,sig1,Z,levels,titleText,showX,showY)
    hold(ax,'on');
    contour(ax,s1,sig1,Z,levels,'k','LineWidth',0.8,...
        'ShowText','on','LabelSpacing',260);
    Zzero=Z;
    Zzero(sig1<=0.04)=NaN;
    contour(ax,s1,sig1,Zzero,[0,0],'r','LineWidth',1.2,...
        'ShowText','on','LabelSpacing',420);
    setupAxes(ax,titleText,showX,showY);
end

function drawContourOnlyPanel(ax,s1,sig1,Z,levels,titleText,showX,showY)
    contour(ax,s1,sig1,Z,levels,'k','LineWidth',0.8,...
        'ShowText','on','LabelSpacing',260);
    setupAxes(ax,titleText,showX,showY);
end

function drawImagPanel(ax,s1,sig1,Z,bottomSig,s,titleText,showX,showY,signFlag)
    hold(ax,'on');
    patch(ax,[s,fliplr(s)],[zeros(size(s)),fliplr(bottomSig)],[0.88,0.88,0.88],...
        'EdgeColor','none');
    if signFlag>0
        contour(ax,s1,sig1,Z,[0.9,1.2],'k','LineWidth',0.8,...
            'ShowText','on','LabelSpacing',320);
    else
        contour(ax,s1,sig1,Z,[-1.2,-1,-0.9],'k','LineWidth',0.8,...
            'ShowText','on','LabelSpacing',320);
    end
    plot(ax,s,bottomSig,'r--','LineWidth',1.2);
    text(ax,0.43,max(bottomSig(s>=0.43))+0.05,'0',...
        'FontName','Times New Roman','FontSize',10);
    setupAxes(ax,titleText,showX,showY);
end

function setupAxes(ax,titleText,showX,showY)
    axis(ax,[0.05,0.55,0,5]);
    title(ax,titleText,'FontName','Times New Roman','FontSize',13,'FontWeight','normal');
    set(ax,'FontName','Times New Roman','FontSize',11,...
        'Box','on','Layer','top','TickDir','both','LineWidth',0.6);
    xticks(ax,0.05:0.10:0.55);
    yticks(ax,0:1:5);
    if showX
        xlabel(ax,'s','FontName','Times New Roman','FontSize',13,'FontAngle','italic');
        xtickangle(ax,45);
    else
        ax.XTickLabel=[];
    end
    if showY
        ylabel(ax,'\sigma','Interpreter','tex','FontName','Times New Roman','FontSize',13,'FontAngle','italic');
    else
        ax.YTickLabel=[];
    end
end
