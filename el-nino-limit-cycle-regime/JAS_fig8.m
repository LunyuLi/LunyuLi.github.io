%% Stability diagram and eigenvalues for the zonally asymmetric equilibrium
clear
clc

s=1/3;
del=1/2;
R=1:0.002:5;
chi=0:0.002:2;

[R1,chi1]=meshgrid(R,chi);
m=chi1-s-1;
n=4*s*(R1-1);
x=m+1+sqrt(m.^2+n);

a=1;
b=-((x-1).*(6*chi1-x-3)./(4*s)+R1-x-del-1);
c=-((x-1).*(4*del*chi1+6*chi1-del.*(x+3)-(3*x+1))./(4*s)+(R1-x).*(del+1)-del);
d=-((x-1).*(4*del*chi1-del.*(3*x+1))./(4*s)+del.*(R1-x));

A=b.^2-3*a*c;
B=b.*c-9*a*d;
C=c.^2-3*b.*d;
DET=B.^2-4*A.*C;

X1=nan(size(R1));
X2=nan(size(R1));
X3=nan(size(R1));
X4=nan(size(R1));
X5=nan(size(R1));

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

topChi=nan(1,numel(R));
for j=1:numel(R)
    idx=find(complexRoot(:,j),1,'last');
    if ~isempty(idx)
        topChi(j)=chi(idx);
    end
end
topChi=fillmissing(topChi,'linear','EndValues','nearest');

limitCycle=complexRoot & X2>0;

fig=figure('Color','w','Units','pixels','Position',[100,100,1020,620]);
t=tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');

axA=nexttile(t,1);
drawStabilityPanel(axA,R1,chi1,limitCycle);

axC=nexttile(t,2);
drawRealPanel(axC,R1,chi1,X2,[-1.5,-1,-0.5,-0.25,0.5,1,1.5,5],...
    '(c) Re(eigenvalue 2)',false,false);

axE=nexttile(t,3);
drawImagPanel(axE,R1,chi1,X4,topChi,R,'(e) Im(eigenvalue 2)',false,false,1);

axB=nexttile(t,4);
drawContourOnlyPanel(axB,R1,chi1,X1,[-0.9,-0.8,-0.7,-0.6],...
    '(b) eigenvalue 1 (Real)',true,true);

axD=nexttile(t,5);
drawRealPanel(axD,R1,chi1,X3,[-1.5,-1,-0.5,-0.25,0.5,1,1.5],...
    '(d) Re(eigenvalue 3)',true,false);

axF=nexttile(t,6);
drawImagPanel(axF,R1,chi1,X5,topChi,R,'(f) Im(eigenvalue 3)',true,false,-1);

function drawStabilityPanel(ax,R1,chi1,limitCycle)
    contourf(ax,R1,chi1,double(limitCycle),[0.5,1.5],...
        'LineStyle','none');
    hold(ax,'on');
    colormap(ax,[1,0.45,0.48]);
    axis(ax,[1,5,0,2]);
    setupAxes(ax,'(a) Stability Diagram',false,true);
    text(ax,2.8,1.66,'Oscillation Stop',...
        'FontName','Times New Roman','FontWeight','bold',...
        'HorizontalAlignment','center','FontSize',10);
    quiver(ax,2.75,1.58,0,-0.10,0,...
        'Color','k','LineWidth',1,'MaxHeadSize',1.5);
    text(ax,2.7,1.13,'Limit Cycle Regime',...
        'FontName','Times New Roman','FontWeight','bold',...
        'HorizontalAlignment','center','FontSize',10);
    text(ax,2.55,0.55,'Hopf Bifurcation',...
        'FontName','Times New Roman','FontWeight','bold',...
        'HorizontalAlignment','center','FontSize',10);
    quiver(ax,2.95,0.66,0.18,0.20,0,...
        'Color','k','LineWidth',1,'MaxHeadSize',1.5);
end

function drawRealPanel(ax,R1,chi1,Z,levels,titleText,showX,showY)
    hold(ax,'on');
    contour(ax,R1,chi1,Z,levels,'k','LineWidth',0.8,...
        'ShowText','on','LabelSpacing',280);
    contour(ax,R1,chi1,Z,[0,0],'r','LineWidth',1.2,...
        'ShowText','on','LabelSpacing',380);
    setupAxes(ax,titleText,showX,showY);
end

function drawContourOnlyPanel(ax,R1,chi1,Z,levels,titleText,showX,showY)
    contour(ax,R1,chi1,Z,levels,'k','LineWidth',0.8,...
        'ShowText','on','LabelSpacing',240);
    setupAxes(ax,titleText,showX,showY);
end

function drawImagPanel(ax,R1,chi1,Z,topChi,R,titleText,showX,showY,signFlag)
    hold(ax,'on');
    patch(ax,[R,fliplr(R)],[topChi,2*ones(size(topChi))],[0.88,0.88,0.88],...
        'EdgeColor','none');
    if signFlag>0
        contour(ax,R1,chi1,Z,[0.5,1,1.5,2],'k','LineWidth',0.8,...
            'ShowText','on','LabelSpacing',320);
    else
        contour(ax,R1,chi1,Z,[-2,-1.5,-1,-0.5],'k','LineWidth',0.8,...
            'ShowText','on','LabelSpacing',320);
    end
    plot(ax,R,topChi,'r--','LineWidth',1.2);
    text(ax,1.75,topChi(find(R>=1.75,1))+0.02,'0',...
        'FontName','Times New Roman','FontSize',10);
    setupAxes(ax,titleText,showX,showY);
end

function setupAxes(ax,titleText,showX,showY)
    axis(ax,[1,5,0,2]);
    title(ax,titleText,'FontName','Times New Roman','FontSize',13,'FontWeight','normal');
    set(ax,'FontName','Times New Roman','FontSize',11,...
        'Box','on','Layer','top','TickDir','both','LineWidth',0.6);
    xticks(ax,1:1:5);
    yticks(ax,0:0.5:2);
    if showX
        xlabel(ax,'R','FontName','Times New Roman','FontSize',14,'FontAngle','italic');
    else
        ax.XTickLabel=[];
    end
    if showY
        ylabel(ax,'\chi','Interpreter','tex','FontName','Times New Roman','FontSize',13,'FontAngle','italic');
    else
        ax.YTickLabel=[];
    end
end
