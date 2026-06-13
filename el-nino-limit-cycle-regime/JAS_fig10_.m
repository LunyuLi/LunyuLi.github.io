%% Fig. 10(b-f): equilibrium SST difference and Jacobian eigenvalues in s-chi coordinates
clear
clc

R=1.5;
sigma=2.5;
s=0.001:0.001:1;
chi=0.001:0.002:3;

[s1,chi1]=meshgrid(s,chi);
xcoord=s1;

m=chi1-s1-1;
n=4*s1*(R-1);
x=m+1+sqrt(m.^2+n);

T1eqnasy=-(x-1).^2./4./s1./R;
T2eqnasy=(1-x.^2)./4./s1./R;
SSTdiff=T1eqnasy-T2eqnasy;

a=1;
b=-((x-1).*(6*chi1-x-3)./(4*s1)+R-x-sigma-1);
c=-((x-1).*(4*sigma*chi1+6*chi1-sigma.*(x+3)-(3*x+1))./(4*s1)+(R-x).*(sigma+1)-sigma);
d=-((x-1).*(4*sigma*chi1-sigma.*(3*x+1))./(4*s1)+sigma.*(R-x));

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

[hopfBoundary,stopBoundary]=extractStabilityBoundaries(xcoord,chi1,X2,complexRoot);

fig=figure('Color','w','Units','pixels','Position',[100,100,1030,620]);
t=tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');

axA=nexttile(t,1);
drawStabilityPanel(axA,sigma,xcoord(1,:),hopfBoundary,stopBoundary);

axC=nexttile(t,2);
drawRealPanel(axC,xcoord,chi1,X2,[-0.7,-0.6,0,0.6,1.2,1.8,5,10,100],...
    '(c) Re(eigenvalue 2)',false,false);

axE=nexttile(t,3);
drawImagPanel(axE,xcoord,chi1,X4,stopBoundary,[0,1,1.3],...
    '(e) Im(eigenvalue 2)',false,false);

axB=nexttile(t,4);
drawSSTPanel(axB,xcoord,chi1,SSTdiff);

axD=nexttile(t,5);
drawRealPanel(axD,xcoord,chi1,X3,[-0.8,-0.7,-0.6,0,0.6,1.2],...
    '(d) Re(eigenvalue 3)',true,false);

axF=nexttile(t,6);
drawImagPanel(axF,xcoord,chi1,X5,stopBoundary,[-1.3,-1,0],...
    '(f) Im(eigenvalue 3)',true,false);

function [hopfBoundary,stopBoundary]=extractStabilityBoundaries(xcoord,chi1,X2,complexRoot)
    xLine=xcoord(1,:);
    hopfBoundary=nan(size(xLine));
    stopBoundary=nan(size(xLine));
    for j=1:numel(xLine)
        idxStop=find(complexRoot(:,j),1,'last');
        if ~isempty(idxStop)
            stopBoundary(j)=chi1(idxStop,j);
        end

        z=X2(:,j);
        chiCol=chi1(:,j);
        valid=isfinite(z);
        crossIdx=find(valid(1:end-1) & valid(2:end) & z(1:end-1).*z(2:end)<=0,1,'last');
        if ~isempty(crossIdx)
            z1=z(crossIdx);
            z2=z(crossIdx+1);
            y1=chiCol(crossIdx);
            y2=chiCol(crossIdx+1);
            if z2~=z1
                hopfBoundary(j)=y1-z1*(y2-y1)/(z2-z1);
            else
                hopfBoundary(j)=y1;
            end
        end
    end
    hopfBoundary=fillmissing(hopfBoundary,'linear','EndValues','nearest');
    stopBoundary=fillmissing(stopBoundary,'linear','EndValues','nearest');
end

function drawStabilityPanel(ax,sigma,xLine,hopfBoundary,stopBoundary)
    hold(ax,'on');
    valid=isfinite(hopfBoundary) & isfinite(stopBoundary) & stopBoundary>hopfBoundary;
    xxReal=xLine(valid);
    hopfReal=hopfBoundary(valid);
    stopReal=stopBoundary(valid);
    patch(ax,[xxReal,fliplr(xxReal)],[hopfReal,fliplr(stopReal)],[1,0.5,0.5],...
        'EdgeColor','none');

    ss=linspace(0,1,400);
    fitX=sigma*ss;
    hopf=0.49*fitX+0.95;
    stop=0.78*fitX+1.17;
    hHopf=plot(ax,ss,hopf,'b-','LineWidth',1.1);
    hStop=plot(ax,ss,stop,'b--','LineWidth',1.1);

    setupAxes(ax,'(a) Stability Diagram',false,true);
    text(ax,0.35,1.55,'Limit Cycle Regime','FontName','Times New Roman',...
        'FontWeight','bold','FontSize',9,'Rotation',28,...
        'HorizontalAlignment','center');
    text(ax,0.25,2.15,'Oscillation Stop','FontName','Times New Roman',...
        'FontWeight','bold','FontSize',9,'Rotation',28);
    quiver(ax,0.28,2.05,0,-0.16,0,'Color','k','LineWidth',0.8,'MaxHeadSize',1.2);
    text(ax,0.30,1.10,'Hopf Bifurcation','FontName','Times New Roman',...
        'FontWeight','bold','FontSize',9,'Rotation',28);
    quiver(ax,0.40,1.25,0.08,0.14,0,'Color','k','LineWidth',0.8,'MaxHeadSize',1.2);
    legend(ax,[hHopf,hStop],{'$\chi=0.49s\sigma+0.95$','$\chi=0.78s\sigma+1.17$'},...
        'Interpreter','latex','Location','southeast','FontSize',7,'Box','on');
end

function drawSSTPanel(ax,xcoord,chi1,SSTdiff)
    Zplot=SSTdiff;
    Zplot(Zplot<0)=0;
    Zplot(Zplot>5)=5;
    levels=linspace(0,5,11);
    contourf(ax,xcoord,chi1,Zplot,levels,'LineStyle','none');
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
    colormap(ax,map);
    caxis(ax,[0,5]);
    cb=colorbar(ax);
    set(cb,'Ticks',0:1:5,'FontSize',11,'FontName','Times New Roman');
    setupAxes(ax,'(b) T1*-T2* Equilibrium',true,true);
end

function drawRealPanel(ax,xcoord,chi1,Z,levels,titleText,showX,showY)
    hold(ax,'on');
    contour(ax,xcoord,chi1,Z,levels,'k','LineWidth',0.8,...
        'ShowText','on','LabelSpacing',260);
    Zzero=Z;
    Zzero(chi1<0.05)=NaN;
    contour(ax,xcoord,chi1,Zzero,[0,0],'r','LineWidth',1.2,...
        'ShowText','on','LabelSpacing',360);
    setupAxes(ax,titleText,showX,showY);
end

function drawImagPanel(ax,xcoord,chi1,Z,stopBoundary,levels,titleText,showX,showY)
    hold(ax,'on');
    xLine=xcoord(1,:);
    stopInterp=interp1(xLine,stopBoundary,xcoord,'linear','extrap');
    valid=isfinite(stopBoundary);
    xx=xLine(valid);
    stop=stopBoundary(valid);
    patch(ax,[xx,fliplr(xx)],[stop,3*ones(size(xx))],[0.88,0.88,0.88],...
        'EdgeColor','none');
    Zmask=Z;
    Zmask(chi1>=stopInterp)=NaN;
    contour(ax,xcoord,chi1,Zmask,levels,'k','LineWidth',0.8,...
        'ShowText','on','LabelSpacing',320);
    plot(ax,xx,stop,'r--','LineWidth',1.2);
    setupAxes(ax,titleText,showX,showY);
end

function setupAxes(ax,titleText,showX,showY)
    axis(ax,[0,1,0,3]);
    title(ax,titleText,'FontName','Times New Roman','FontSize',13,'FontWeight','normal');
    set(ax,'FontName','Times New Roman','FontSize',11,...
        'Box','on','Layer','top','TickDir','both','LineWidth',0.6);
    xticks(ax,0:0.2:1);
    yticks(ax,0:0.5:3);
    if showX
        xlabel(ax,'s','Interpreter','tex',...
            'FontName','Times New Roman','FontSize',13,'FontAngle','italic');
    else
        ax.XTickLabel=[];
    end
    if showY
        ylabel(ax,'\chi','Interpreter','tex',...
            'FontName','Times New Roman','FontSize',13,'FontAngle','italic');
    else
        ax.YTickLabel=[];
    end
end
