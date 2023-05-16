clc,close all,clear all;
%% Juan Camilo Alzate Urrego Cod. 316505
%% Taller 2 Informatica II
% Ecuación movimiento uniforme
    t=0:1:3600; %Segundos
    t1=0:1:10;  %Segundos
    Xo=0;
    V=5; %Velociad 5 m/s
    x=Xo+V*t;
    V1=(x-Xo)./t;
    V1(1)=0;
    a1=zeros(size(t));
% Ecuación movimiento uniformemente acelerado
    a=2.5; %m/s2
    Vo=0;  %Velocidad inicial
    x1=Xo+Vo*t1+(1/2)*a*t1.^2;
    V=Vo+a*t1;
    a=(V-Vo)./t1;
%Ecuación movimiento armoico simple
    T=1000;    %Periodo
    f=1/T; %frecuencia
    w=2*pi()*f; %Frecuencia Angular
    A=2; %amplitud
    y=A*sin(w*t);
% Graficas
    figure('name','Ejercicio 1','NumberTitle','off','Color','w')
%MU
    subplot(1,3,1)
    plot(t,x,'--sb','LineWidth',0.5,'Markersize',0.5,'MarkerEdgeColor','b','MarkerFacecolor','r')
    title('Movimiento Uniforme','Fontname','TimesNewRoman','Fontsize',12)
    ylabel('Unidades','Fontname','TimesNewRoman','Fontsize',15,'color','red','FontWeight','bold')
    xlabel('Tiempo','Fontname','Times New Roman','Fontsize',10,'FontWeight','bold','rotation',15),hold on
    plot(t,V1,'--pg','LineWidth',0.5,'Markersize',0.5,'MarkerEdgeColor','g','MarkerFacecolor','r')
    plot(t,a1,':^y','LineWidth',0.5,'Markersize',0.5,'MarkerEdgeColor','r','MarkerFaceColor','b')
    legend('posición','Velocidad','Aceleración','location','best')
    grid minor
%MUA
    subplot(1,3,2)
    plot(t1,x1,'--*r','LineWidth',0.5,'Markersize',0.5,'MarkerEdgeColor','g','MarkerFaceColor','g')
    title('Movimiento uniformemente acelerado','Fontname','Times New Roman','Fontsize',5)
    xlabel('Tiempo','Fontname','Arial','Fontsize',20,'color','g')
    ylabel('Unidades','Fontname','Calibri','Fontsize',25,'rotatio',15)
    hold on
    plot(t1,a,'-og','LineWidth',1,'MarkerSize',1,'MarkerEdgeColor','r','MarkerFaceColor','b')
    plot(t1,V,'-.dc','LineWidth',1,'MarkerSize',1,'MarkerEdgeColor','r','MarkerFaceColor','b')
    legend('posición','Aceleración','Velocidad','location','best')
    grid minor
    %MAS
    subplot(1,3,3)
    plot(t,y,'m','LineWidth',1,'MarkerSize',1,'MarkerEdgeColor','r','MarkerFaceColor','b')
    hold on
    z=zeros(length(y));
    plot((0:t(end)),(z));
    title('Movimiento Armonico Simple','FontName','chiller','FontSize',20)
    ylabel('Desplazamiento','FontName','Algebrian','FontSize',40,'color',[0.3 0.5 0.7])
    xlabel('Tiempo','FontName','Times New Roman','Fontangle','italic','FontSize',30,'FontWeight','bold','color',[0.5 0.7 0.3])
    grid minor
%     ginput 

