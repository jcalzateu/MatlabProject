clc,close all,clear all;
%% Juan Camilo Alzate Urrego Cod. 316505
%% Taller 2 Informatica II
%% 17 de marzo de 2020
figure('name','COVID-19','NumberTitle','off','Color','w')
subplot(1,2,1)
x=[40 1 8 5 5 2 3 1 1 1 7 1];
c=categorical({'Bogotá D.C','Cundinamarca','Antioquica','Valle del Cauca','Bolivar','Atlántico','Norte de Santander','Santander','Caldas','Risaralda','Huila','Meta'});
bar(c,x,0.7,'FaceColor','r','EdgeColor','b','LineWidth',1.5)
title('17/Marzo/2020','FontSize',40,'FontName','Times New Roman','FontWeight','bold','color',[0.2 0.3 0.4],'Fontangle','italic')
xlabel('Ciudades','FontSize',25,'FontName','Gabriola','color','r')
ylabel('Número de Casos','FontSize',20,'FontName','chiller','color',[0.5 0.3 0.9],'rotation',25,'Fontangle','italic')
%18 de marzo de 2020
subplot(1,2,2)
y=[45 3 8 13 7 2 3 2 1 5 2 8 2 1];
explode=ones(size(y));
labels={'Bogotá D.C','Cundinamarca','Antioquica','Valle del Cauca','Bolivar','Atlántico','Norte de Santander','Santander','Caldas','Risaralda','Quindío','Huila','Tolima','Meta'};
pie(y,explode)
title('18/Marzo/2020','FontSize',25,'FontName','Tempus Sans ITC','FontWeight','bold','color','red')
legend(labels,'Location','southwest','Orientation','vertical')