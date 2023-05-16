close all; clear all; clc;
tic
%% Reacción catalítica de SO2 a SO3 en presencia del catalizador V2O5
% Rxn: SO2(g) + 1/2 O2(g) --> SO3(g)
%        A    + 1/2 B     --> C         N2=I
% Datos
%% Reactor

P       = 2.0265;                            %[Bar]
To      = 680;                               %Kelvin
Fac_ef  = 1;                                 %Factor de efectividad
Di      = 10.7;                              %[m]
A       = pi*Di^2/4;                         %Área de sección transversal [m^2]
coef    = [-1 -0.5 1 0];
%% Fracción molar

x_A     = 8.36/100;                           
x_B     = 8.9/100;
x_I     = 82.74/100;
x_C     = 0;
Fracc   = [x_A x_B x_I x_C];
%% Peso molecular

PMA     = 64.066;                            %[Kg/Kmol]
PMB     = 16;                                %[Kg/Kmol]
PMC     = 80.066;                            %[Kg/Kmol]
PMI     = 28.013;                            %[Kg/Kmol]
PM      = [PMA PMB PMC PMI];                 %[Kg/Kmol]
%% Flujos de alimentación 

FT      = 5000;                              %[Kmol/h]
m       = FT*(x_A*PMA+x_B*PMB+x_I*PMI);      %[Kg/h]
FAo     = FT*x_A;                            %[Kmol/h]
FBo     = FT*x_B;                            %[Kmol/h]
FCo     = FT*x_C;                            %[Kmol/h]
FIo     = FT*x_I;                            %[Kmol/h]

%% Relaciones de entrada

Theta_B = FBo/FAo;
Theta_I = FIo/FAo;
Theta_A = FAo/FAo;
Theta_C = FCo/FAo;
Theta_T = [Theta_A Theta_B Theta_C Theta_I];
%% Propiedades del catalizador

W_f     = 1.5e4;                               %[Kg] --> Peso del Catalizador
rho_c   = 33.8*(1/2.2)*(1/0.3048^3);         %[kg/m^3] Catalizador
phi     = 0.43;                              % Porosidad
dp      = 0.013;                             % m
%% Propiedades termodinamica

Cp      = [29.637 0.034735 9.2903E-6 -2.9885E-8 1.0937E-11
        29.526 -8.8999E-3 3.8083E-5 -3.2629E-8 -8.8607E-12
        22.466 0.11981 -9.0842E-5 2.5503e-8 -7.9208E-13
        29.342 -3.5395e-3 1.0076e-5 -4.3116e-9 2.5935e-13];          % KJ/Kmol~K  ---> Yaws

%% Temperatura

T_final_Cur       = 800+273.15;                            %K  Temperatura para hacer las curvas de equilibrio
T_inicial_Cur     = 600+273.15;                            %K  Temperatura para hacer las curvas de equilibrio
T_Cur             = T_inicial_Cur:T_final_Cur;             %K  
%% Viscosidad de cada sustancia
 
    mu_Ec=[-11.103 0.502 -1.08E-4
        44.224 0.562  -1.13E-4
        -12.039 0.543 -1.6E-4
        42.606 0.475 -9.88E-5];                     %[micropoise]
    
%% SOLUCIÓN

    % Balance molar y de energía
    [W_sol,X_sol]=ode45(@(W,Y) fun_PFR_ODE(W,Y,Fac_ef,FAo,Theta_B,Theta_I,Cp,coef,m,A,dp,phi,rho_c,FIo,PM,T_Cur(1),P),linspace(0,W_f,10000),[0 To P]);
    FA=FAo*(1-X_sol(:,1));
    FB   = FAo*(Theta_B-((1/2)*X_sol(:,1)));                              %[Kmol/h]
    FC   = FAo*X_sol(:,1);
    [sal,sal1,sal2,sal3]=Fun_llamar(1);
    disp(sal);
    disp(sal1);
    disp(sal2);
    disp(sal3);
    disp('__________________________________________________________________');
    disp('RESULTADOS');
    sal = ['La conversión máxima es: ' num2str(X_sol(end,1))]; disp(sal); 
    sal = ['La temperatura máxima es: ' num2str(X_sol(end,2)) 'Kelvin']; disp(sal);
    sal = ['El flujo de A es: ' num2str(FA(end)) 'Kmol/h']; disp(sal);
    sal = ['El flujo de B es: ' num2str(FB(end)) 'Kmol/h']; disp(sal);
    sal = ['El flujo de C es: ' num2str(FC(end)) 'Kmol/h']; disp(sal);
    %% Graficas
figure('name','Balances','NumberTitle','off','color','w')
[ax]=plotyy(W_sol,X_sol(:,1),W_sol,X_sol(:,2));
title('Temperatura y conversión')
xlabel('Masa del catalizador [Kg]','FontSize',12,'FontName','Times New Roman')
set(get(ax(1),'YLabel'), 'String', 'Conversión','FontSize',12,'FontName','Times New Roman')
set(get(ax(2),'YLabel'), 'String', 'Temperatura [K]','FontSize',12,'FontName','Times New Roman')
grid minor
legend('Conversión','Temperatura','location','best')
figure('name','Flujos','NumberTitle','off','color','w')
plot(W_sol,FA,W_sol,FB,W_sol,FC)
grid minor
xlabel('Masa del catalizador [Kg]','FontSize',12,'FontName','Times New Roman')
ylabel('$$ Flujo \left[ \frac{Kmol}{h}\right]  $$','Interpreter','latex', 'FontSize' , 12,'Fontname','Times new roman')
legend('SO_2','O_2','SO_3','location','best')
toc

%% FUNCIONES

% Calculo de la entalpia de reacción
    function [H]=H_T(T,P,Cp,w,Tc,Pc,coef,Hf)
%% Calculo de la entalpia de la reacción estado estandar
Hrxn=sum(coef.*Hf);    %[j/mol]
%% Correción por temperatura
Hrxn_T=0;
DeltaCp1=0;
for i=1:length(coef)
    FCp=@(x) (Cp(i,1)+Cp(i,2)*x+Cp(i,3)*x.^2+Cp(i,4)*x.^3+Cp(i,5)*x.^4);
    ICp=integral(FCp,298.15,T); %[j/mol]
    Hrxn_T=Hrxn_T+(coef(i)*ICp);   %[j/mol]
end
%% Correción por presión
    R=8.31447e-3;                               %
    Hrxn_tp=0;
for i=1:length(coef)
        k=0.37464+1.54226*w(i)-0.26992*(w(i))^2;
        alpha=(1+k*(1-sqrt(T/Tc(i))))^2;
        a=0.45724*((R^2*Tc(i)^2)/(Pc(i)))*alpha;
        b=0.07780*(R*Tc(i))/(Pc(i));       
        A=(a*P)/(R^2*T^2);
        B=(b*P)/(R*T);
        da_dt=-0.45724*((R^2*(Tc(i))^2)/Pc(i))*k*sqrt(alpha/(T*Tc(i)));
        Z=[1,-(1-B),(A-(3*B^2)-(2*B)),-(A*B-B^2-B^3)];
        Zr=roots(Z);
    for j=1:length(coef)-1
        if isreal(Zr(j))==0
            Zr(j)=0;
        end
    end
        Zp=max(Zr); 
        Y=(R*T*(Zp-1)+((((T*da_dt)-a)/(2*sqrt(2)*b))*(log((Zp+(1+sqrt(2))*B)/(Zp+(1-sqrt(2))*B)))));
        Hrxn_tp=Hrxn_tp+Y*coef(i);
end
    H=(Hrxn+Hrxn_T+Hrxn_tp)*(1/abs(coef(1)));
    end

% Constante cinética 1
function [K1]=fun_K1(T)
K1=exp(12.160-(5473./T)); 
end
% Constante cinética 2
function [K2]=fun_K2(T)
K2=exp(-9.953+(8619./T)); 
end
% Constante cinética 3
function [K3]=fun_K3(T)
K3=exp(-71.745+(52596./T));
end
% Constante de equilibrio
function [Kp]=fun_Kp(T) 
 Kp = exp(11300/T-10.68);
end

%% Ley de velocidad
function [r1]=fun_r1(X,Theta_B,Theta_I,P,T,FAo)
    K1 = fun_K1(T);
    K2 = fun_K2(T);
    K3 = fun_K3(T);
    Kp = fun_Kp(T);
    FT=FAo*(1+Theta_B+Theta_I-((1/2)*X));                        %[lbmol/h]
    FA=FAo*(1-X);                                                %[lbmol/h]
    FB=FAo*(Theta_B-((1/2)*X));                                  %[lbmol/h]
    FC=FAo*X;                                                    %[lbmol/h]
    r1=-(K1*(P*FB/FT)*(P*FA/FT)*(1-((P*FC/FT))/((P*FA/FT)*(P*FB/FT)^(1/2)*Kp)))/(22.414*(1+(K2*P*FA/FT)+(K3*P*FC/FT))^2);
end

%% Balances

function [dYdw]=fun_PFR_ODE(~,Y,Fac_ef,FAo,Theta_B,Theta_I,Cp,coef,m,A,dp,phi,rho_c,FI,PM,To,Po)
    X=Y(1); T=Y(2); P=Y(3);
    %% Propieades termodinamicas 
        w       = [0.245 0.022 0.422 0.04];                                  % Factor acéntrico   ---> Yaws
        Pc      = [78.84 50.43 82.07 33.94];                                 % [Bar]              ---> Yaws
        Tc      = [430.75 154.58 490.85 126.1];                              % K                  ---> Yaws  
        Hf      = [-296800 0 -395700 0];                                     % KJ/Kmol ---> Yaws
        Gf      = [-300100 0 -371100 0];                                     % KJ/Kmol ---> Yaws
        Sf      = (Gf-Hf)/298.15;                                            % KJ/Kmol~K---> Yaws
        Hrxn     =  H_T(T,P,Cp,w,Tc,Pc,coef,Hf);
        r1   = fun_r1(X,Theta_B,Theta_I,P,T,FAo);
        FT   = FAo*(1+Theta_B+Theta_I-((1/2)*X));                    %[Kmol/h]
        FA   = FAo*(1-X);                                            %[Kmol/h]
        FB   = FAo*(Theta_B-((1/2)*X));                              %[Kmol/h]
        FC   = FAo*X;                                               %[Kmol/h]
    
    %% Viscosidad de cada sustancia y mezcla
        x_A     = FA/FT;                           
        x_B     = FB/FT;
        x_I     = FI/FT;
        x_C     = FC/FT;
        Fracc   = [x_A x_B x_I x_C];

        mu_Ec=[-11.103 0.502 -1.08E-4
            44.224 0.562  -1.13E-4
            -12.039 0.543 -1.6E-4
            42.606 0.475 -9.88E-5];                     %[micropoise]
        for j=1:length(coef)
            mu_fun =@(x) mu_Ec(j,1)+mu_Ec(j,2)*x+mu_Ec(j,3)*x^2;
            mu_eval(j) = feval(mu_fun,T);
        end
        % Viscosidad de la mezcla
        for i=1:length(coef)
            for j=1:length(coef)
                    if j ~= i
                        F_I(i,j)=(1+(mu_eval(i)/mu_eval(j))^0.5*(PM(j)/PM(i))^(1/4))^2/(8*(1+PM(i)/PM(j) ))^0.5;
                    else
                        F_I(i,j)=1;
                    end
            end
        end

        mu_mix_Eval = 0;
            for j=1:length(coef)
                mu_mix_Eval = (Fracc(j)*mu_eval(j))/(sum(Fracc.*F_I(j,:))) + mu_mix_Eval; % [micropoise]
            end
        mu_mix_Eval=mu_mix_Eval*1e-7;                        %% Kg/m~s
      
    %% Densidad 
    R_rho  = 0.08205;                               %[atm*m3/Kmol~K]   
    PM_mix = sum(PM.*Fracc);
    rho_o    = (P/(PM_mix*R_rho*T));
    % Balance molar
    dXdw = (r1*Fac_ef)/(-FAo);                      
   % Evaluación de las capacidades calorificas 
        for i=1:length(coef)
            fCp  = @(x)(Cp(i,1)+Cp(i,2)*x+Cp(i,3)*x.^2+Cp(i,4)*x.^3+Cp(i,5)*x.^4);
            Cp_Val(i) = feval(fCp,T);                      % KJ/Kmol~K
        end
    
    dTdw = (-r1*(-Hrxn))/(FA*Cp_Val(1)+FB*Cp_Val(2)+FC*Cp_Val(3)+FI*Cp_Val(4));    %---> Balance de energia

%    dPdw = -(m/A)*(1/dp)*((1-phi)/(phi^3))*((1-0.5*X)/(rho_o))*(Po/P)*(T/To)*((150*mu*(1-phi)/dp)+(1.75*m/A))*(1/(rho_c*(1-phi)*A));   %---> Caida de presión 
    dPdw = 0;
    dYdw = [dXdw dTdw dPdw]';
end
function [sal,sal1,sal2,sal3]=Fun_llamar(A)
sal  = 'Grupo 1';
sal1 = 'Juan Camilo Alzate Urrego Cod:316505';
sal2 = 'Erika Alejandra Pardo Acosta Cod:316554';
sal3 = 'Estefania Rojas Zuleta Cod:315557';
end