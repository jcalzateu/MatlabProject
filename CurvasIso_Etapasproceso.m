clc, close all, clear all;
tic
%% Curvas de isovelocidad
% Rxn : SO2 + 0.5 O2 ---> O3
%        A  + 0.5 B  --->  C  N2=I
% Datos
%% Reactor

P       = 2.0265;                            %[Atm]
To      = 683;                               %Kelvin
Fac_ef  = 1;                                 %Factor de efectividad
Di      = 10.7;                              %[m]
A       = pi*Di^2/4;                         %Área de sección transversal [m^2]

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

W_f     = 40000;                               %[Kg] --> Peso del Catalizador
rho_c   = 33.8*(1/2.2)*(1/0.3048^3);         %[kg/m^3] Catalizador
phi     = 0.43;                              % Porosidad
dp      = 0.013;                             % m

%% Propiedades termodinamica

Cp      = [29.637 0.034735 9.2903E-6 -2.9885E-8 1.0937E-11
        29.526 -8.8999E-3 3.8083E-5 -3.2629E-8 -8.8607E-12
        22.466 0.11981 -9.0842E-5 2.550e-8 -7.9208E-13
        29.342 -3.5395e-3 1.0076e-5 -4.3116e-9 2.5935e-13];          % KJ/Kmol~K  ---> Yaws
    
w       = [0.245 0.022 0.422 0.04];                                  % Factor acéntrico   ---> Yaws
Pc      = [77.809 49.7705 80.9968 33.4961];                          % [Atm]              ---> Yaws
Tc      = [430.75 154.58 490.85 126.1];                              % K                  ---> Yaws
coef    = [-1 -0.5 1 0];                                             % Coeficientes estequiométricos   
Hf      = [-296800 0 -395700 0];                                     % KJ/Kmol ---> Yaws
Gf      = [-300100 0 -371100 0];                                     % KJ/Kmol ---> Yaws
Sf      = (Gf-Hf)/298.15;                                            % KJ/Kmol~K---> Yaws
Hrxn_o  = sum(Hf.*coef);                                             % Kj/Kmol

%% Curvas de isovelocidad

T_final_Cur       = 800+273.15;                            %K  Temperatura para hacer las curvas de equilibrio
T_inicial_Cur     = 350+273.15;                            %K  Temperatura para hacer las curvas de equilibrio
T_Cur             = T_inicial_Cur:T_final_Cur;             %K

%% Condiciones para el método numérico

n  = 200;
v  = 1e-9;
xo = 0;

%% Solución 

for i=1:length(T_Cur)
[Conversion(i),iter,error,fun]=N_R1(@(X)fun_r1(X,Theta_B,Theta_I,P,T_Cur(i),FAo,0),xo,v,n);
end

for i=1:length(T_Cur)
[Conversion1(i),iter,error,fun]=N_R1(@(X)fun_r1(X,Theta_B,Theta_I,P,T_Cur(i),FAo,1e-6),xo,v,n);
end

for i=1:length(T_Cur)
[Conversion2(i),iter,error,fun]=N_R1(@(X)fun_r1(X,Theta_B,Theta_I,P,T_Cur(i),FAo,1e-5),xo,v,n);
end

for i=1:length(T_Cur)
[Conversion3(i),iter,error,fun]=N_R1(@(X)fun_r1(X,Theta_B,Theta_I,P,T_Cur(i),FAo,1e-4),xo,v,n);
end

for i=1:length(T_Cur)
[Conversion4(i),iter,error,fun]=N_R1(@(X)fun_r1(X,Theta_B,Theta_I,P,T_Cur(i),FAo,1e-3),xo,v,n);
end

for i=1:length(T_Cur)
[Conversion5(i),iter,error,fun]=N_R1(@(X)fun_r1(X,Theta_B,Theta_I,P,T_Cur(i),FAo,1e-2),xo,v,n);
end

for i=1:length(T_Cur)
[Conversion6(i),iter,error,fun]=N_R1(@(X)fun_r1(X,Theta_B,Theta_I,P,T_Cur(i),FAo,3.5e-2),xo,v,n);
end

for i=1:length(T_Cur)
[Conversion7(i),iter,error,fun]=N_R1(@(X)fun_r1(X,Theta_B,Theta_I,P,T_Cur(i),FAo,0.05),xo,v,n);
end

for i=1:length(T_Cur)
[Conversion8(i),iter,error,fun]=N_R1(@(X)fun_r1(X,Theta_B,Theta_I,P,T_Cur(i),FAo,8.5e-2),xo,v,n);
end

for i=1:length(T_Cur)
[Conversion9(i),iter,error,fun]=N_R1(@(X)fun_r1(X,Theta_B,Theta_I,P,T_Cur(i),FAo,0.1),xo,v,n);
end

for i=1:length(T_Cur)
[Conversion10(i),iter,error,fun]=N_R1(@(X)fun_r1(X,Theta_B,Theta_I,P,T_Cur(i),FAo,0.15),xo,v,n);
end

for i=1:length(T_Cur)
[Conversion11(i),iter,error,fun]=N_R1(@(X)fun_r1(X,Theta_B,Theta_I,P,T_Cur(i),FAo,0.3),xo,v,n);
end
% Parte real
Conversion1=real(Conversion1);
Conversion2=real(Conversion2);
Conversion3=real(Conversion3);
Conversion4=real(Conversion4);
Conversion5=real(Conversion5);
Conversion6=real(Conversion6);
Conversion7=real(Conversion7);
Conversion8=real(Conversion8);
Conversion9=real(Conversion9);
Conversion10=real(Conversion10);
Conversion11=real(Conversion11);

for i=1:length(Conversion)
    if Conversion1(i)<0  
        Conversion1(i)=NaN;
    elseif Conversion1(i) == Inf
        Conversion1(i)=NaN;
    end
end
for i=1:length(Conversion)
    if Conversion2(i)<0  
        Conversion2(i)=NaN;
    elseif Conversion2(i) == Inf
        Conversion2(i)=NaN;
    end
end
for i=1:length(Conversion)
    if Conversion3(i)<0  
        Conversion3(i)=NaN;
    elseif Conversion3(i) == Inf
        Conversion3(i)=NaN;
    end
end
for i=1:length(Conversion)
    if Conversion4(i)<0  
        Conversion4(i)=NaN;
    elseif Conversion4(i) == Inf
        Conversion4(i)=NaN;
    end
end
for i=1:length(Conversion)
    if Conversion5(i)<0  
        Conversion5(i)=NaN;
    elseif Conversion5(i) == Inf
        Conversion5(i)=NaN;
    end
end
for i=1:length(Conversion)
    if Conversion6(i)<0  
        Conversion6(i)=NaN;
    elseif Conversion6(i) == Inf
        Conversion6(i)=NaN;
    end
end
for i=1:length(Conversion)
    if Conversion7(i)<0  
        Conversion7(i)=NaN;
    elseif Conversion7(i) == Inf
        Conversion7(i)=NaN;
    end
end
for i=1:length(Conversion)
    if Conversion8(i)<0  
        Conversion8(i)=NaN;
    elseif Conversion8(i) == Inf
        Conversion8(i)=NaN;
    end
end
for i=1:length(Conversion)
    if Conversion9(i)<0  
        Conversion9(i)=NaN;
    elseif Conversion9(i) == Inf
        Conversion9(i)=NaN;
    end
end
for i=1:length(Conversion)
    if Conversion10(i)<0  
        Conversion10(i)=NaN;
    elseif Conversion10(i) == Inf
        Conversion10(i)=NaN;
    end
end
for i=1:length(Conversion)
    if Conversion11(i)<0  
        Conversion11(i)=NaN;
    elseif Conversion11(i) == Inf
        Conversion11(i)=NaN;
    end
end

%% Progresión óptima de temperatura

    Val_Max1 = max(Conversion1);
    posi     = find(Conversion1==Val_Max1);
    T_Max1   = T_Cur(posi);
    Val_Max2 = max(Conversion2);
    posi     = find(Conversion2==Val_Max2);
    T_Max2   = T_Cur(posi);
    Val_Max3 = max(Conversion3);
    posi     = find(Conversion3==Val_Max3);
    T_Max3   = T_Cur(posi);
    Val_Max4 = max(Conversion4);
    posi     = find(Conversion4==Val_Max4);
    T_Max4   = T_Cur(posi);
    Val_Max5 = max(Conversion5);
    posi     = find(Conversion5==Val_Max5);
    T_Max5   = T_Cur(posi);
    Val_Max6 = max(Conversion6);
    posi     = find(Conversion6==Val_Max6);
    T_Max6   = T_Cur(posi);
    Val_Max7 = max(Conversion7);
    posi     = find(Conversion7==Val_Max7);
    T_Max7   = T_Cur(posi);
    Val_Max8 = max(Conversion8);
    posi     = find(Conversion8==Val_Max8);
    T_Max8   = T_Cur(posi);
    Val_Max9 = max(Conversion9);
    posi     = find(Conversion9==Val_Max9);
    T_Max9   = T_Cur(posi);
    Val_Max10= max(Conversion10);
    posi     = find(Conversion10==Val_Max10);
    T_Max10  = T_Cur(posi);
    Val_Max11= max(Conversion11);
    posi     = find(Conversion11==Val_Max11);
    T_Max11  = T_Cur(posi);
    Val_Max  = [1 Val_Max1 Val_Max2 Val_Max3 Val_Max4 Val_Max5 Val_Max6 Val_Max7 Val_Max8 Val_Max9 Val_Max10 Val_Max11];
    T_Max    = [T_Cur(1) T_Max1 T_Max2 T_Max3 T_Max4 T_Max5 T_Max6 T_Max7 T_Max8 T_Max9 T_Max10 T_Max11];

%% Polinomio para la linea de POT
    [CR]  = polyfit(T_Max,Val_Max,5);
    CR_EQ = @(X) CR(1)*X^3+CR(2)*X^2+CR(3)*X+CR(4);
    XCR   = linspace(T_Cur(1),T_Cur(end),450);
    YCR   = polyval(CR,XCR);
    
    for i=1:length(YCR)
        if YCR(i)<0
        YCR(i) = NaN;
        end
    end
    
%% Reactor 1 Solución Para parada en POT

% Linea de operación adibatica
    % Ractor 1
        To1    = 400+273.15;                                                          %Kelvin
        [TF1]  = N_R1(@(Y) fun_OB(Y,P,Theta_T,To1,Cp,coef,0,CR),800,1e-12,400);
        TR1    = linspace(To1,TF1,100);
        for i=1:length(TR1)
            [X_BER1(i)] = fun_X_BE(TR1(i),P,Theta_T,To1,Cp,coef);
            [rA1(i)]    = fun_r1(X_BER1(i),Theta_B,Theta_I,P,TR1(i),FAo,0);
            rA1_INV(i)  = 1/rA1(i);
        end

%% Desde POT hasta punto b (regla de la palanca)

    [CR_r0] = polyfit(T_Cur,Conversion,5);
    XCR0    = linspace(T_Cur(1),T_Cur(end),450);
    YCR0    = polyval(CR_r0,XCR);
    [TF1_r0]=N_R1(@(Y) fun_OB(Y,P,Theta_T,To1,Cp,coef,0,CR_r0),800,1e-12,400);

    TR1_b   = linspace(TF1,TF1_r0,100);
    for i=1:length(TR1_b)
        [X_BER1_b(i)] = fun_X_BE(TR1_b(i),P,Theta_T,To1,Cp,coef);
        [rA1_b(i)]    = fun_r1(X_BER1_b(i),Theta_B,Theta_I,P,TR1_b(i),FAo,0);
        rA1_INV_b(i)  = 1/rA1_b(i);
    end
    
%% Metodo de los trapecios
    % Etapa 1
        I1 = 0;
        h  = (X_BER1(end)-X_BER1(1))/length(X_BER1);
        for i=1:length(X_BER1)-1
            I1 = (h*(rA1_INV(i)+rA1_INV(i+1))/2)+I1;
        end
        error = 1;
        tol   = 1e-1;
        I1_b  = 0;
        h     = (X_BER1_b(end)-X_BER1_b(1))/(10);

            for i=1:length(X_BER1_b)-1
                I1_b =(h*(rA1_INV_b(i)+rA1_INV_b(i+1))/2)+I1_b;
                error=abs(abs(I1_b)-abs(I1))/abs(I1_b);
                if error<tol
                    break
            break
                end
            end

        X_BER1_b = X_BER1_b(1:i);
        TR1_b    = TR1_b(1:i);
        rA1_b    = rA1_b(1:i);
    
%% Enfriamiento   

    FAo1=FAo*(1-X_BER1_b(end));
    FBo1=FAo*(Theta_B-0.5*X_BER1_b(end));
    FCo1=FAo*X_BER1_b(end);
    FIo1=FIo;
    %% Curva de corte enfriamiento
    
        for i=1:length(T_Cur)
        [Conversion_etapa1(i),iter,error,fun]=N_R1(@(X)fun_r1(X,Theta_B,Theta_I,P,T_Cur(i),FAo,rA1_b(end)),xo,v,n);
        end
Conversion_etapa1=real(Conversion_etapa1);
        for i=1:length(Conversion_etapa1)
        if Conversion_etapa1(i)<0  
            Conversion_etapa1(i)=NaN;
        elseif Conversion_etapa1(i) == Inf
            Conversion_etapa1(i)=NaN;
        end
        end
[T2o]=N_R1(@(Y) fun_r1_OBT(Y,X_BER1_b(end),Theta_B,Theta_I,P,FAo,rA1_b(end)),728.15,1e-2,50);
%% 

% Etapa 2 
    Ec_Etaps = polyfit(TR1,X_BER1,1);
    Ec_Etaps(2) = X_BER1_b(end)-T2o*Ec_Etaps(1);
    Theta_T2 = [FAo1/FAo1 FBo1/FAo1 FCo1/FAo1 FIo1/FAo1];

    [TF2]  = N_R1(@(Y) fun_objEtap2(Y,CR,Ec_Etaps),600,1e-8,400);
    TR2      = linspace(T2o,TF2,450);
    X_BER2   = polyval(Ec_Etaps,TR2);
    
    
%% Metodo de la palanca
    [CR_r0] = polyfit(T_Cur,Conversion,5);
    XCR0    = linspace(T_Cur(1),T_Cur(end),450);
    YCR0    = polyval(CR_r0,XCR);
    [TF2_r0]=N_R1(@(Y) fun_OB(Y,P,Theta_T2,720,Cp,coef,0,CR_r0),700,1e-2,50);
    
    TR2_b      = linspace(TF2,TF2_r0,100);
    X_BER2_b   = polyval(Ec_Etaps,TR2_b);
   
        for i=1:length(TR2)
            [rA2(i)]    = fun_r1(X_BER2(i),Theta_T2(2),Theta_T2(4),P,TR2(i),FAo1,0);
            rA2_INV(i)  = 1/rA2(i);
        end
    
    for i=1:length(TR2_b)
        [rA2_b(i)]    = fun_r1(X_BER2_b(i),Theta_T2(2),Theta_T2(4),P,TR2_b(i),FAo1,0);
        rA2_INV_b(i)  = 1/rA2_b(i);
    end
    
    %% Metodo de los trapacios
    % Etapa 2
        I1 = 0;
        h  = (X_BER2(end)-X_BER2(1))/length(X_BER2);
        for i=1:length(X_BER2)-1
            I1 = (h*(rA2_INV(i)+rA2_INV(i+1))/2)+I1;
        end
        error = 1;
        tol   = 1e-1;
        I1_b  = 0;
        h     = (X_BER2_b(end)-X_BER2_b(1))/(100);

            for i=1:length(X_BER1_b)-1
                I1_b =(h*(rA1_INV_b(i)+rA1_INV_b(i+1))/2)+I1_b;
                error=abs(abs(I1_b)-abs(I1))/abs(I1_b);
                if error<tol
                    break
            break
                end
            end

        X_BER2_b = X_BER2_b(1:i);
        TR2_b    = TR2_b(1:i);
        rA2_b    = rA2_b(1:i);
%% curva de enfriamiento

%% Enfriamiento   

    FAo2=FAo1*(1-X_BER2_b(end));
    FBo2=FAo1*(Theta_T2(2)-0.5*X_BER2_b(end));
    FCo2=FAo1*X_BER2_b(end);
    FIo2=FIo;
    %% Curva de corte enfriamiento
    
        for i=1:length(T_Cur)
        [Conversion_etapa2(i),iter,error,fun]=N_R1(@(X)fun_r1(X,Theta_T2(2),Theta_T2(4),P,T_Cur(i),FAo1,rA2_b(end)),xo,v,n);
        end
Conversion_etapa2=real(Conversion_etapa2);
        for i=1:length(Conversion_etapa2)
        if Conversion_etapa2(i)<0  
            Conversion_etapa2(i)=NaN;
        elseif Conversion_etapa2(i) == Inf
            Conversion_etapa2(i)=NaN;
        end
        end
[T3o]=N_R1(@(Y) fun_r1_OBT(Y,X_BER2_b(end),Theta_T2(2),Theta_T2(4),P,FAo2,rA2_b(end)),710,1e-2,50);

% Etapa 3
    Ec_Etaps = polyfit(TR2,X_BER2,1);
    Ec_Etaps(2) = X_BER2_b(end)-T3o*Ec_Etaps(1);
    Theta_T2 = [FAo2/FAo2 FBo2/FAo2 FCo2/FAo2 FIo2/FAo2];

    [TF3]  = N_R1(@(Y) fun_objEtap2(Y,CR,Ec_Etaps),750,1e-8,400);
    TR3      = linspace(T3o,TF3,450);
    X_BER3   = polyval(Ec_Etaps,TR3);
%% Metodo de la palanca
    [CR_r0] = polyfit(T_Cur,Conversion,5);
    XCR0    = linspace(T_Cur(1),T_Cur(end),450);
    YCR0    = polyval(CR_r0,XCR);
    [TF3_r0]=N_R1(@(Y) fun_OB(Y,P,Theta_T2,715,Cp,coef,0,CR_r0),600,1e-2,50);
    
    TR3_b      = linspace(TF3,TF3_r0,100);
    X_BER3_b   = polyval(Ec_Etaps,TR3_b);
   
        for i=1:length(TR3)
            [rA3(i)]    = fun_r1(X_BER2(i),Theta_T2(2),Theta_T2(4),P,TR3(i),FAo1,0);
            rA3_INV(i)  = 1/rA2(i);
        end
    
    for i=1:length(TR2_b)
        [rA3_b(i)]    = fun_r1(X_BER2_b(i),Theta_T2(2),Theta_T2(4),P,TR2_b(i),FAo1,0);
        rA3_INV_b(i)  = 1/rA2_b(i);
    end
    
    %% Metodo de los trapacios
    % Etapa 3
        I1 = 0;
        h  = (X_BER3(end)-X_BER3(1))/length(X_BER3);
        for i=1:length(X_BER3)-1
            I1 = (h*(rA3_INV(i)+rA3_INV(i+1))/2)+I1;
        end
        error = 1;
        tol   = 1e-1;
        I1_b  = 0;
        h     = (X_BER3_b(end)-X_BER3_b(1))/(100);

            for i=1:length(X_BER1_b)-1
                I1_b =(h*(rA1_INV_b(i)+rA1_INV_b(i+1))/2)+I1_b;
                error=abs(abs(I1_b)-abs(I1))/abs(I1_b);
                if error<tol
                    break
            break
                end
            end

        X_BER3_b = X_BER3_b(1:i);
        TR3_b    = TR3_b(1:i);
        rA3_b    = rA3_b(1:i);

 clc;
    [sal,sal1,sal2,sal3]=Fun_llamar(1);
    disp(sal);
    disp(sal1);
    disp(sal2);
    disp(sal3);
    disp('__________________________________________________________________');
    disp('RESULTADOS');
sal = ['La temperatura a la entrada de la etapa 1 es:  ' num2str(To1) ' Kelvin']; disp(sal);
sal = ['La temperatura a la salida de la etapa  1 es:  ' num2str(TR1_b(end)) ' Kelvin']; disp(sal);
sal = ['La conversión  a la salida de la etapa  1 es:  ' num2str(X_BER1_b(end))]; disp(sal); 
sal = ['La temperatura a la entrada de la etapa 2 es:  ' num2str(T2o) ' Kelvin']; disp(sal);
sal = ['La temperatura a la salida del la etapa 2 es:  ' num2str(TR2_b(end)) ' Kelvin']; disp(sal); 
sal = ['La conversión  a la salida de la etapa  2 es:  ' num2str(X_BER2_b(end))]; disp(sal);
sal = ['La temperatura a la entrada de la etapa 3 es:  ' num2str(T3o) ' Kelvin']; disp(sal);
sal = ['La temperatura a la salida del la etapa 3 es:  ' num2str(TR3_b(end)) ' Kelvin']; disp(sal); 
sal = ['La conversión  a la salida de la etapa  3 es:  ' num2str(X_BER3_b(end))]; disp(sal);

%% Gráficas
figure('name','Curvas de isovelocidad','NumberTitle','off','Color','w')
plot(T_Cur,Conversion)
hold on
plot(T_Cur,Conversion1)
plot(T_Cur,Conversion2)
plot(T_Cur,Conversion3)
plot(T_Cur,Conversion4)
plot(T_Cur,Conversion5)
plot(T_Cur,Conversion6)
plot(T_Cur,Conversion7)
plot(T_Cur,Conversion8)
plot(T_Cur,Conversion9)
plot(T_Cur,Conversion10)
plot(T_Cur,Conversion11)
plot(T_Max,Val_Max,'--g')
a = [0.3 0.52];
b = [0.62 0.63];
annotation('textarrow',a,b,'String','POT','Color','black','Fontsize',12)
title('$$ Curvas \: de \: isovelocidad \:  $$','Interpreter','latex', 'FontSize' , 12,'Fontname','Times new roman')
xlabel('$$ Temperatura \: \left[K\right] $$','Interpreter','latex', 'FontSize' , 12,'Fontname','Times new roman')
ylabel('$$ Conversion  $$','Interpreter','latex', 'FontSize' , 12,'Fontname','Times new roman')
grid minor
legend('r=0','r=1e-6','r=1e-5','r=1e-4','r=1e-3','r=1e-2','r=3.5e-2','r=0.05','r=0.085','r=0.1','r=0.15','r=0.3','location','best');
%% POT
    figure('name','Línea POT','NumberTitle','off','Color','w')
    P1 = plot(T_Max,Val_Max,'--g');
    hold on
    P2 = plot(T_Cur,Conversion);
    P4 = plot(XCR,YCR);
%% Linea de operación etapa 1
    plot(TR1,X_BER1,'color','black')
    plot(TR1_b,X_BER1_b,'color','black')
    plot(T_Cur(100:250),Conversion_etapa1(100:250),'--r');
    P3 = plot([T2o TR1_b(end)],[X_BER1_b(end) X_BER1_b(end)],'--black');
    a = [0.25 0.41];
    b = [0.5 0.5];
    annotation('textarrow',a,b,'String','Etapa 1','Color','black','Fontsize',10)
    a = [0.3 0.41];
    b = [0.6 0.74];
    annotation('textarrow',a,b,'String','Enfriamiento','Color','black','Fontsize',10)
%% Linea de operación etapa 2
    plot(TR2,X_BER2,'color','black');
    plot(TR2_b,X_BER2_b,'color','black');
    plot([T3o TR2_b(end)],[X_BER2_b(end) X_BER2_b(end)],'--black')
    plot(T_Cur(80:150),Conversion_etapa2(80:150),'--r');
    a = [0.25 0.35];
    b = [0.8 0.8];
    annotation('textarrow',a,b,'String','Etapa 2','Color','black','Fontsize',10)
    grid minor
    title('$$ Linea \: de \: Operacion \: adiabatica \:  $$','Interpreter','latex', 'FontSize' , 12,'Fontname','Times new roman')
    xlabel('$$ Temperatura \: \left[K\right] $$','Interpreter','latex', 'FontSize' , 12,'Fontname','Times new roman')
    ylabel('$$ Conversion  $$','Interpreter','latex', 'FontSize' , 12,'Fontname','Times new roman')
%     legend('POT','r=0','regresion','location','best')
%% Linea de operación etapa 3
    plot(TR3,X_BER3,'color','black')
    plot(TR3_b,X_BER3_b,'color','black');
    a = [0.25 0.3];
    b = [0.9 0.9];
    annotation('textarrow',a,b,'String','Etapa 3','Color','black','Fontsize',10)
    legend([P1 P2 P3 P4],{'POT' 'r=0' 'Enfriamiento' 'regresión'},'location','best')
    toc
    
%% Funciones 
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
Kp=10.^((5022./T)-4.765); 
end

% Ley de velocidad 1
function [r1]=fun_r1(X,Theta_B,Theta_I,P,T,FAo,R_Val)
    K1=fun_K1(T);
    K2=fun_K2(T);
    K3=fun_K3(T);
    Kp=fun_Kp(T);
    FT=FAo*(1+Theta_B+Theta_I-((1/2)*X));         %[lbmol/h]
    FA=FAo*(1-X);                                 %[lbmol/h]
    FB=FAo*(Theta_B-((1/2)*X));                   %[lbmol/h]
    FC=FAo*X;                                     %[lbmol/h]
    r1=-R_Val+(K1.*(P/FT).^2*FB.*FA.*(1-((P.*FC)./FT)./((P./FT).^(3/2).*FA.*(FB).^(1/2).*Kp)))./(22.414.*(1+P./FT.*(FA.*K2+FC.*K3))^2); 
end

%% Método

function[x1,iter,error,fun]=N_R1(f,xo,v,n)
% f:función
% xo:estimado inicial
% v:paso
% n:#max de iteraciones
iter=1;
error=1;
tol=1e-5;
while error>tol
    iter=iter+1;
    fun=feval(f,xo);
    fun_aum=feval(f,xo+v);
    der=(fun_aum-fun)/v;
    x1=xo-fun/der;
    error=abs((x1-xo)/x1);
    xo=x1;
    if iter>n
        break
    end
end
end

% Linea de operación adiabatica

function [XB_E]=fun_OB(Y,P,Theta_T,To,Cp,coef,X_menos1,CR)
    T=Y(1);
    w       = [0.245 0.022 0.422 0.04];                                  % Factor acéntrico   ---> Yaws
    Pc      = [78.84 50.43 82.07 33.94];                                 % [Bar]              ---> Yaws
    Tc      = [430.75 154.58 490.85 126.1];                              % K                  ---> Yaws  
    Hf      = [-296800 0 -395700 0];                                     % KJ/Kmol ---> Yaws
    Gf      = [-300100 0 -371100 0];                                     % KJ/Kmol ---> Yaws
    Sf      = (Gf-Hf)/298.15;                                            % KJ/Kmol~K---> Yaws
    Hrxn     =  H_T(T,P,Cp,w,Tc,Pc,coef,Hf); 
    
for i=1:length(coef)
        fCp  = @(x)(Cp(i,1)+Cp(i,2)*x+Cp(i,3)*x.^2+Cp(i,4)*x.^3+Cp(i,5)*x.^4);
         Cp_Val(i) = feval(fCp,T);                      % KJ/Kmol~K
         Cp_Valo(i)= feval(fCp,To);
%         Cp_Val(i) = integral(fCp,298.15,T);                      % KJ/Kmol~K
end
 Valorpol=polyval(CR,T);
 ValorBalance= ((sum(Cp_Val.*Theta_T)*T)/(-Hrxn))-(sum(Cp_Valo.*Theta_T)*To)/(-Hrxn);
    XB_E = ValorBalance-Valorpol-X_menos1;
end

function [r1]=fun_r1_OBT(Y,X,Theta_B,Theta_I,P,FAo,R_Val)
    T=Y(1);
    K1=fun_K1(T);
    K2=fun_K2(T);
    K3=fun_K3(T);
    Kp=fun_Kp(T);
    FT=FAo*(1+Theta_B+Theta_I-((1/2)*X));         %[lbmol/h]
    FA=FAo*(1-X);                                 %[lbmol/h]
    FB=FAo*(Theta_B-((1/2)*X));                   %[lbmol/h]
    FC=FAo*X;                                     %[lbmol/h]
    r1=-R_Val+(K1.*(P/FT).^2*FB.*FA.*(1-((P.*FC)./FT)./((P./FT).^(3/2).*FA.*(FB).^(1/2).*Kp)))./(22.414.*(1+P./FT.*(FA.*K2+FC.*K3))^2); 
end

function [T_etap]=fun_objEtap2(Y,CR,Ec_Etaps)
T=Y(1);
T_etap=polyval(CR,T)-polyval(Ec_Etaps,T);
end

function [XB_E]=fun_X_BE(T,P,Theta_T,To,Cp,coef)
    w       = [0.245 0.022 0.422 0.04];                                  % Factor acéntrico   ---> Yaws
    Pc      = [78.84 50.43 82.07 33.94];                                 % [Bar]              ---> Yaws
    Tc      = [430.75 154.58 490.85 126.1];                              % K                  ---> Yaws  
    Hf      = [-296800 0 -395700 0];                                     % KJ/Kmol ---> Yaws
    Gf      = [-300100 0 -371100 0];                                     % KJ/Kmol ---> Yaws
    Sf      = (Gf-Hf)/298.15;                                            % KJ/Kmol~K---> Yaws
    Hrxn     =  H_T(T,P,Cp,w,Tc,Pc,coef,Hf); 
    
for i=1:length(coef)
        fCp  = @(x)(Cp(i,1)+Cp(i,2)*x+Cp(i,3)*x.^2+Cp(i,4)*x.^3+Cp(i,5)*x.^4);
         Cp_Val(i) = feval(fCp,T);                      % KJ/Kmol~K
         Cp_Valo(i)= feval(fCp,To);
 end
    XB_E = ((sum(Cp_Val.*Theta_T)*T)/(-Hrxn))-(sum(Cp_Valo.*Theta_T)*To)/(-Hrxn);
end

function [sal,sal1,sal2,sal3]=Fun_llamar(A)
sal  = 'Grupo 1';
sal1 = 'Juan Camilo Alzate Urrego Cod:316505';
sal2 = 'Erika Alejandra Pardo Acosta Cod:316554';
sal3 = 'Estefania Rojas Zuleta Cod:315557';
end

%% Correción de la entalpia
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
