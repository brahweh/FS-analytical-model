close all
clearvars
clc


%% generate temperature function in 'x', add conductivity

%temperature model constants
alpha = 0.143E-6;
lambda_p = 84E-3;
lambda_f = 598E-3;
d = 10e-6;
del = (1.625e-3)*(1); %%tune this to match experimental behavior if nec
lu = 5000e-6; %1900 um is electroactive area according to FEM
ld = lu;
P = 20e-3 %heater power
b = 930e-6;
% b = 880e-6;
lh = 400e-6;
% A = 1.7e-5; %half area of shunt tubing

C_fact=.5; %correction coefficient for effect of reduced heat transfer in y (1 is no change)(k in write up)

% ydist=0.631; %temp of interest is the avg in y. we assume linear T distribution in y
ydist = 0.482; %see discussion of temperature weighted by current distribution in y

W=250e-6; %1/2 heater width
% W=465e-6;

T_amb = 20;

vel_start=10e-6; %m/s
vel_end=1000e-6;
vel_number=20; %number of points to create velocity data
xx=linspace(vel_start,vel_end,vel_number); %velocity space for use vs avg T
xx2=xx/6.69e-4*20; %velocity space in mL/hr
xx3=xx/6.69e-4*333; %uL/min
xx_int=1:length(xx);

T_avg_A=zeros(1,length(xx));
T_avg_B=zeros(1,length(xx));
T_avg_C=zeros(1,length(xx));
T_avg_total=zeros(1,length(xx));

d_US = 50;
d_heat = 50;
d_DS = 50;
d_tot = d_US + d_heat + d_DS;

%ionic conductivity constants
%Na:1, Cl:2, HCO3:3, H:4, K:5, Ca:6, Mg:7, PO4:8, H2PO4:9, HPO4:10
z = [1 -1 -1 1 1 2 2 3 -1 -2];
q = 1.602*10^-19; %C
F = 96485.33212; %C/mol
Evac=8.854*10^-12;
kb = 1.381e-23;
R_gas=8.314;

%CSF electrolyte concentrations: CSF:1, PBS(Corning 46-013):2
conc(:,1) = [141 119.1 21.5 10^-(7-3) 2.9 2.3/2 2.25/2 0.168 0 0]; %OG unit is mMol, convert to mol/m3
%denominators are to correct equivalency unit in source [Irani 2009]
conc(:,2) = [1.572 1.400 0 0 0.0445 0 0 0 0.0176 0.101]*10^2; %OG unit is mol/L, 1000x to mol/m3, 1/10x for dilution from 10x to 1x PBS

%i:velocity, j:x distance, k:ion
kappa=zeros(vel_number,d_tot);
for i=1:length(xx)
    
    %calculate T vs x for different velocities
    nu = xx(i);
    
    r1 = (nu*del + sqrt(del^2*nu^2+C_fact*8*alpha^2))/(2*alpha*del);
    r2 = (nu*del - sqrt(del^2*nu^2+C_fact*8*alpha^2))/(2*alpha*del);
    
    Th=ydist*del*P/(lambda_f*b*(C_fact*2*W+del^2*(r1-r2)));
    
    %T(W)-Th, T(BC1)=0 as BC1 approaches inf. this is Lammerink sol'n (correct)
    c1=Th/exp(r2*W);
    c2=0;
    
    x1 = linspace(-lu,-W,d_US);
    T1 = Th.*exp(r1.*(x1+W));
    
    x2 = linspace(-W,W,d_heat);
    T2 = Th*ones(1,d_heat);
    
    x3 = linspace(W,ld,d_DS);
    T3 = Th.*exp(r2.*(x3-W));
    
    x=[x1 x2 x3];
    T=[T1 T2 T3];
%     x=linspace(20,50,150);
%     T=linspace(0,30,150);
    
    figure(1)
    plot(x*10^3,T+20)
    hold on
    xlabel('Distance(mm)','FontSize',16)
    ylabel('T(^oC)','FontSize',16)
    xticks([-5 -2.5 0 2.5 5])
    yticks([20 21 22 23 24])
    
    T_abs=T+T_amb; %convert overheat temperature to absolute temp
    T(75)
    
    
    %weighted average of temperature upstream,above,downstream of heater
    %analytical expressions for the average
    T_avg_A(i) = Th*(1-exp(r1*(W-lu)))/(r1*(lu+ld)); %function avg x<-W
    T_avg_B(i) = 2*W*Th/(lu+ld); %function avg -W<x<W
    T_avg_C(i) = Th*(exp(r2*(ld-W))-1)/(r2*(lu+ld)); %function avg x>W
    
    T_avg_total(i) = T_avg_A(i)+T_avg_B(i)+T_avg_C(i);
    
    %convert fluid temperature to conductivity
    for j=1:d_tot
%         if j>d_US && j<d_US+d_heat %speed up for constant T above heater
%             kappa(i,j) = kappa(i,d_US);
%         else
            [kappatotal, lambdaequiv, lambda0]=a20200325_BernardandAnderko_conductivity(T_abs(j));
            kappa(i,j)=kappatotal;
%         end
    end
    
    figure(2)
    plot(x*10^3,kappa(i,:))
    hold on
    xlabel('Distance(mm)','FontSize',16)
    ylabel('Conductivity(S/m)','FontSize',16)
end

% writematrix(kappa(1,:), 'M.csv')

%normalize temperature in each volume block, first unweight each block
T_avg_A_norm = T_avg_A*(lu+ld)/(lu-W)+T_amb;
T_avg_A_norm = T_avg_A_norm/T_avg_A_norm(1);
T_avg_B_norm = T_avg_B*(lu+ld)/(2*W)+T_amb;
T_avg_B_norm = T_avg_B_norm/T_avg_B_norm(1);
T_avg_C_norm = T_avg_C*(lu+ld)/(ld-W)+T_amb;
T_avg_C_norm = T_avg_C_norm/T_avg_C_norm(1);
T_avg_total_norm = T_avg_total+T_amb;
T_avg_total_norm = T_avg_total_norm/T_avg_total_norm(1);



%average conductivity by distance block (upstream, heater, downstream)
for i=1:length(xx)
    kappa_avg_A(i)=mean(kappa(i,1:d_US));
    kappa_avg_B(i)=mean(kappa(i,d_US+1:d_US+d_heat));
    kappa_avg_C(i)=mean(kappa(i,d_US+d_heat+1:d_tot));
end

kappa_avg_total = (lu-W)/(lu+ld)*kappa_avg_A+(W+W)/(lu+ld)*kappa_avg_B+(ld-W)/(lu+ld)*kappa_avg_C;

%% calculate apparent impedance

[kappa_amb, lambdaequiv, lambda0]=a20200325_BernardandAnderko_conductivity(T_amb);

Z_app_amb = Za(kappa_amb); %scalar

Z_app_A = Za(kappa_avg_A);
Z_app_B = Za(kappa_avg_B);
Z_app_C = Za(kappa_avg_C);
Z_app = Za(kappa_avg_total); %vector length is velocity points

Z_norm_A = 100*(Z_app_A-Z_app_amb)/Z_app_amb;
Z_norm_B = 100*(Z_app_B-Z_app_amb)/Z_app_amb;
Z_norm_C = 100*(Z_app_C-Z_app_amb)/Z_app_amb;
Z_norm = 100*(Z_app-Z_app_amb)/Z_app_amb;


%% additional plots

%T vs v in 3 subplots in x (upstream, heater, downstream)
figure(4)
subplot(2,2,1)
plot(xx3,T_avg_total)
% ylim([0 8])
xlabel('Flow Rate(uL/min)','FontSize',14)
ylabel('T(^oC)','FontSize',14)
title('Weighted Average','FontSize',16)

subplot(2,2,2)
plot(xx3,T_avg_A*(lu+ld)/(lu-W))
xlabel('Flow Rate(uL/min)','FontSize',14)
ylabel('T(^oC)','FontSize',14)
title('Upstream Fluid','FontSize',16)

subplot(2,2,3)
plot(xx3,T_avg_B*(lu+ld)/(2*W))
xlabel('Flow Rate(uL/min)','FontSize',14)
ylabel('T(^oC)','FontSize',14)
title('Fluid Above Heater','FontSize',16)

subplot(2,2,4)
plot(xx3,T_avg_C*(lu+ld)/(ld-W))
xlabel('Flow Rate(uL/min)','FontSize',14)
ylabel('T(^oC)','FontSize',14)
title('Downstream Fluid','FontSize',16)

%resistivity
figure(5)
plot(x,Za(kappa(1,:)))
% xlabel('Temperature (deg C)')
% ylabel('Apparent Impedance (ohm)')
xticks([20 30 40 50])
yticks([4000 5000 6000 7000])

%normalized T and conductivity
figure(6)
subplot(2,2,1)
plot(xx3,T_avg_total_norm)
hold on
plot(xx3,Z_app/Z_app(1))
xlabel('Flow Rate(uL/min)','FontSize',14)
ylabel('AU','FontSize',14)
title('Weighted Average','FontSize',16)

subplot(2,2,2)
plot(xx3,T_avg_A_norm)
hold on
plot(xx3,Z_app_A/Z_app_A(1))
xlabel('Flow Rate(uL/min)','FontSize',14)
ylabel('AU','FontSize',14)
legend('Temperature','Apparent Impedance')
title('Upstream Fluid','FontSize',16)

subplot(2,2,3)
plot(xx3,T_avg_B_norm)
hold on
plot(xx3,Z_app_B/Z_app_B(1))
xlabel('Flow Rate(uL/min)','FontSize',14)
ylabel('AU','FontSize',14)
title('Fluid Above Heater','FontSize',16)

subplot(2,2,4)
plot(xx3,T_avg_C_norm)
hold on
plot(xx3,Z_app_C/Z_app_C(1))
xlabel('Flow Rate(uL/min)','FontSize',14)
ylabel('AU','FontSize',14)
title('Downstream Fluid','FontSize',16)




figure(8)
subplot(2,2,1)
plot(xx3,Z_norm)
title('total')
xlabel('Flow Rate(uL/min)')
ylabel('Normalized Impedance Drop')

subplot(2,2,2)
plot(xx3,Z_norm_A)
title('Upstream')
xlabel('Flow Rate(uL/min)')
ylabel('Normalized Impedance Drop')

subplot(2,2,3)
plot(xx3,Z_norm_B)
title('above heater')
xlabel('Flow Rate(uL/min)')
ylabel('Normalized Impedance Drop')

subplot(2,2,4)
plot(xx3,Z_norm_C)
title('downstream')
xlabel('Flow Rate(uL/min)')
ylabel('Normalized Impedance Drop')


writematrix([xx3' Z_norm' Z_norm_A' Z_norm_B' Z_norm_C'],'20 mW dump.xlsx')


%% helper functions
function etaval = waterviscosity(T)
%input T(K)

etaval=0;
aa=[280.68 511.45 61.131 .45903];
b=[-1.9 -7.7 -19.6 -40];

for i=1:4
    etaval=etaval+aa(i)*((T+0)/300)^b(i);
end
etaval=etaval*10^-6;
end

function [E0val,Einfval] = Ecalc(T)
%input T(K)

E0val=77.66-103.3*(1-300/(T+0));
Einfval=0.066*E0val;

end

function D0 = D0solve(lambda0,q,t,T) %outputs in weird units (1/10 of Scm^2/mol)
D0 = T*(lambda0./q).*(1./t);
end

function E1val = E1(gamma, sigmaij)

fun = @(r) exp(-gamma*r)./r;

E1val = integral(fun,sigmaij,1e-4);

end

function Zaval = Za(kappa) %calculate apparent impedance using Randalls circuit

Rs = 2*log(4)./(kappa*pi*sqrt(1.44e-8)); %area in m2

Zc = 1/(2*pi*50000*63*1.44e-4*10^-6); %cap/area of Pt elec, elec area in cm2, uF to F conversion

Zaval = Rs + 2*Zc;

end


