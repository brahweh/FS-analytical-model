function [kappatotal, lambdaequiv, lambda0]=a20200325_BernardandAnderko_conductivity(T)
%constants
F = 96485.33212; %C/mol
Evac=8.854*10^-12;
kb = 1.381e-23;
R_gas=8.314;
e_val = 1.602*10^-19; %C, fundamental charge
NA = 6.022*10^23;

%inputs
% T=20;
T = 273.15 + T; %use Kelvin
z = [1 -1];
e = z*e_val;

%empirical fitting constants for lambda to D conversion [Matsuyama 2018]
q_D0 = [exp(20.2) exp(19.4)]; %1:Na 2:Cl
t_D0 = [1 0.95];

eta = waterviscosity(T);%kg/(m*s) or Pa*s

sigma = [0.98 1.81]*10^-10; %crystallographic ionic radii
sigma = sigma*2; %diameter

conc = [141 119.1]; %mol/m3 , number density of solutes
rho = conc*NA;


%% conductivity at infinite dilution calc

A = [-3.3594 -3.4051]; %1:Na 2:Cl
B = [75.492 216.03]; %1:Na 2:Cl

lambda0 = exp(A + B./T)./eta; %S cm^2/mol
lambda0 = lambda0/10^4; %S m^2 / mol

%% relaxation term 1

%building blocks

D0 = D0solve(lambda0,q_D0,t_D0,T);

[E0val,Einfval] = Ecalc(T);

epsilon = 4*pi*Evac*E0val;

sigmaij = 0;
for i=1:length(sigma)
    sigmaij = (sigmaij+sigma(i))/2;
end



%mid range equations
kappaq = sqrt(4*pi/(epsilon*kb*T)*(rho(1)*e(1)^2*D0(1)+rho(2)*e(2)^2*D0(2))/(D0(1)+D0(2)));

i0 = sinh(kappaq*sigmaij)/(kappaq*sigmaij);

i1 = cosh(kappaq*sigmaij)/(kappaq*sigmaij) - sinh(kappaq*sigmaij)/(kappaq*sigmaij)^2;


%integral equations and their support

alpha = sqrt(4*pi*e_val^2/(epsilon*kb*T));

Gammaiterate = 0;
Gamma=5.933000333579967e+08;
Gammaprev = 0;
while abs(Gamma - Gammaprev) > 1e0
    
    Gammaprev = Gamma;
    
    a = alpha^2./(2.*Gamma).*z./(1+Gamma.*sigma);
    
    Deltasum = 0;
    for i=1:length(z)
        Deltasum = Deltasum + rho(i)*sigma(i)^3;
    end
    Delta = 1 - pi/6*Deltasum;
    
    Omegasum = 0;
    for i=1:length(z)
        Omegasum = Omegasum + rho(i)*sigma(i)^3/(1+Gamma*sigma(i));
    end
    Omega = 1 + pi/(2*Delta)*Omegasum;
    
    Pnsum = 0;
    for i=1:length(z)
        Pnsum = Pnsum + rho(i)*sigma(i)*z(i)/(1+Gamma*sigma(i));
    end
    Pn = 1/Omega*Pnsum;
    
    Gammasum = 0;
    for i=1:length(z)
        Gammasum = Gammasum + rho(i)*((z(i)-pi/(2*Delta)*Pn*sigma(i)^2)/(1+Gamma*sigma(i)))^2;
    end
    Gamma = sqrt( 4*pi*e_val^2/(epsilon*kb*T)*Gammasum  / 4);
    
    Gammaiterate = Gammaiterate + 1;
end

Aij = e(1)*e(2)/(epsilon*kb*T*(1+Gamma*sigma(1))*(1+Gamma*sigma(2)));

I1sum = 0;
for i=1:length(z)
    I1sum = I1sum + rho(i)*a(i)^2*exp(-kappaq*sigmaij);
end
I1 = -kappaq*Aij*exp(-kappaq*sigmaij)/(kappaq^2 + 2*Gamma*kappaq + 2*Gamma^2 - (2*Gamma^2/alpha^2)*I1sum);

%main
relax1 = -kappaq^2/3*(i0-epsilon*kb*T/(e(1)*e(2))*kappaq*sigmaij^2*i1)*I1;

%% relaxation term 2

%supporting expressions

sigmabar1sum = 0;
for i=1:length(z)
    sigmabar1sum = sigmabar1sum + rho(i)*a(i)^2*sigma(i)^1;
end
sigmabar1 = 1/alpha^2*sigmabar1sum;

% sigmabar2sum = 0;
% for i=1:length(z)
%     sigmabar2sum = sigmabar2sum + rho(i)*a(i)^2*sigma(i)^2;
% end
% sigmabar2 = 1/alpha^2*sigmabar2sum;

% X = 2*Gamma*(1+Gamma*sigmabar1)/(1-Gamma^2*sigmabar2)
X = 2*Gamma*(1 + Gamma*sigmabar1);

Cij = cosh(kappaq*sigmaij) + X/kappaq*sinh(kappaq*sigmaij);

%Bij = e(1)*e(2)/(epsilon*kb*T*(1+Gamma*sigma(1))*(1+Gamma*sigma(2))*(1-Gamma^2*sigmabar2));


%main
relax2 = -kappaq^2/3*(i0 - epsilon*kb*T/(e(1)*e(2))*kappaq*sigmaij^2*i1) *...
    Aij^2*((X^2 + kappaq^2)/(X^2 - kappaq^2) * (kappaq^2/(4*X^2)*exp(2*X*sigmaij)*E1(2*X+kappaq,sigmaij)+...
    exp(-kappaq*sigmaij)/(4*X^2*sigmaij^2) * (1 + (2*X-kappaq)*sigmaij))+kappaq*Cij/(X - kappaq)*((X^2-2*kappaq^2)/(2*kappaq*(X+kappaq))*...
    exp(X*sigmaij) * E1(X+2*kappaq,sigmaij)-exp(-2*kappaq*sigmaij)/(2*kappaq*(X+kappaq)*sigmaij^2)*(1+X*sigmaij))+...
    (1+X*sigmaij)*((X^2-kappaq^2)^2/(4*X^2*kappaq^2)*exp(X*sigmaij)*E1(X+kappaq,sigmaij)-...
    X^2/(4*kappaq^2)*exp((X-kappaq)*sigmaij)*E1(X,sigmaij)*(1 + kappaq*sigmaij) +...
    exp(-kappaq*sigmaij)/(4*X^2*sigmaij^2)*(1+(X-kappaq)*sigmaij+X^3*sigmaij^2/kappaq)));


%% hydrodynamic relaxation

%main

hydro = -4*Gamma^2*Aij*kb*T/(48*pi*eta*(D0(1)+D0(2)))*(i0 - epsilon*kb*T/(e(1)*e(2))*kappaq*sigmaij^2*i1)*...
    ((1 + X*sigmaij + X^2*sigmaij^2/3)*(X^2/kappaq^2*(exp((X-kappaq)*sigmaij))*E1(X,sigmaij)*(1+...
    kappaq*sigmaij) - X^2*exp(-kappaq*sigmaij)/(kappaq*(X+kappaq)) - X^2/kappaq^2*E1(X+kappaq,sigmaij)+...
    ((2*X^2-kappaq^2)/X^2)*exp(X*sigmaij)*E1(X+kappaq,sigmaij) - exp(-kappaq*sigmaij)/(X^2*sigmaij^2)*(1+(X-...
    kappaq)*sigmaij) - X*exp(-kappaq*sigmaij)/(X+kappaq)) + exp(-kappaq*sigmaij)/(X^2*sigmaij^2)*(1+(2*X-kappaq)*sigmaij)-...
    (4*pi^2-kappaq^2)/X^2*exp(2*X*sigmaij)*E1(2*X+kappaq,sigmaij));

%% electrophoretic correction
v0 = zeros(1,length(z));
dv1 = zeros(1,length(z));
I=zeros(1,length(z));
J = zeros(1,length(z));
        
for ii=1:length(z)
    E=1;
    
    %denominator
    v0(ii) = e(ii)*E*D0(ii)/(kb*T);
    
    %first order
    dv1sum = 0;
    for i=1:length(z)
        dv1sum = dv1sum + rho(ii)*z(ii)*sigmaij^2;
    end
    dv1(ii) = -e(ii)*E/(3*pi*eta)*(Gamma/(1+Gamma*sigma(ii))+pi/(2*Delta)*...
        Pn*sigma(ii)/(z(ii)*(1+Gamma*sigma(ii)))+pi/z(ii)*dv1sum);
    
    %second order
    I(ii) = e(ii)*E*kappaq^2*Aij/(24*pi*eta*(X^2-kappaq^2)*sigmaij^2*(1+Gamma*sigma(1))*(1+Gamma*sigma(2)))*...
        ((1+2*X*sigmaij)-Cij*exp(-kappaq*sigmaij)*(1+(X-kappaq)*sigmaij)-...
        2*X^2*sigmaij^2*exp(2*X*sigmaij)*E1(2*X,sigmaij)+...
        Cij*(X^2+kappaq^2)*sigmaij^2*exp(X*sigmaij)*E1(X+kappaq,sigmaij));
    

    Jsum = 0;
    for i=1:length(z)
        Jsum = Jsum + rho(ii)*a(ii)^2*exp(-kappaq*sigma(ii));
    end
    J(ii) = e(ii)*E/(12*pi*eta*(1+Gamma*sigma(1))*(1+Gamma*sigma(2)))*(cosh(kappaq*sigmaij)-sinh(kappaq*sigmaij)/(kappaq*sigmaij))*...
        kappaq^2*exp(-kappaq*sigmaij)/(kappaq^2+2*Gamma*kappaq+2*Gamma^2-(2*Gamma^2)/alpha^2*Jsum);
end

%% final conductivity value

relaxtotal = relax1 + relax2 + hydro;

EPtotal = (dv1 + I + J)./v0;

lambdaequiv = lambda0.*(1+EPtotal).*(1+relaxtotal);

lambdatotal = sum(lambdaequiv); %S*m^2/mol

kappa = conc.*abs(z).*lambdaequiv;

kappatotal = sum(kappa); %S/m



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
% 
end
