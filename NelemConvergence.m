% minimize wing spar weight subject to stress constraints at manuever
clear all;
close all;
clc

global std_stress
global mean_stress

% carbon fiber values from http://www.performance-composites.com/carbonfibre/mechanicalproperties_2.asp
L = 7.5; % semi-span in meters
rho = 1600; % density of standard carbon fiber, kg/m^3
yield = 600e6; % tensile strength of standard carbon fiber, Pa
E = 70e9; % Young's modulus, Pa
W = 0.5*500*9.8; % half of the operational weight, N
nominal_mass = 29.32;

Nelem = 10;

data_fval = [];
data_opt_dvs = [];
mean_stresses = [];
stdevs_stresses = [];
data_tip_stress = [];
t_end = []; % comp times for each opt (with differing quad points)


for i = Nelem:10:60

Nelem = i;

f_nom = (2*(2.5*W)/(L^2))*[L:-L/Nelem:0].'; % loading at manueuver

% define function and constraints
fun = @(r) SparWeight(r, L, rho, Nelem);
lb = 0.01*ones(2*(Nelem+1),1);
up = 0.05*ones(2*(Nelem+1),1);
A = zeros(Nelem+1,2*(Nelem+1));
b = -0.0025*ones(Nelem+1,1);
for k = 1:(Nelem+1)
    A(k,k) = 1.0;
    A(k,Nelem+1+k) = -1.0;
end

% define initial guess (the nominal spar)
r0 = zeros(2*(Nelem+1),1);
r0(1:Nelem+1) = 0.0415*ones(Nelem+1,1);
r0(Nelem+2:2*(Nelem+1)) = 0.05*ones(Nelem+1,1);

options = optimset('GradObj','on','GradConstr','on', 'TolCon', 1e-4, ...
    'TolX', 1e-8, 'Display','iter','Algorithm','sqp'); %, 'DerivativeCheck','on');


% t_start = tic; % start timer for calculating computation time
nQuadPts = 3; % Number of Hermite Quadrature Points
nonlcon = @(r) WingConstraints(r, L, E, f_nom, yield, Nelem,nQuadPts);
[ropt,fval,exitflag,output] = fmincon(fun, r0, A, b, [], [], lb, up, ...
    nonlcon, options);

% plot optimal radii
r_in = ropt(1:Nelem+1);
r_out = ropt(Nelem+2:2*(Nelem+1));
x = [0:L/Nelem:L].';
figure(2*i - 1)
plot(x, r_in, '-ks');
hold on;
plot(x, r_out, '--ks');
title_str = 'Spar Geometry for Uncertain Loading with 3 Quadrature Points';
title(sprintf(title_str,i))
xlabel('Position Along Spar (m)')
ylabel('Position from Spar centerline (m)')

% display weight and stress constraints
[f,~] = fun(ropt);
[c,~,~,~] = nonlcon(ropt);

data_fval(i/10) = fval;
data_tip_stress(i/10) = real(mean_stress(end));
end

Nelem = 40;
L = 7.5;
nominal_mass = 29.32;
x = [0:L/Nelem:L].';

load('opt_Nelem_convergence_data')
n_Elems = 10:10:40;
figure(10)
plot(n_Elems,abs(data_fval-nominal_mass)/nominal_mass,'-o')
ylabel('Percent Spar Mass Improvement')
xlabel('Number of Spar Elements')
title('Effect of Number of Elements on Spar Mass Optimization')
ylim([0.7069,0.709])

figure(11)
plot(x,real(mean_stress))
hold on
plot(x,(real(mean_stress) + 6*real(std_stress)))
hold on
plot(x,(real(mean_stress) - 6*real(std_stress)))
title_str = 'Mean Spar Stress with Uncertainty for 3 Quadrature Points and 40 Elements';
title(sprintf(title_str,i))
xlabel('Position Along Spar (m)')
ylabel('Spar Stress (Pa)')
legend('mean stress','mean stress - 6*\sigma','mean stress + 6*\sigma','location','southwest')
ylim([-3.5e8, 7.5e8])

figure(12)
plot(x, ropt(1:Nelem+1), '-ks');
hold on;
plot(x, ropt(Nelem+2:end), '--ks');
title_str = 'Spar Geometry for Uncertain Loading with 3 Quadrature Points and 40 Elements';
title(sprintf(title_str,3))
xlabel('Position Along Spar (m)')
ylabel('Position from Spar centerline (m)')

%save opt_Nelem_convergence_data data_fval ropt mean_stress std_stress
