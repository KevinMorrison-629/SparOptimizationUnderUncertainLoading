% minimize wing spar weight subject to stress constraints at manuever
clear all;
close all;
clc

global std_stress
global mean_stress

% carbon fiber values from http://www.performance-composites.com/carbonfibre/mechanicalproperties_2.asp
Nelem = 15;
L = 7.5; % semi-span in meters
rho = 1600; % density of standard carbon fiber, kg/m^3
yield = 600e6; % tensile strength of standard carbon fiber, Pa
E = 70e9; % Young's modulus, Pa
W = 0.5*500*9.8; % half of the operational weight, N
f_nom = (2*(2.5*W)/(L^2))*[L:-L/Nelem:0].'; % loading at manueuver
nominal_mass = 29.32;
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

data_fval = [];
data_opt_dvs = [];
mean_stresses = [];
stdevs_stresses = [];
t_end = []; % comp times for each opt (with differing quad points)

for i = 1:7 % up to 7 quadrature points
t_start = tic; % start timer for calculating computation time
nQuadPts = i; % Number of Hermite Quadrature Points
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
title_str = 'Spar Geometry for Uncertain Loading with %d Quadrature Points';
title(sprintf(title_str,i))
xlabel('Position Along Spar (m)')
ylabel('Position from Spar centerline (m)')
dim = [.6 .5 .3 .3];
str = sprintf('mass = %d',fval);
annotation('textbox',dim,'String',str,'FitBoxToText','on');


% display weight and stress constraints
[f,~] = fun(ropt);
[c,~,~,~] = nonlcon(ropt);

% Plot optimal stress
figure(2*i)
plot(x,real(mean_stress))
hold on
plot(x,(real(mean_stress) + 6*real(std_stress)))
hold on
plot(x,(real(mean_stress) - 6*real(std_stress)))
title_str = 'Mean Spar Stress with Uncertainty for %d Quadrature Points';
title(sprintf(title_str,i))
xlabel('Position Along Spar (m)')
ylabel('Spar Stress (Pa)')
legend('mean stress','mean stress - 6*\sigma','mean stress + 6*\sigma')

data_fval(i) = f;
data_opt_dvs(:,i) = ropt;
mean_stresses(:,i) = real(mean_stress);
stdevs_stresses(:,i) = real(std_stress);

t_end(i) = toc(t_start);

end

% Plot effect of quad points on computation time and percent improvement
n_quadPts = 1:1:length(t_end);
figure(25)
yyaxis left
plot(n_quadPts(2:end),t_end(2:end))
title('Effect of Number of Quadrature Points on Computation Time')
xlabel('Number of Quadrature Points')
ylabel('Computation Time (seconds)')
yyaxis right
plot(n_quadPts(2:end),abs(data_fval(2:end)-nominal_mass)/nominal_mass)
ylabel('Percent Spar Mass Improvement')
