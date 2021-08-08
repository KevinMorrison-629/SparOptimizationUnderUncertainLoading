% minimize wing spar weight subject to stress constraints at manuever
clear all;
close all;
clc

global std_stress
global mean_stress

% carbon fiber values from http://www.performance-composites.com/carbonfibre/mechanicalproperties_2.asp
Nelem = 50;
L = 7.5; % semi-span in meters
rho = 1600; % density of standard carbon fiber, kg/m^3
yield = 600e6; % tensile strength of standard carbon fiber, Pa
E = 70e9; % Young's modulus, Pa
W = 0.5*500*9.8; % half of the operational weight, N
fnom = (2*(2.5*W)/(L^2))*[L:-L/Nelem:0].'; % loading at manueuver

nQuadPts = 7; % Number of Hermite Quadrature Points

% compute the displacements and the stresses
r_in = 0.0415*ones(1,Nelem+1).';
r_out = 0.05*ones(1,Nelem+1).';
Iyy = CalcSecondMomentAnnulus(r_in, r_out);

x = linspace(0,L,Nelem+1);
mu = 0;
sigman = @(n) fnom(1)/(10*n);

if nQuadPts == 1
    % n=1 point
    xi = [0];
    wts = [sqrt(pi)]./sqrt(pi);
elseif nQuadPts == 2
    % n=2 points
    xi = [0.707107 -0.707107];
    wts = [0.886227 0.886227]./sqrt(pi);
elseif nQuadPts == 3
    % n=3 points
    xi = [-1.22474487 0 1.22474487];
    wts = [0.295408975 1.1816359 0.295408975]./sqrt(pi);
elseif nQuadPts == 4
    % n=4 points
    xi = [-1.6506801 -0.52464762 0.52464762 1.6506801];
    wts = [0.081312835 0.80491409 0.80491409 0.081312835]./sqrt(pi);
elseif nQuadPts == 5
    % n=5 points
    xi = [0 -0.958572 0.958572 -2.02018 2.02018];
    wts = [0.945309 0.393619 0.393619 0.0199532 0.0199532]./sqrt(pi);
elseif nQuadPts == 6
    % n=6 points
    xi = [-2.3506050 -1.3358491 -0.43607741 0.43607741 1.3358491 2.3506050];
    wts = [0.0045300099 0.15706732 0.72462960 0.72462960 0.15706732 0.0045300099]./sqrt(pi);
elseif nQuadPts == 7
    % n=7 points
    xi = [-2.651961357, -1.673551629, -0.8162878829, 0, 0.816287883, 1.673551629, 2.651961357];
    wts = [9.71781245E-4, 0.05451558282, 0.4256072526, 0.8102646176, 0.4256072526, 0.0545155828, 9.71781245E-4]./sqrt(pi);
else
    % Error
end

mean_stress = 0;
mean_squared_stress = 0;
f = @(x,xi1,xi2,xi3,xi4) fnom.' + xi1*cos((2*1-1)*pi.*x./(2*L)) + xi2*cos((2*2-1)*pi.*x./(2*L)) + xi3*cos((2*3-1)*pi.*x./(2*L)) + xi4*cos((2*4-1)*pi.*x./(2*L));

for i1 = 1:size(xi,2)
    pt1 = sqrt(2)*sigman(1)*xi(i1) + mu;
    for i2 = 1:size(xi,2)
        pt2 = sqrt(2)*sigman(2)*xi(i2) + mu;
        for i3 = 1:size(xi,2)
            pt3 = sqrt(2)*sigman(3)*xi(i3) + mu;
            for i4 = 1:size(xi,2)
                pt4 = sqrt(2)*sigman(4)*xi(i4) + mu;
                mean_f = f(x,pt1,pt2,pt3,pt4);
                u = CalcBeamDisplacement(L, E, Iyy, mean_f, Nelem);
                stress = CalcBeamStress(L, E, r_out, u, Nelem);
                mean_stress = mean_stress + wts(i1)*wts(i2)*wts(i3)*wts(i4).*stress;
                mean_squared_stress = mean_squared_stress + wts(i1)*wts(i2)*wts(i3)*wts(i4).*(stress.^2);
            end
        end
    end
end

std_stress = sqrt(mean_squared_stress - (mean_stress.^2));

figure(1)
plot(x,mean_stress)
hold on
plot(x,mean_stress + 6*std_stress)
hold on
plot(x,mean_stress - 6*std_stress)

xlabel('Position along Spar (m)')
ylabel('Stress (Pa)')
legend('mean stress','mean stress - 6*\sigma','mean stress + 6*\sigma')
title('Nominal Stress Distribution Along Spar')



