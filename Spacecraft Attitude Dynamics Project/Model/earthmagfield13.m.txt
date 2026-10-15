function [B_N]  = earthmagfield13(r_N, t, g, h, alpha_G_0, n)
%
% Earth Magnetic Field with Spherical Harmonics
%
%DESCRIPTION:
%This code computes the components of the Earth's Magnetic Field using
%spherical harmonics from order 1 until order 13.
%
%Gaussian Coefficients from IGRF 13th Gen. 2020 are transformed from the
%Schmidt Quasi-Normalization to Gaussian Normalization using recursive
%formulas [1].
%
%--------------------------------------------------------------------------
% INPUTS:
%   r_N        [3x1]       S/C Pos. (Inertial Frame) [km]
%   t          [1x1]       Day Time                  [sec]
%   g          [14x13]     Gauss Coff. g_nm Matrix   [T]
%   h          [14x13]     Gauss Coff. h_nm Matrix   [T]
%   alpha_G_0  [1x1]       Greenwich Init. Longitude [rad]
%   n          [1x1]       Order Desired             [-]
%--------------------------------------------------------------------------
% OUTPUTS:
%   B_N        [3x1]       Mag. Field Components     [T]
%--------------------------------------------------------------------------
%
%NOTES:
% - Matrices containing coefficients g_nm and h_nm are contained in script
%   called "IGRF13.m"
% - In order to define the Model Order (n) tke into account that the lower 
%   is the orbit altitude the higher has to be the Model Order, therefore:
%     n = 1: Dipole Model (accurate enough only for orbits with altitudes above 20000 Km)
%     n < 4: for medium altitude orbits
%     n > 4: for low altitude orbits
% - Magnetic Field components and spacecraft position are in Earth centered
%   inertial frame
% - The Gaussian Normalization enables to save up to 7% of computational
%   effort [1]
%{
%NOTES ABOUT CELESTIAL COORDINATES
%The whole script works with just the s/c position (given by r_N).
%The Wertz notation defines:
% - phi = Longitude from Greenwich Meridian
% - alpa_G = Greenwich Merdian Longitude (note that this angle is updated
%            while orbiting in order to take into account the Earth
%            rotation)
% - delta = Elevation (but it is named also Declination) and here, due to
%           the fact that we do not take into account the Earth oblateness,
%           it is also the Latitude
% - theta = Co-Elevation (or Co-Declination) (note that in the global model
%           this variable is used for True Anomaly, but here, for the only
%           reason of being in accordance with Wertz's notation, it has
%           been renamed
% - alpha = Right Ascension, it is the angle between the projection of the
%           s/c position on the Equatorial Plane and the Vernal Axis (which
%           we assume to be the x-axis of the Inertial Ref. Frame);
%}
%
%CALLED FUNCTIONS:
% (none)
%
%UPDATES:
% (none)
%
%REFERENCES:
% [1] James R. Wertz, "Spacecraft Attitude Determination and Control". 
%     Kluwer Academic Publishers, Spinger, 1993.
%
%AUTHOR(s):
%Luigi De Maria, 2020
%
%% Pre-Settings & Memory Allocation
%Matrix of Coefficients S_n,m
S_nm = zeros(n+1,n);
S_00 = 1;
%Matrix of Coefficients P^n,m
P_nm = zeros(n+1,n);
P_00 = 1;
%Matrix of Coefficients Derivatives d_P^n,m
d_P_nm = zeros(n+1,n);
d_P_00 = 0;
%Matrix of Coefficients K^n,m
K_nm = zeros(n+1,n);
%Matrix of Coefficients g_nm, h_nm
g_nm = zeros(n+1,n);
h_nm = zeros(n+1,n);
%% Topocentric Horizon Angles
z = r_N(3,1);
x = r_N(1,1);
y = r_N(2,1);
ab = norm(r_N);
%Earth Radius [km]
Re = 6371;
%Earth Angular Velocity [rad/sec]
om_E = 7.292124*1e-5;
%Elevation (or Declination) = Latitude (no J2 effect)
delta = asin(z/ab);
%Coelevation
theta = pi/2 - delta;
%Greenwich Meridian Longitude
alpha_G = alpha_G_0 + t * om_E;
%Right Ascension (Local Sidereal Time)
if (y/ab) > 0
    alpha = acos((x/ab)/cos(delta));
else
    alpha = 2*pi - acos((x/ab)/cos(delta));
end
%Longitude (from Greenwich)
phi = alpha - alpha_G;
%% Gauss Normalization
%Coefficients S_nm
%OBS i = m and j = n
%First Row (S_n,0)
for j = 1 : n
    if j == 1
        S_nm(1,j) = S_00 * ((2*j - 1)/j);
    elseif j >= 2
        S_nm(1,j) = S_nm(1,j-1)*((2*j - 1)/j);
    end
end
%Remaining Components (S_n,m)
for j = 1 : n
    for i = 2 : j+1
            S_nm(i,j) = S_nm(i-1,j)*sqrt(((kronDel(1,i-1) + 1)*(j - (i-1) + 1))...
                                    /(j + (i -1)));
    end
end
%Coefficients K^n,m
for j = 1 : n
    for i = 1 : j+1
        if j == 1
            K_nm(i,j) = 0;
        elseif j > 1
            K_nm(i,j) = ((j - 1)^2 - (i-1)^2) / ((2*j - 1)*(2*j - 3));
        end
    end
end
%Coefficients P^n,m
P_nm(1,1) = cos(theta)*P_00;        %P_1,0
P_nm(2,1) = sin(theta)*P_00;        %P_1,1
for j = 2 : n
    for i = 1 : j+1
        if i-1 == j                             %P_n,n
            P_nm(i,j) = sin(theta)*P_nm(j,j-1);
        else
            if j < 3                            %P_n,m
                P_nm(i,j) = cos(theta)*P_nm(i,j-1);
            else
                P_nm(i,j) = cos(theta)*P_nm(i,j-1) - K_nm(i,j)*P_nm(i,j-2);
            end
        end
    end
end
%Derivative of d_P^n,m (wrt theta)
d_P_nm(1,1) = cos(theta)*d_P_00 - sin(theta)*P_00; %d_P_1,0
d_P_nm(2,1) = sin(theta)*d_P_00 + cos(theta)*P_00; %d_P_1,1
for j = 2 : n
    for i = 1 : j+1
        if i-1 == j
            d_P_nm(i,j) = sin(theta)*d_P_nm(j,j-1) + cos(theta)*P_nm(j,j-1);
        else
            if j < 3
                d_P_nm(i,j) = cos(theta)*d_P_nm(i,j-1) - sin(theta)*P_nm(i,j-1);
            else
                d_P_nm(i,j) = cos(theta)*d_P_nm(i,j-1) - ...
                        sin(theta)*P_nm(i,j-1) - K_nm(i,j)*d_P_nm(i,j-2);
            end
        end
    end
end
%% Guassian Coefficients Normalization
for j = 1 : n
    for i = 1 : j+1
        g_nm(i,j) = S_nm(i,j)*g(i,j);
        h_nm(i,j) = S_nm(i,j)*h(i,j);
    end
end
%% Magnetic Field Components - Spherical Coordinates
%Reset Handle Variables
dx = 0;
pr = 0;
%Radial Component
for j = 1 : n
    sx = (Re/norm(r_N))^(j+2) * (j + 1);
    for i = 1 : j+1
        dx = dx + (g_nm(i,j)*cos((i-1)*phi) + h_nm(i,j)*sin((i-1)*phi)) * P_nm(i,j);
    end
    pr = pr + sx*dx;
    dx = 0;
end
B_r = pr;
%Reset Handle Variables
dx = 0;
pr = 0;
%Coelevation Component
for j = 1 : n
    sx = ((Re/norm(r_N))^(j+2));
    for i = 1 : j+1
        dx = dx + (g_nm(i,j)*cos((i-1)*phi) + h_nm(i,j)*sin((i-1)*phi)) * d_P_nm(i,j);
    end
    pr = pr + sx*dx;
    dx = 0;
end
B_th = -pr;
%Reset Handle Variables
dx = 0;
pr = 0;
%Azimuth Component
for j = 1 : n
    sx = ((Re/norm(r_N))^(j+2));
    for i = 1 : j+1
        dx = dx + (i*(-g_nm(i,j)*sin((i-1)*phi) + h_nm(i,j)*cos((i-1)*phi)) * P_nm(i,j));
    end
    pr = pr + sx*dx;
    dx = 0;
end
B_phi = -1/sin(theta) * pr;
%% Magnetic Field Components - Inertial Ref. Frame
B_1 = (B_r*cos(delta) + B_th*sin(delta))*cos(alpha) - B_phi*sin(alpha);
B_2 = (B_r*cos(delta) + B_th*sin(delta))*sin(alpha) + B_phi*cos(alpha);
B_3 = (B_r*sin(delta) - B_th*cos(delta));
B_N = [B_1; B_2; B_3];
end
%% Kronecker Delta Function
%
%Do to the unavailability of the function "kroneckerDelta" in Simulink, we
%had to implement it.
%
function d = kronDel(j,k)
if j == k
    d = 1;
else
    d = 0;
end
end