clear all
close 
clc

%% DATA

%Spacecraft
Ix = 4.208e11*10^(-3-6);        %Original data in g*mm^2, then converted in kg/m^2
Iy = 4.922e12*10^(-3-6);        %Inertias taken from CAD, from the fully deployed configuration.
Iz = 4.399e12*10^(-3-6);        %Iz major axis of rotation  
I = diag([Ix,Iy,Iz]);

%Sensors
Mount = [0 0 1 0 1 0];          %this is the mounting vector for the earth sensor
                                % the first 3 numbers identify the pitch axis of
                                % the sensor and the last 3 the roll axis, theese
                                % have different FOV limits, which are defined in
                                % the simulink

%Actuators (TEMP)
alpha = pi/6;
l=1.2   ;
A_act=[l*sin(alpha) 0 0;0 l*cos(alpha) 0;0 0 (0.25*sin(alpha)-0.25*cos(alpha))];
B=[1 1 -1 -1 ;-1 1 1 -1;1 -1 1 -1];
C=A_act*B;
q=null(C);
q=abs(q(:,1));
u_c = 1;

%Geometry (TEMP)

N_b =  [0,0,+1; 0,0,-1; 0,0,1; 0,0,-1;1,0,0 ; -1,0,0; 0,1,0; 0,-1,0; 0,0,1; 0,0,-1];
rho_s =  [0.5; 0.5; 0.5; 0.5; 0.5; 0.5; 0.1; 0.1; 0.1; 0.1];
rho_d = 0.1.*ones(10,1);
Area_vect = [3;3;3;3;1;1;1;1;1;1];
R_f = [2.5,0,0; 2.5,0,0; -2.5,0,0; -2.5,0,0; 1,0,0; -1,0,0; 0,1,0; 0,-1,0; 0,0,1; 0,0,-1];

n_face = length(rho_s);
c_g= 76.595*10^-3;
R_f_disp = [0,c_g,0];
for i=1:10
    R_f(i,:)=R_f(i,:)-R_f_disp;
end

%Orbit

a = 6899750;
e = 0.0005895;
i = deg2rad(5.0097);
OM = deg2rad(308.7725);
om = deg2rad(280.9595);
e_vers = [cos(i)*cos(OM) ; cos(i)*sin(OM) ; sin(i)];
G = 6.67e-11;
M = 5.97219e24;
mu = G*M;
p = a * (1 - e^2);
n = sqrt(mu/a^3);
theta0 = 0;
R_sun = 149597870700;

%Enviroment

    %MAG

    w_earth = 2*pi/(24*3600);
    J_b = [1 ; 18 ; 5]; % random orientation, magnitude assuming class 2 s/c
    
    [g, h] = IGRF13;

    %DRAG
    rho = 3.13e-13; %kg/m^3
    C_d = 2.2; % chosen assuming flat plate, range is 1.5 to 2.6
    
    rotOM = [cos(OM) sin(OM) 0; -sin(OM) cos(OM) 0; 0 0 1];
    roti  = [1 0 0; 0 cos(i) sin(i); 0 -sin(i) cos(i)];
    rotom = [cos(om) sin(om) 0; -sin(om) cos(om) 0; 0 0 1];
    Rot = rotOM * roti * rotom;

    %GG
    Ky = (Iz-Iy)/Ix;
    Kr = (Iz-Ix)/Iy;
    Kp = (Iy-Ix)/Iz;

    if Iy>Iz && Iz>Ix
        fprintf("The Gravity Gradient has a stabilizing effect on the Z axis, that is the minimum inertia axis")
    end

%%
%Preliminary Maximum disturbances
q_pan= 0.25; %assumed to be a non coated surface that has an average between 20% and 30%
q_body= 0.9; %took an extreme case of a unpolished aluminium body 
c_ps= 0; %the panels are symmetrical to each others
c_pa= c_ps;
c_g= 76.595*10^-3;
I_M= Iy;
I_m= Ix;
A_pan= (3+3); %area of two solar panel faces
A_body= 1; %area of one face of the body 
P_SRP= 4.54*10^-6; %solar radiation oressure for 1AU 
B_max= 2.47 * 10^-5; %maximum magnetic strength in Tesla
D_max= norm(J_b); %value taken from NASA report SP-8018, assumed Class 2


h_orbit = sqrt(mu*a*(1-e^2));
V_drag = mu/h_orbit*(1+e);
T_max_GG = ((3*mu)/(2*a^3)) * abs(I_M - I_m);
T_max_SRP = P_SRP * A_pan * ((1+q_pan)*(c_ps - c_g)) + P_SRP * A_body * ((1+q_body)*(c_ps - c_g)); %first part is the torque created by the panels and the second by the body  
T_max_mag = D_max * B_max;
T_max_Drag = 0.5*rho*V_drag^2*(A_pan+A_body)*C_d*(c_pa - c_g);

%% Regulators data

    %Detumbling
    k_det = 110;

    %Slew
    K_slew = -[0.3,0.3,0.3,0.2,0.2,0.2].*50;
    K_slew = diag(K_slew);
    filter_tau = 500;
    filter_a = 6/filter_tau;

    %Pointing
    tau = 200;
    


%% Manoeuvres

iterations = 100;

%% Detumbling
%set_param('Detumbling', 'FastRestart', 'on');



for jj = 1:iterations

    j =rand(1)*0.001;   %starting value randomization
    f =rand(1)*0.2;     %Main starting rotation
    k =rand(1)*0.001;

    w0x = f;
    w0y = j;
    w0z = k;

    j=pi/20*rand(1);
    A=[1,0,0;0,cos(j),sin(j);0,-sin(j),cos(j)];
    A=A*[1,0,0;0,cos(i),sin(i);0,-sin(i),cos(i)];
    w_LN = [0 0 n];

    out_det = sim("Detumbling.slx");
    w_out_det(:,jj) = out_det.w_out_det(end,:);
    w_det_residual(jj) = rad2deg(norm(out_det.w_out_det(end,:)));
    convergence_time(jj) = length(out_det.w_out_det)*0.1;
    theta0_det(jj) = out_det.theta_det(end); 
    DCM_out_det(:,:,jj) = out_det.real_DCM_det(:,:,end);

    w0x = w_out_det(1,jj);
    w0y = w_out_det(2,jj);
    w0z = w_out_det(3,jj);

    theta0 = theta0_det(jj);

    A = zeros(3,3);
    A = A +  DCM_out_det(:,:,jj);


    out_slew = sim("slew_manneouvre.slx");
    DCM_out_slew(:,:,jj) = out_slew.real_DCM_slew(:,:,end);
    A_LN_target(:,:,jj) = out_slew.A_LN_target(:,:,end);
    w_out_slew(:,jj) = out_slew.w_out_slew(end,:);
    theta0_slew(jj) = out_slew.theta0_slew(end);
    A_e_slew(:,:,jj) = out_slew.A_e(:,:,end);
end

%set_param('Detumbling', 'FastRestart', 'off');
%% stat analysis of det

% Calculate the mean
residual_w_mean = mean(w_det_residual)
mean_time = mean(convergence_time)
[max_residual_w] = max(w_det_residual)

% Calculate the variance
residual_w_variance = var(w_det_residual)
maximum_time = max(convergence_time);

close all
figure;
histogram(w_det_residual, 'Normalization', 'probability', 'BinWidth', 0.05);
title('Occurrence of Values');
xlabel('Value');
ylabel('Frequency');
grid on;
axis tight

%% stat analysis of slew

for jj=1:iterations
    eu_angl(jj) = rad2deg(acos((trace(DCM_out_slew(:,:,jj)*A_LN_target(:,:,jj)')-1)/2));
    %A_e = DCM_out_slew(:,:,jj)*A_LN_target(:,:,jj)';
    x_angl(jj) = abs(A_e_slew(2,3,jj));
    y_angl(jj) = abs(A_e_slew(3,1,jj));
    z_angl(jj) = abs(A_e_slew(1,2,jj));
end

% Calculate the mean
eu_angl_mean = mean(eu_angl)

% Calculate the variance
eu_angl_variance = var(eu_angl)
maximum_euler_slew = max(eu_angl)


%% Tracking

jj = randi([1 length(theta0_slew)]);

w0x = w_out_slew(1,jj);
w0y = w_out_slew(2,jj);
w0z = w_out_slew(3,jj);

theta0 = theta0_slew(jj);

A = zeros(3,3);
A = A + DCM_out_slew(:,:,jj);

out_tracking = sim("Tracking.slx");
A_LN_tracking = out_tracking.A_LN;
DCM_out_tracking = out_tracking.real_DCM_tracking;

%% stat analysis tracking

for jj=1000:length(A_LN_tracking)
    eu_angl_track(jj-999) = rad2deg(acos((trace(DCM_out_tracking(:,:,jj)*A_LN_tracking(:,:,jj)')-1)/2));
    A_e = DCM_out_tracking(:,:,jj)*A_LN_tracking(:,:,jj)';
    x_angl(jj-999) = rad2deg(abs(A_e(2,3)));
    y_angl(jj-999) = rad2deg(abs(A_e(3,1)));
    z_angl(jj-999) = rad2deg(abs(A_e(1,2)));
end

% Calculate the mean
x_angl_mean = (mean(x_angl));
y_angl_mean = (mean(y_angl));
z_angl_mean = (mean(z_angl));
% Calculate the variance
x_angl_variance = var(x_angl);
y_angl_variance = var(y_angl);
z_angl_variance = var(z_angl);


maximum_euler_ = max(eu_angl_track);
max_x_angl = max(x_angl) ;
max_y_angl = max(y_angl) ;
max_z_angl = max(z_angl) ;
frequency_over_1 = sum(eu_angl_track(:)>1)/length(eu_angl_track)
%%
% Plot the occurrence of the values

close all
figure;
histogram(x_angl, 'Normalization', 'probability', 'BinWidth', 0.001);
title('Occurrence of Values');
xlabel('Value');
ylabel('Frequency');
grid on;
axis tight

figure;
histogram(y_angl, 'Normalization', 'probability', 'BinWidth', 0.001);
title('Occurrence of Values');
xlabel('Value');
ylabel('Frequency');
grid on;
axis tight

figure;
histogram(z_angl, 'Normalization', 'probability', 'BinWidth', 0.001);
title('Occurrence of Values');
xlabel('Value');
ylabel('Frequency');
grid on;
axis tight
