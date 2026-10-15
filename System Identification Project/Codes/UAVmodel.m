function [A,B,C,D] = UAVmodel(Xu,Xq,Mu,Mq,Xd,Md,sample_time)
% UAVMODEL for SANT-X Quadrotor - Course 2025/2026
% States: [u, q, theta]
% Inputs: [delta_lon]
% Outputs: [q, ax]

    g = 9.81; % Gravity acceleration [m/s^2]

    % State-transition matrix (A)
    A = [Xu,  Xq, -g;
         Mu,  Mq,  0;
         0,   1,   0];

    % Input-to-state matrix (B)
    B = [Xd; 
         Md; 
         0];

    % Output matrix (C) - SIMO (Single Input Multiple Outputs)
    % Row 1: q (pitch rate) -> [0 1 0]
    % Row 2: ax (longitudinal acceleration) -> [Xu Xq -g]
    % Note: ax is defined as u_dot
    C = [0,  1,  0;   
         Xu, Xq, 0]; 

    % Direct feedthrough matrix (D)
    % q does not depend on the input (0)
    % ax depends on the input through Xd
    D = [0;           
         Xd];
end