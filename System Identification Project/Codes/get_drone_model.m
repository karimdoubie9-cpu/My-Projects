function sys_tf = get_drone_model(theta)
    % Unpack parameters
    [A,B,C,D] = uav_longitudinal_grey(theta, 1);
    sys=ss(A,B,C,D);
    sys_tf=[tf(sys(1, 1));tf(sys(2, 1))];
    % Xu = theta(1);
    % Xq = theta(2);
    % Mu = theta(3);
    % Mq = theta(4);
    % Xd = theta(5);
    % Md = theta(6);
    % 
    % % Fixed gravity
    % g = 9.81;
    % 
    % % State Space Matrices (from your code)
    % % States: [u; q; theta]
    % A = [Xu, Xq, -g; 
    %      Mu, Mq,  0; 
    %       0,  1,  0];
    % 
    % B = [Xd; Md; 0];
    % 
    % % Outputs defined in your code:
    % % 1: u (Velocity)
    % % 2: q (Pitch Rate)
    % % 3: theta (Angle)
    % % 4: ax (Acceleration = Xu*u + Xq*q + Xd*delta)
    % C = [ 0, 1, 0;
    %      Xu, Xq, 0]; 
    % 
    % D = [ 0; Xd];
    % 
    % % Create State Space System
    % sys_ss = ss(A, B, C, D);
    % % 
    % % Extract Transfer Function for the specific output requested
    % % (e.g., output_idx = 2 for Pitch Rate)
    % %sys_tf = tf(sys_ss(output_idx, 1)); 
    % sys_tf=[tf(sys_ss(1, 1));tf(sys_ss(2, 1))];
end