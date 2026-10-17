function [T_final, rho_initial, rho_final] = adiabatic_expansion(P_initial, P_final, T_initial, gamma, R)
    % ADIABATIC_EXPANSION Computes final temperature and densities during adiabatic expansion.
    %
    %   [T_final, rho_initial, rho_final] = adiabatic_expansion(P_initial, P_final, T_initial, gamma, R)
    %
    %   Inputs:
    %       P_initial  - Initial pressure (Pa)
    %       P_final    - Final pressure (Pa)
    %       T_initial  - Initial temperature (K)
    %       gamma      - Heat capacity ratio (Cp/Cv)
    %       R          - Specific gas constant (J/kg·K)
    %
    %   Outputs:
    %       T_final    - Final temperature (K)
    %       rho_initial - Initial density (kg/m^3)
    %       rho_final   - Final density (kg/m^3)
    %
    % Compute final temperature using the adiabatic relation
    T_final = T_initial * (P_final / P_initial)^((gamma - 1) / gamma);
    
    % Compute initial density using the ideal gas law: rho = P / (R * T)
    rho_initial = P_initial / (R * T_initial);
    
    % Compute final density using the ideal gas law
    rho_final = P_final / (R * T_final);
end
