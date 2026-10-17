function w = gasAdiabaticVelocity(P0, P, T0, k, Rstar)
    % gasAdiabaticVelocity Computes the speed of a perfect gas flow.
    %
    %   w = perfect_gas_speed(P0, P, T0, k, Rstar) computes the speed
    %   using the given parameters.
    %
    %   Inputs:
    %       P0     - Stagnation pressure (Pa)
    %       P      - Static pressure (Pa)
    %       T0     - Stagnation temperature (K)
    %       k      - Heat capacity ratio (Cp/Cv)
    %       Rstar  - Specific gas constant (J/kg·K)
    %
    %   Output:
    %       w      - Speed of the gas (m/s)
    %
    % Compute the speed of the gas
    w = sqrt( (2 * k / (k - 1)) * Rstar * T0 * (1 - (P0 / P)^((1 - k) / k)) );
end
