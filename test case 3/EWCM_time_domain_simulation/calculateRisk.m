function Risk = calculateRisk(X0, x, M_p, margin_percent)
    % Function to calculate risk based on transient simulation results
    % Inputs:
    %   X0 - Initial condition vector
    %   x - Simulation results matrix
    %   M_p - Number of pipes (used for indexing)
    %   margin_percent - Allowed margin percentage (e.g., 0.01 for 1%)
    % Output:
    %   Risk - Computed risk value
    
    % Extract the relevant portion of the initial conditions
    X0_risk = X0(M_p+1:end);
    x_risk = x(:, M_p+1:end);
    
    % Ensure row vector format
    X0_risk = X0_risk(:)';
    
    % Compute upper and lower bounds
    upper_bound = X0_risk .* (1 + margin_percent);
    lower_bound = X0_risk .* (1 - margin_percent);
    
    % Replicate bounds for matrix operations
    upper_bound_mat = repmat(upper_bound, size(x_risk, 1), 1);
    lower_bound_mat = repmat(lower_bound, size(x_risk, 1), 1);
    
    % Initialize deviation matrix
    x2 = zeros(size(x_risk));
    
    % Identify values exceeding the bounds
    above_threshold = x_risk > upper_bound_mat;
    below_threshold = x_risk < lower_bound_mat;
    
    % Compute deviations from the allowed range
    x2(above_threshold) = x_risk(above_threshold) - upper_bound_mat(above_threshold);
    x2(below_threshold) = lower_bound_mat(below_threshold) - x_risk(below_threshold);
    
    % Compute the total risk
    Risk = sum(sum(x2));
end