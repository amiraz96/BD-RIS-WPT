function z1 = gaussian_randomization_cholesky(X, num_trials, K1)
    % Input: X - the higher rank matrix (should be positive semidefinite)
    %        num_trials - the number of randomization trials to run
    % Output: X1 - the best rank-one approximation matrix
    elem_num = size(X, 1);
    % Step 1: Cholesky factorization
    L = chol(X, 'lower');  % Perform Cholesky decomposition: X = L * L^T

    % Step 2: Initialize the best rank-one solution
    best_obj = inf;
    X1 = zeros(size(X));
    z1 = [];
    % Step 3: Perform Gaussian randomization
    for i = 1:num_trials
        % Generate a random Gaussian vector
        v = randn(elem_num, 1) + 1j*randn(elem_num, 1);  % Random complex vector
        v = v / norm(v);  % Normalize u_i to be on the unit sphere
        
        % Step 4: Compute z = L * v
        z = L * v;

        % Step 5: Construct the rank-one matrix
        X_trial = z * z';

        % Step 6: Evaluate the objective function (for example, maximize trace(X'X))
        obj_value = trace(X_trial * K1);
        
        if i == 1
            [V, F] = eigs(X);
            eigenvalues = diag(F);
            [~, idx] = max(eigenvalues);
            dominant_eigenvector = V(:, idx);
            dominant_eigenvalue = eigenvalues(idx);
            z = sqrt(dominant_eigenvalue) * dominant_eigenvector;
            X_trial = z * z';
            obj_value = trace(X_trial * K1);
        end

        % Step 7: Update the best solution if this trial is better
        if obj_value < best_obj
            best_obj = obj_value;
            X1 = X_trial;
            z1 = z;
        end
    end
end