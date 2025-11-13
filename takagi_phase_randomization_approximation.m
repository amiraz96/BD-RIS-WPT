function best_U = takagi_phase_randomization_approximation(A, num_samples, H_R, H_I, w0, K2, K4)
    % A: The complex symmetric matrix to approximate
    % num_samples: Number of random samples to generate

    % Ensure A is symmetric
%     A = (A + A') / 2;
    
    % Size of matrix A
    matrix_size = size(A, 1);
    
    % Initialize variables
    best_U = [];
    min_diff = 0;
    
    % Takagi Decomposition: A = U * Sigma * U'
    [U, ~] = svd(A);  % SVD for complex symmetric matrix
    [~, Sigma] = eig(A);
    % Diagonal elements of Sigma
    singular_values = diag(Sigma);
    N_f = size(H_R, 2);
    
    for i = 1:num_samples
        % Generate random phases for the singular values
        random_phases = exp(1i * 2 * pi * rand(matrix_size, 1));
        
        % Construct the diagonal matrix with random phases
        if i == 1
            Lambda = diag(singular_values./abs(singular_values));
        else
            Lambda = diag(random_phases);
        end
        
        % Reconstruct the symmetric unitary matrix
        U_new = U * Lambda * U.';
        
%         equiv_chan = zeros(N_f, 1);
%         for n = 1:N_f
%             equiv_chan(n) = H_R(:, n).'*U_new*H_I(:, n);
%         end
        equiv_chan = sum(H_R .* (U_new * H_I), 1).';
        [~, f0value] = aprox_fixed_w_new(w0, equiv_chan, K2, K4, N_f, 1, 1);
        diff = real(f0value);
        
        % Update the best approximation if the current one is better
        if diff > min_diff
            min_diff = diff;
            best_U = U_new;
        end
    end
end