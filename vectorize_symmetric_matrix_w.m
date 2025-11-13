function [a] = vectorize_symmetric_matrix_w(A)
    % Ensure A is symmetric
    if ~issymmetric(A)
        error('Input matrix A must be symmetric.');
    end
    
    % Get the size of A
    n = size(A, 1);
    
    % Extract unique elements from the lower triangular part of A
    a = zeros(n*(n+1)/2, 1);
    idx = 1;
    for j = 1:n
        for i = j:n
            a(idx) = A(i, j);
            idx = idx + 1;
        end
    end
end