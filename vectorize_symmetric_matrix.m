function [a, P] = vectorize_symmetric_matrix(A)
    % Ensure A is symmetric
%     if ~issymmetric(A)
%         error('Input matrix A must be symmetric.');
%     end
    
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
    
    % Create permutation matrix P
    P = zeros(n^2, n*(n+1)/2);
    
    count = 0;
    for col = 1:n
        for row = col:n
            count = count + 1;
            P((col-1)*n + row, count) = 1; % Adjusted indexing to match Vec(A)
            if row ~= col
                P((row-1)*n + col, count) = 1; % Symmetric pair
            end
        end
    end
end