A = sparse(1000000, 1000000); % Create sparse matrix
% Example: Add non-zero elements
A(1, 1) = 5; 
A(2, 2) = 10;

B = rand(1000000, 1); % Dense right-hand side
X = A * B; % Solve the linear system
