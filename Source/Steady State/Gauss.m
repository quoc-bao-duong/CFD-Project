function X=Gauss(A,B)
%{
    AX=B
    Input : A,B
    Output: X

%}
    n=length(A); % A, Matrix Size : n x n

    % Step 1: Forward Elimination
    % k
    % i: row 
    % j: column 
    
    for k = 1 : n-1
        % Elimination A(i,k)
        for i = k + 1 : n
            factor = A(i,k)/A(k,k);
            for j = k + 1 : n
                A(i,j) = A(i,j) - factor*A(k,j);
            end
            B(i) = B(i) - factor*B(k);
        end
    end
    
    % Step 2: Back Substitution

    X(n) = B(n)/A(n,n);
    for i = n-1:-1:1
        sum = B(i);
        for j = i + 1 : n
            sum = sum - A(i,j)*X(j);
        end
        X(i) = sum/A(i,i);
    end
X=X';
end