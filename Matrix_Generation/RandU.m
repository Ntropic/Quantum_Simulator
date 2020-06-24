function [U,verify]= RandU(n)
%Generate random unitary matrices
    X = complex(rand(n),rand(n))/sqrt(2);
    [Q,R] = qr(X);
    R = diag(diag(R)./abs(diag(R)));
    U = Q*R;
    verify = diag(ones(n,1))-U*U';
end