function [H]= RandH(n,s)
%Generate random Hamiltonian matrices
if nargin==1
    X = complex(rand(n),rand(n))/sqrt(2);
    H=(X+X')/2;
else
    if s==1 %Create sparse matrix
        X=complex(rand(n),rand(n))/sqrt(2);
        X(abs(X)<0.8)=0;
        X=round(X);
        H=(X+X')/2;
    else
        X = complex(rand(n),rand(n))/sqrt(2);
        H=(X+X')/2;
    end
end
end