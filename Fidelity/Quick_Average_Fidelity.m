function [ F_av ] = Quick_Average_Fidelity( A_ex,A_ap )
%QUICK_AVERAGE_FIDELITY calculates the average fidelity over all possible states
%between an approximation A_ap and an exact A_ex but without subspaces via
%indexes
if size(A_ex,1)~=size(A_ex,2)
    error('Matrices must be square (A_ex)')
end
if size(A_ap,1)~=size(A_ap,2)
    error('Matrices must be square (A_ap)')
end
if size(A_ex,1)~=size(A_ap,1)
    error('Matrices must be of equal size')
end

s1=(log2(length(A_ex)))/2;
B1=A_ex'*A_ap;

F_av=1/(2^(2*s1)*(2^(2*s1)+1))*(trace(B1*B1')+abs(trace(B1))^2);    %Fidelity of quantum operations
end

