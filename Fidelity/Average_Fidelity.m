function [ F_av1,F_av2 ] = Average_Fidelity( A_ex,A_ap,indexes )
%AVERAGE_FIDELITY calculates the average fidelity over all possible states
%between an approximation A_ap and an exact A_ex -> indexing works for
%unitary matrices -> apart from whcih this method needs a more refined
%approach
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

F_av1=1/(2^(2*s1)*(2^(2*s1)+1))*(trace(B1*B1')+abs(trace(B1))^2);    %Fidelity of quantum operations


%Squeeze matrices to subspace with indexes
if nargin==3
    A_ex=A_ex(indexes,indexes);
    A_ap=A_ap(indexes,indexes);
end
s2=(log2(length(A_ex)))/2;
B2=A_ex'*A_ap;
     
F_av2=1/(2^(2*s2)*(2^(2*s2)+1))*(trace(B2*B2')+abs(trace(B2))^2);    %Fidelity of quantum operations

end

