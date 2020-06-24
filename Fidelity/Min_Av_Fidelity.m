function [ min_F, F_av ] = Min_Av_Fidelity( A_ex,A_ap,indexes )
%MIN_AV_FIDELITY calculates the minimum and average fidelity 
%over all possible states between an approximation A_ap and an exact A_ex 
%-> indexing works for unitary matrices -> apart from whcih this method 
%   needs a more refined approach
if size(A_ex,1)~=size(A_ex,2)
    error('Matrices must be square (A_ex)')
end
if size(A_ap,1)~=size(A_ap,2)
    error('Matrices must be square (A_ap)')
end
if size(A_ex,1)~=size(A_ap,1)
    error('Matrices must be of equal size')
end

%Squeeze matrices to subspace with indexes
if nargin==3
    A_ex=A_ex(indexes,indexes);
    A_ap=A_ap(indexes,indexes);
end
s=(log2(length(A_ex)))/2;
B=A_ex'*A_ap;
     
F_av=1/(2^(2*s)*(2^(2*s)+1))*(trace(B*B')+abs(trace(B))^2);    %Fidelity of quantum operations

%Compute minimum Fidelity
lambda=eig(B);
lambda=mod(angle(lambda)+2*pi,2*pi);
lambda=sort(lambda);
d_lambda=diff([lambda;lambda(1)]);
m_d=max(d_lambda);
if m_d<pi
    min_F=0;
else
    min_F=cos(m_d/2)^2;
end

end

