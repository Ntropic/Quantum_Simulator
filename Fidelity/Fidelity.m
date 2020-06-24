function [ F ] = Fidelity( A_ex,A_ap,indexes,phi )
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
if length(indexes)>0
    A_ex=A_ex(indexes,indexes);
    A_ap=A_ap(indexes,indexes);
end

B=A_ex'*A_ap;

%% Compute how many Fidelities of random Wave vectors
%xs=linspace(1-bins_int(2),1-bins_int(1),bins);
how_many=size(phi,2);
for i=1:how_many
    phi1=phi(:,i);
    F(i)=abs(phi1'*B*phi1)^2;
end

end

