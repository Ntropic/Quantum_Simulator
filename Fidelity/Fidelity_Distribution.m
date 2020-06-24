function [ min_F, F_av, F_bins, xs ] = Fidelity_Distribution( A_ex,A_ap,indexes,how_many,bins,bins_int )
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
%% Average Fidelity
s=length(A_ex);
B=A_ex'*A_ap;
     
F_av=1/(s*(s+1))*(trace(B*B')+abs(trace(B))^2);    %Fidelity of quantum operations

%% Compute minimum Fidelity
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

%% Compute how many Fidelities of random Wave vectors
%xs=linspace(1-bins_int(2),1-bins_int(1),bins);
bins_int=log10([bins_int(2),bins_int(1)]);
xs=10.^(linspace(bins_int(1),bins_int(2),bins));
fs=zeros(how_many,1);
for i=1:how_many
    phi=RandWave(s);
    fs(i)=abs(phi'*B*phi)^2;
end
F_bins=hist(1-fs,xs(end:-1:1));
F_bins=F_bins(:,end:-1:1);
end

