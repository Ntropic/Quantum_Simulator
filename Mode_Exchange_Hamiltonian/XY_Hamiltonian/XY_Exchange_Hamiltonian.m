function [ H ] = XY_Exchange_Hamiltonian(n)
%XY_Exchange_Hamiltonian Create n qubit exchange matrix for N modes with
%exchange between modes ij and ij+1

N=n-1;  %Number of photons
index=2.^(0:n-1)+1;

Ni=N:-1:1;
in=1:N;
Jx=sqrt(in).*sqrt(Ni);

H=sparse(2^n,2^n);
for i=1:length(Jx)
    H(index(i),index(i+1))=Jx(i);
    H(index(i+1),index(i))=Jx(i);
end

end