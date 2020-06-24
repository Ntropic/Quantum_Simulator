function [ H_MWG ] = H2MWG( H,coupled_guides,coupling_strengths )
%H2MWG Creates a multi Waveguide Hamiltonian from a single Waveguide
%Hamiltonian
%   - coupled_guides     -  is a M*2 list of M different interactions with each 2
%                           interaction partners (or alternatively a scalar
%                           defining the number of waveguides)
%   - coupling_strengths -  the strengths of couplings between the
%                           waveguides

%Determine number of waveguides
if length(coupled_guides)>1
    N=max(coupled_guides(:));
else 
    N=coupled_guides;
    %Now define nearest neighbor coupling
    coupled_guides=[(1:N-1)',(2:N)'];
end
n=size(H,1);
if size(H,1)~=size(H,2)
    error('Hamiltonians need to be symmetric')
end
s=log2(n);
if s~=round(s)
    error('Hamiltonian needs to be given in a qubit coding, with dimensions 2^n');
end

if nargin==2
    %Coupling strengths weren't defined -> choose equivalent strength of 1
    coupling_strengths=ones(size(coupled_guides,1),1);
end

%Now create larger version of Hamiltonians
H_MWG=zeros(2^(s/2*N));
for i=1:length(coupling_strengths)
    %Determine indexes
    indexes=[s/2*(coupled_guides(i,1)-1)+(1:s/2),s/2*(coupled_guides(i,2)-1)+(1:s/2)];
    %Use Gate EMbedding algorithm to create Hamiltonian  
    H_MWG=H_MWG+coupling_strengths(i)*Embed_Gate(H,indexes,s/2*N);
end
end

