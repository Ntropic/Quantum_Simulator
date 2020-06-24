function [ U_MWG ] = U2MWG_Trotter( U,coupled_guides,l )
%H2MWG_Trotter Creates a multi Waveguide Unitarian from a single Waveguide
%Unitarian
%   - coupled_guides     -  is a M*2 list of M different interactions with each 2
%                           interaction partners (or alternatively a scalar
%                           defining the number of waveguides)
couples=coupled_guides;
%Determine number of waveguides
if length(couples)>1
    N=max(couples(:));
else 
    N=couples;
    %Now define nearest neighbor coupling
    couples=[(1:N-1)',(2:N)'];
end
n=size(U,1);
if size(U,1)~=size(U,2)
    error('Hamiltonians need to be symmetric')
end
s=log2(n);
if s~=round(s)
    error('Hamiltonian needs to be given in a qubit coding, with dimensions 2^n');
end


%Now create larger version of Hamiltonians
U_MWG=eye(2^(s/2*N));
if length(coupled_guides)>1
    for i=1:size(couples,1)
        %Determine indexes
        indexes=[s/2*(couples(i,1)-1)+(1:s/2),s/2*(couples(i,2)-1)+(1:s/2)];
        %Use Gate EMbedding algorithm to create Hamiltonian  
        U_MWG=Embed_Gate(U{i},indexes,s/2*N)*U_MWG;
    end
else
    for i=1:size(couples,1)
        %Determine indexes
        indexes=[s/2*(couples(i,1)-1)+(1:s/2),s/2*(couples(i,2)-1)+(1:s/2)];
        %Use Gate EMbedding algorithm to create Hamiltonian  
        U_MWG=Embed_Gate(U,indexes,s/2*N)*U_MWG;
    end
end
U_MWG=U_MWG^l;
end

