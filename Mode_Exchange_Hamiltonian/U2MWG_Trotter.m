function [ U ] = U2MWG_Trotter ( H,t,l,coupled_guides,coupling_strengths )
%H2MWG Creates a multi Waveguide Hamiltonian from a 2 Waveguide
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

if isa(H,'cell')==0
    n=size(H,1);
    if size(H,1)~=size(H,2)
        error('Hamiltonians need to be symmetric')
    end
else
    n=size(H{1},1);
    if size(H{1},1)~=size(H{1},2)
        error('Hamiltonians need to be symmetric')
    end
end
s=log2(n);
if s~=round(s)
    error('Hamiltonian needs to be given in a qubit coding, with dimensions 2^n');
end

if nargin==4
    coupling_strengths=ones(size(coupled_guides,1),1);
end

%Now create larger version of Hamiltonians
if nargin==4
    U=eye(2^(s/2*N));
    if isa(H,'cell')
        Uij=eye(size(H{1}));
        for i=1:length(H)
            Uij=expm(-1i*H{i}*t/l)*Uij;
        end
    else
        Uij=expm(-1i*H*t/l);
    end
    for i=1:length(coupling_strengths)
        %Determine indexes
        indexes=[s/2*(coupled_guides(i,1)-1)+(1:s/2),s/2*(coupled_guides(i,2)-1)+(1:s/2)];
        %Use Gate EMbedding algorithm to create Hamiltonian  
        U=Embed_Gate(Uij,indexes,s/2*N)*U;
    end
else
    U=eye(2^(s/2*N));
    for i=1:length(coupling_strengths)
        if isa(H,'cell')
            Uij{i}=eye(size(H));
            for j=1:length(H)
                Uij{i}=expm(-1i*H{j}*t/l)*Uij;
            end
        else
            Uij{i}=expm(-1i*coupling_strengths(i)*H*t/l);
        end
        
    end
    for i=1:length(coupling_strengths)
        %Determine indexes
        indexes=[s/2*(coupled_guides(i,1)-1)+(1:s/2),s/2*(coupled_guides(i,2)-1)+(1:s/2)];
        %Use Gate EMbedding algorithm to create Hamiltonian  
        U=Embed_Gate(Uij{i},indexes,s/2*N)*U;
    end
end
U=U^l;
end

