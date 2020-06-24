function [ U ] = U2MWG_Higher_Trotter( H,t,l,pattern,coupled_guides,coupling_strengths )
%H2MWG_Strang_Trotter Creates a multi Waveguide Unitarian from a single Waveguide
%Unitarian
%   - coupled_guides     -  is a M*2 list of M different interactions with each 2
%                           interaction partners (or alternatively a scalar
%                           defining the number of waveguides)
%   - U={U_half_step,U_full_step}...
[seq,d_seq]=Higher_Trotter_Time_Steps(pattern);

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

if nargin==5
    coupling_strengths=ones(size(coupled_guides,1),1);
end

%Now create larger version of Hamiltonians
if nargin==5
    U=eye(2^(s/2*N));
    for j=1:length(d_seq)
        if isa(H,'cell')
         	Uij{1}=eye(size(H{1}));
            Uij{2}=eye(size(H{1}));
            for i=1:length(H)
                Uij{1}=expm(-1i*H{i}*d_seq(j)*t/l/2)*Uij{1};
            end
            for i=length(H):-1:1
                Uij{2}=expm(-1i*H{i}*d_seq(j)*t/l/2)*Uij{2};
            end

            for i=1:size(coupled_guides,1)
                %Determine indexes
                indexes=[s/2*(coupled_guides(i,1)-1)+(1:s/2),s/2*(coupled_guides(i,2)-1)+(1:s/2)];
                %Use Gate EMbedding algorithm to create Hamiltonian  
                U=Embed_Gate(Uij{1},indexes,s/2*N)*U;
            end
            for i=size(coupled_guides,1):-1:1
                %Determine indexes
                indexes=[s/2*(coupled_guides(i,1)-1)+(1:s/2),s/2*(coupled_guides(i,2)-1)+(1:s/2)];
                %Use Gate EMbedding algorithm to create Hamiltonian  
                U=Embed_Gate(Uij{2},indexes,s/2*N)*U;
            end
        else
            Uij={expm(-1i*H*d_seq(j)*t/l/2),expm(-1i*H*d_seq(j)*t/l)};

            for i=1:size(coupled_guides,1)-1
                %Determine indexes
                indexes=[s/2*(coupled_guides(i,1)-1)+(1:s/2),s/2*(coupled_guides(i,2)-1)+(1:s/2)];
                %Use Gate EMbedding algorithm to create Hamiltonian  
                U=Embed_Gate(Uij{1},indexes,s/2*N)*U;
            end
            i=size(coupled_guides,1);
            indexes=[s/2*(coupled_guides(i,1)-1)+(1:s/2),s/2*(coupled_guides(i,2)-1)+(1:s/2)];
            U=Embed_Gate(Uij{2},indexes,s/2*N)*U;
            for i=size(coupled_guides,1)-1:-1:1
                %Determine indexes
                indexes=[s/2*(coupled_guides(i,1)-1)+(1:s/2),s/2*(coupled_guides(i,2)-1)+(1:s/2)];
                %Use Gate EMbedding algorithm to create Hamiltonian  
                U=Embed_Gate(Uij{1},indexes,s/2*N)*U;
            end
        end
    end
else
    U=eye(2^(s/2*N));
    for j=1:length(d_seq)
        if isa(H,'cell')
            for i=1:size(coupled_guides,1)-1
                Uij{i,1}=eye(size(H{1}));
                Uij{i,2}=eye(size(H{1}));

                for i=1:length(H)
                    Uij{i,1}=expm(-1i*d_seq(j)*H{i}*t/l/2)*Uij{1};
                end
                for i=length(H):-1:1
                    Uij{i,2}=expm(-1i*d_seq(j)*H{i}*t/l/2)*Uij{2};
                end
            end
            for i=1:size(coupled_guides,1)
                %Determine indexes
                indexes=[s/2*(coupled_guides(i,1)-1)+(1:s/2),s/2*(coupled_guides(i,2)-1)+(1:s/2)];
                %Use Gate EMbedding algorithm to create Hamiltonian  
                U=Embed_Gate(Uij{i,1},indexes,s/2*N)*U;
            end
            for i=size(coupled_guides,1):-1:1
                %Determine indexes
                indexes=[s/2*(coupled_guides(i,1)-1)+(1:s/2),s/2*(coupled_guides(i,2)-1)+(1:s/2)];
                %Use Gate EMbedding algorithm to create Hamiltonian  
                U=Embed_Gate(Uij{i,2},indexes,s/2*N)*U;
            end
        else
            for i=1:length(coupled_guides)
                Uij{i,1}=expm(-1i*coupling_strengths(i)*d_seq(j)*H*t/l/2);
                Uij{i,2}=expm(-1i*coupling_strengths(i)*d_seq(j)*H*t/l);
            end

            for i=1:size(coupled_guides,1)-1
                %Determine indexes
                indexes=[s/2*(coupled_guides(i,1)-1)+(1:s/2),s/2*(coupled_guides(i,2)-1)+(1:s/2)];
                %Use Gate EMbedding algorithm to create Hamiltonian  
                U=Embed_Gate(Uij{i,1},indexes,s/2*N)*U;
            end
            i=size(coupled_guides,1);
            indexes=[s/2*(coupled_guides(i,1)-1)+(1:s/2),s/2*(coupled_guides(i,2)-1)+(1:s/2)];
            U=Embed_Gate(Uij{i,2},indexes,s/2*N)*U;
            for i=size(coupled_guides,1)-1:-1:1
                %Determine indexes
                indexes=[s/2*(coupled_guides(i,1)-1)+(1:s/2),s/2*(coupled_guides(i,2)-1)+(1:s/2)];
                %Use Gate EMbedding algorithm to create Hamiltonian  
                U=Embed_Gate(Uij{i,1},indexes,s/2*N)*U;
            end
        end
    end
end
U=U^l;
end

