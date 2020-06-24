function [ K,p,ind,chan ] = Pauli_Twirling( n,t_step,T1,T2 )
%PAULI_TWIRLING creates a list of Pauli Twirling Channels for n qubits, 
%time steps of length t_step, relaxation time T1 and decoherence time T2
%depth specifies how many pauli channels are maximally present within one 
%Twirl (default n)
%
%Model details: rho0=[1-rho11,rho01
%                     rho01* ,rho11]
%   t=t_step    ---> [1-rho11*e^(-t/T1),rho01*e^(-t/T2)
%                     rho01'*e^(-t/T2) ,rho11*e^(-t/T1)]
%       -> Dephasing noise propto 1/f
   
if nargin==1
    t_step=sym('t_step');
    T1=sym('T1');
    T2=sym('T2');
end

%Probabilities for one qubit
px=1/4*(1-exp(-t_step/T1));
py=px;
pz=1/2-px-1/2*exp(-t_step/T2);
pI=1-px-py-pz;

%Pauli matrix
x=[0 1; 1 0];
y=[0 -1i; 1i 0];
z=[1 0; 0 -1];
I=[1 0; 0 1];

sig={x,y,z};
pr=[px,py,pz];

K={};
ind={};
chan={};
%Use Kronecker product (kron) | Tensor product
p(1)=pI^n;
for j=1:n
    p0=pI^(n-j);
    how_many=factorial(n)/factorial(n-j)/factorial(j);
    index_list=Pick_J_of_N(j,n); %Ordered list
    comb=Comb_J_of_N(j,3);
    for k=1:how_many
        i=index_list(k,:);
        for m=1:size(comb,1)
            mat=1;
            p=[p;prod(pr(comb(m,:)))*p0];
            for l=1:length(i)
                if l==1
                    mat=kron(diag(ones(2^(i(l)-1),1)),mat);
                else
                    mat=kron(diag(ones(2^(i(l)-i(l-1)-1),1)),mat);
                end
                mat=kron(sig{comb(m,l)},mat);
            end
            mat=kron(diag(ones(2^(n-i(l)),1)),mat);
            K={K{:},mat};
            ind={ind{:},i};
            chan={chan{:},comb(m,:)};
        end
    end
end

end

