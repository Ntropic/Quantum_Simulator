function [ K,comb ] = Kraus_Twirling( n,t_step,T1,T2 )
%KRAUS_TWIRLING creates a list of Kraus Twirling Channels for n qubits, 
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
g=1-exp(-t_step/T1);
l=exp(-t_step/T1)-exp(-2*t_step/T2);
p1=exp(-2*t_step/T2);

if l<0
    error('Decoherence time must be smaller than 2 times the decay time: T2<2*T1')
end

%Pauli matrix
E1=[1 0;0 sqrt(p1)];
E2=[0 sqrt(g); 0 0];
E3=[0 0; 0 sqrt(l)];

sig={E1,E2,E3};

K={};
ind={};
chan={};
%Use Kronecker product (kron) | Tensor product

comb=Comb_J_of_N(n,3);
for m=1:size(comb,1)
    mat=1;
    for l=1:size(comb,2)
        mat=kron(sig{comb(m,l)},mat);
    end
    K={K{:},mat};
end

end

