%test_Gray_Hamiltonians.m
%In binary encoding
clc;
clear all;
close all;

p=pwd;
if any(strfind(p,'\'));
    elem=strsplit(p,'\');
else
    elem=strsplit(p,'/');
end
shortened=fullfile(elem{1:end-3});
addpath(genpath(shortened));

%Number of qubits per mode
test=2;

N=5;
n=2*N;
sizle=2^n;

%Load gates
load('../../../Gates_Table/elem_gates.mat','-mat')
load('../../../Gates_Table/comp_gates.mat','-mat')


%Create/Load Controlled Unitaries for different sizes
if test==1
    H=Gray_Exchange_Hamiltonian(N);

    mode={'fock_gray',2^N-1,2};
    PColorMat(H,mode)

elseif test==2
    H=Gray_Exchange_Hamiltonian_Overfull(N);

    mode={'fock_gray',2^N-1,2};
    PColorMat(H,mode)

elseif test==3
    H=Single_Particle_Gray_Exchange_Hamiltonian(N,[3,2,-1]);

    mode={'fock_gray',2^N-1,2};
    PColorMat(H,mode)

end