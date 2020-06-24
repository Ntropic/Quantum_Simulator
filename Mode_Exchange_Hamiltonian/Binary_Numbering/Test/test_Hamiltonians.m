%test_Hamiltonians.m
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

N=3;
n=2*N;
sizle=2^n;

%Load gates
load('../../../Gates_Table/elem_gates.mat','-mat')
load('../../../Gates_Table/comp_gates.mat','-mat')


%Create/Load Controlled Unitaries for different sizes
if test==1
    H=Exchange_Hamiltonian(N);

    mode={'fock',2^N-1,2};
    PColorMat(H,mode)

elseif test==2
    H=Exchange_Hamiltonian_Overfull(N);

    mode={'fock',2^N-1,2};
    PColorMat(H,mode)

elseif test==3
    H=Single_Particle_Exchange_Hamiltonian(N,[2,1,1]);

    mode={'fock',2^N-1,2};
    PColorMat(H,mode)

end