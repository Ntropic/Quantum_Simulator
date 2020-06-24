%Fidelity_Compare_Photon_Numbers_Sym_Gauss_Jordan.m
clear all;
close all;
clc;

p=pwd;
if any(strfind(p,'\'));
    elem=strsplit(p,'\');
else
    elem=strsplit(p,'/');
end
shortened=fullfile(elem{1:end-3});
addpath(genpath(shortened));
%Load gates
load('../../../Gates_Table/elem_gates.mat','-mat')
load('../../../Gates_Table/comp_gates.mat','-mat')

%Photon number
n_m=1:7;

%Trotter Iteration number
steps=10;

N_m=2:4;

%Iteration time
t=pi/4;

name=['mwg_sym_compare_fid_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_s_' num2str(steps) '_t_' num2str(t) '_N_' num2str(max(N_m)) '_pauli_gray'];
F_av_ind=zeros(length(N_m),length(n_m),steps);

for h=1:length(N_m)
    N=N_m(h);
    fprintf([num2str(N) ' Waveguides\n'])
    
    connections=[(1:N-1)',(2:N)'];
    weights=ones(N-1,1);

    for i=1:length(n_m)
        fprintf(['  ' num2str(i) ' Photons\n'])
        n=n_m(i);

      	[H]=Gray_Exchange_Hamiltonian_Particles(n);
        %% Exact decomposition
        Hij=H2MWG(H,N);
        %FockPrint(Hij)
        U_exact=expm(-1i*Hij*t);

        H_parts=Gray_Hamiltonian_Steps(n);
        for l=1:steps
            fprintf(['    ' num2str(l) ' Trotter steps\n'])
            
            U_approx=U2MWG_Strang_Trotter(H_parts,t,l,N);
    
            indexes=Gray_Indexes(n,N);
            [F_av(h,i,l) F_av_ind(h,i,l)]=Average_Fidelity(U_approx,U_exact,indexes);
            
            %Plotting 
            o=Order_Of_Convergence(1-F_av_ind(h,i,1:l));
            plot(o)
            ylabel('Order of Covergence');
            xlabel('Trotter steps')
            title([num2str(N) 'WGs, ' num2str(n) 'Photons'])
            drawnow;
        end 
    end
end

save([name '.mat'],'F_av_ind','F_av');
