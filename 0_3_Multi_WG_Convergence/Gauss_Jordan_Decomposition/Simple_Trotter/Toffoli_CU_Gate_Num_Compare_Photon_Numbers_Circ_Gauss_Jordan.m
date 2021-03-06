%Toffoli_CU_Gate_Num_Compare_Photon_Numbers_Circ_Gauss_Jordan.m
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

%Iteration time
t=pi/4;

name=['twowg_trot_compare_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_s_' num2str(steps) '_t_' num2str(t) '_gj'];
F_av_ind=zeros(length(n_m),steps);
gate_num=zeros(length(n_m),steps);
s_q1=zeros(length(n_m),steps);
s_q2=zeros(length(n_m),steps);
time_l=zeros(length(n_m),steps);

    N=2;
    fprintf([num2str(N) ' Waveguides\n'])
    
    connections=[(1:N-1)',(2:N)'];
    weights=ones(N-1,1);

    for i=1:length(n_m)
        fprintf(['  ' num2str(i) ' Photons\n'])
        n=n_m(i);

      	[H]=Gray_Exchange_Hamiltonian_Particles(n);
        %% Exact decomposition
        %Hij=H2MWG(H,N);
        %FockPrint(Hij)
        %U_exact=expm(-1i*Hij*t);

        for l=1:steps
            fprintf(['    ' num2str(l) ' Trotter steps\n'])

            U_2=expm(-1i*H*t/l);
            gj_dec_gate=Fast_Gauss_Jordan_Decomposition(U_2);
            other_gates2=Comp_Gate_Merger(comp_gates,gj_dec_gate);
            
            [toffolis,cus]=Some_Gate_Counter(gj_dec_gate);

            s_q1(i,l)=toffolis;
            s_q2(i,l)=cus;
            drawnow;
        end 
    end

save([name '_toff_cu.mat'],'s_q1','s_q2');

%% Create a plot of the number of gates required to perform an exact decomposition of a single Trotter step for 2 interacting waveguides depending on the time evolution per step
t2=t;
for l=1:steps
    t(l)=t2/l;
end
toffoli_gates=s_q1(:,:);
cu_gates=s_q2(:,:);
total_gates=toffoli_gates+cu_gates;
surf(t,n_m,total_gates)
axis tight;
ylabel('N_\text{max}');
xlabel('t')

figure()
plot(n_m,total_gates(:,1),'b')
hold on;
plot(n_m,toffoli_gates(:,1),'r')
axis tight;


