%Higher_Order_Circ_XY_All_Interaction.m
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
load('../XY_quantum_gates.mat','-mat')

%% Change steps definition for different types of higher order interactions
n=5;
N=n-1;%Photon number

%Trotter step number
steps=5;

%Higher Trotter recursion order
rec_order=1;

%Iteration time
t=pi/4;

Hij=XY_Exchange_Hamiltonian(n);
U_exact=expm(-1i*t*Hij);

%Subspace for analysis of algorithm
index=2.^(0:n-1)+1;

Ni=N:-1:1;
in=1:N;
Jx=sqrt(in).*sqrt(Ni);

%Gate A
sizle=floor(n/2)*2;
anc_size=0;
nameA=['A_n_' num2str(n)];
circ_string=['A_{XY}^{n=' num2str(n) '}'];
gates_indexes={};
param_steps={};
gates_names={};
for i=1:2:2*floor(n/2)-1
    gates_indexes={gates_indexes{:},[i,i+1]};
    param_steps={param_steps{:},[{sym('phi'),sym('phi')*Jx(i)}]};
end
[gates_names{1:length(gates_indexes)}]=deal('XX_YY');
A=Create_Comp_Gate(nameA,sizle,anc_size,gates_names,gates_indexes,param_steps);
A=Generate_Gate_Circuit(A,circ_string,1:sizle,[],[]);

%Gate B
sizle=floor((n-1)/2)*2;
anc_size=0;
nameB=['B_n_' num2str(n)];
circ_string=['B_{XY}^{n=' num2str(n) '}'];
gates_indexes={};
param_steps={};
gates_names={};
for i=1:2:2*floor((n-1)/2)-1
    gates_indexes={gates_indexes{:},[i,i+1]};
    param_steps={param_steps{:},[{sym('phi'),sym('phi')*Jx(i+1)}]};
end
[gates_names{1:length(gates_indexes)}]=deal('XX_YY');
B=Create_Comp_Gate(nameB,sizle,anc_size,gates_names,gates_indexes,param_steps);
B=Generate_Gate_Circuit(B,circ_string,1:sizle,[],[]);   
    
for l=1:steps
    fprintf([num2str(l) ' Trotter steps\n'])
    timer=tic();

    %Create Trotter steps
    circ_string=['XY_{\text{approx, higher}}^{n=' num2str(n) '}'];
    circ_string_step=['XY_{\text{step, sym}}^{n=' num2str(n) '}'];
    gates_names={nameA,nameB};
    gates_indexes={[1:A.size],1+[1:B.size]};
    other_gates1=Comp_Gate_Merger(comp_gates,xy_gates,A,B);
    %This command creates a Trotter gate for all the steps -> Great for
    %creating a gate circuit
    pattern=ones(1,rec_order);
    [seq,d_seq,gate,step_gate]=Higher_Trotter(t,l,pattern,'',circ_string,circ_string_step,n,anc_size,gates_names,gates_indexes);
    %if l==1
    %    Gates2Circ2Preview(gate,elem_gates,other_gates1,'xy_gate','gate=expansion',1,30);
    %end
    
    %Do Trotter steps
    other_gates=Comp_Gate_Merger(other_gates1,gate,step_gate);
    other_gates=Gates_Tables_Prepare(other_gates,elem_gates);

    ind=Gate_Index_by_Name(gate.names,elem_gates,other_gates);
    U_approx=Expansion2Matrix(other_gates(ind(2)));
    %average_fidelity
    [F_av(l) F_av_ind(l)]=Average_Fidelity(U_approx,U_exact,index);
    time_l(l)=toc(timer);
end 

f1_cmp=1-F_av_ind;

figure()
%subplot(3,1,1:2);
plot(1:steps,f1_cmp,'r')
%xlabel('Trotter steps')
ylabel('1-Average Fidelity')
legend({'\mathrm{O}(t^2)'})
axis tight;
%matlab2tikz('Average_Fidelity_Convergence1_50.tex')
drawnow;


figure()
n=2:steps;
logn=log((n-1)./n);
u1_diff=log(f1_cmp(2:end))-log(f1_cmp(1:end-1));
semilogx(n,u1_diff./logn,'r')
xlabel('Trotter steps')
ylabel('Convergence Order')
legend({'\mathrm{O}(t^2)','\mathrm{O}(t^3)'})
axis tight;
%matlab2tikz('Average_Fidelity_Convergence_order1_50.tex')
%drawnow;
%xlabel('Trotter steps')
%ylabel('Convergence Order')
%legend({'\mathrm{O}(t^2)','\mathrm{O}(t^3)'})
%axis tight;
%matlab2tikz('Average_Fidelity_Convergence_order1_50.tex')
%drawnow;