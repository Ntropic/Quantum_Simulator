%Circ_Pauli_Gray_All_Interaction.m
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
load('cxy_gates.mat','-mat')
other=Comp_Gate_Merger(comp_gates,cxy_gates);

%Photon number
n=4;

%Trotter Iteration number
steps=2;

%Iteration time
t=pi/4;
circing=0;

%% Exact decomposition
[Hij,s]=Gray_Exchange_Hamiltonian_Particles(n); %s is qubits per mode 2s is total qubit number
%FockPrint(Hij)
U_exact=expm(-1i*Hij*t);

    [init main connector outit gate_names]=Split_Trot_Pauli_Gray_Gates(n);

    inds=main.size+main.anc_size;

    other_gates=Comp_Gate_Merger(other,init,main,connector,outit);
    other_gates=Gates_Tables_Prepare(other_gates,elem_gates);   

    sizle=init.size;
    anc_size=init.anc_size;
    for l=1:steps
        fprintf(['  ' num2str(l) ' Trotter steps\n'])
        timer=tic();

        trott_gate=Trotter_InOut(t,l,'','',sizle,anc_size,gate_names);

        if l==1
            if circing==1
                Gates2Circ2Preview(trott_gate,elem_gates,other_gates,'test_trott_gate','gate=expansion',2,150);
            end
        end

        other_gates2=Comp_Gate_Merger(other_gates,trott_gate);
        other_gates2=Gates_Tables_Prepare(other_gates2,elem_gates);   

        ind=Gate_Index_by_Name(trott_gate.names,elem_gates,other_gates2);
        U_approx=Expansion2Matrix(other_gates2(ind(2)));
        %average_fidelity
        indexes=Gray_Indexes(n);
        [F_av(l) F_av_ind(l)]=Average_Fidelity(U_approx,U_exact,indexes);
    end 
    
f1_cmp=1-F_av_ind;

figure()
plot(1:steps,f1_cmp(1,:),'r')
%xlabel('Trotter steps')
ylabel('1-Average Fidelity')
axis tight;
%matlab2tikz('Average_Fidelity_Convergence1_50.tex')
drawnow;


%n=2:steps;
%logn=log((n-1)./n);
%u1_diff=log(f1_cmp(2:end))-log(f1_cmp(1:end-1));
%semilogx(n,u1_diff./logn','r')
%title('N=10')
%xlabel('Trotter steps')
%ylabel('Convergence Order')
%legend({'\mathrm{O}(t^2)','\mathrm{O}(t^3)'})
%axis tight;
%matlab2tikz('Average_Fidelity_Convergence_order1_50.tex')
%drawnow;