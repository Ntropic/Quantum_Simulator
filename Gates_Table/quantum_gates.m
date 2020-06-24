%quantum_gates.m
clc;
clear all;
close all;

addpath(genpath('..\'))

%The basic one qubit gates
H=1/sym(sqrt(2))*[1 1; 1 -1];        %Hadamard gate
X=[0 1; 1 0];                   %Pauli-X gate
Y=[0 -1i; 1i 0];                %Pauli-Y gate
Z=[1 0; 0 -1];                  %Pauli-Z gate

phi=sym('phi');
P01=sym(exp(1i*phi)*[1 0; 0 1]);  %Global Phase shift gate -> no cost, just for prettier results -> only gate with no error!
P1=sym([1 0; 0 exp(1i*phi)]);  
P0=sym([exp(1i*phi) 0; 0 1]); 
EX=sym(expm(-1i*phi/2*[0 1; 1 0]));%*exp(-1i*(2*pi-phi/2));
EY=sym(expm(-1i*phi/2*[0 -1i; 1i 0]));%*exp(-1i*(2*pi-phi/2));
EZ=sym(expm(-1i*phi/2*[1 0; 0 -1]));%*exp(-1i*(2*pi-phi/2));

%The basic two qubit gates
%iSWAP=sym([1 0 0 0; 0 cos(phi) -1i*sin(phi) 0; 0 -1i*sin(phi) cos(phi) 0;
%0 0 0 1]); ->Multi qubit implementation required
%SQRTiSWAP=subs(iSWAP,phi,pi/2);
CZ=sym([1 0 0 0; 0 1 0 0 ; 0 0 1 0; 0 0 0 -1]);  %Controlled Pauli-Z gate
CZ_phi=sym([1 0 0 0; 0 1 0 0 ; 0 0 1 0; 0 0 0 exp(1i*phi)]);  %Controlled Pauli-Z gate

%% Create Structure with elementary qbit gates
elem_gate_names={'H',{'X','NOT'},'Y','Z',{'P01','Phi'},{'P1','P'},'P0','EX','EY','EZ','CZ',{'CZ_phi','CPHASE1','C1PHASE1','CPHASE'}};
elem_gate_circs={'\gate{H}','\gate{X}','\gate{Y}','\gate{Z}','','\gate{e^{-i\frac{\phi}{2}z}}','\gate{e^{i\frac{\phi}{2}z}}',...
    '\gate{e^{-i\frac{\phi}{2}x}}',...
    '\gate{e^{-i\frac{\phi}{2}y}}',...
    '\gate{e^{-i\frac{\phi}{2}z}}',...
    {[1,2;1,1],'\gate{Z}','\ctrl{_1_}'},{[1,2;1,1],'\gate{e^{-i\frac{\phi}{2}z}}','\ctrl{_1_}'}};
    %{[1,2;1,1],'\qswap','qswap'},...
    %{[1,2;1,1],'\targ','\ctrl{_1_}'},...
    
elem_gate_size=[1,1,1,1,1,1,1,1,1,1,2,2];
elem_gate_pta={'','','','','0','','','','','','cz','cz'}; %Error model ->think about changing error parameters for cz_phi

elem_gates=struct('names',elem_gate_names,'circuit',elem_gate_circs,'twirling',elem_gate_pta);
for i=1:length(elem_gate_names);
    elem_gates(i).size=elem_gate_size(i);
    if iscell(elem_gates(i).names)
        elem_gates(i).matrix=eval(elem_gates(i).names{1});
    else
        elem_gates(i).matrix=eval(elem_gates(i).names);
    end
    %Create matlabFunction's
    if isa(elem_gates(i).matrix,'sym')
        a=symvar(elem_gates(i).matrix);
        if length(a)>0
            char_a={};
            for j=1:length(a)
                char_a={char_a{:},char(a(j))};
            end
            elem_gates(i).fun_mat=matlabFunction(elem_gates(i).matrix,'Vars',a);
            elem_gates(i).fun_vars=char_a;
        else %Create double 
            elem_gates(i).fun_mat=double(elem_gates(i).matrix);
            elem_gates(i).fun_vars=[];
        end
    else
        elem_gates(i).fun_mat=double(elem_gates(i).matrix);
        elem_gates(i).fun_vars=[];
    end
end

save('elem_gates.mat','elem_gates');

%% Create structure with composite qubits
n=1;
comp_gates_names={'ZZ','ZZ_min','XX','XX_min','YY','YY_min',{'CX','CNOT1','CNOT','CNOT1_12'},'SWAP',{'CNOT0','CNOT0_12'},'C0PHASE1','C1PHASE0','C0PHASE0','TOFFOLI11','TOFFOLI01','TOFFOLI00',{'TOPHILI','CCZ_Phi'},'U',{'CU','C1U'},'C0U'};
comp_gates=struct('names',comp_gates_names);
[comp_gates(:).matrix]=deal([]);

%% ZZ
comp_gates(n).size=2;
comp_gates(n).anc_size=0;
comp_gates(n).step_num=6;
comp_gates(n).steps.index={1,[1,2],1,2,[1,2],2};
comp_gates(n).steps.gates={'X','CZ_phi','X','X','CZ_phi','X'};
comp_gates(n).steps.param={[],[],[],[],[],[]};
comp_gates(n).circuit={[1,2;1,1],'\multigate{1}{e^{-i\frac{\phi}{2}z z}}','\ghost{e^{-i\frac{\phi}{2}z z}}';...
        [1,2;1,1],'\gate{e^{-i\frac{\phi}{2}z z}}','\sgate{e^{-i\frac{\phi}{2}z z}}{_1_}'};
n=n+1;
%% ZZ_min
comp_gates(n).size=2;
comp_gates(n).anc_size=0;
comp_gates(n).step_num=6;
comp_gates(n).steps.index={1,[1,2],1,2,[1,2],2};
comp_gates(n).steps.gates={'X','CZ_phi','Y','Y','CZ_phi','X'};
comp_gates(n).steps.param={[],[],[],[],[],[]};
comp_gates(n).circuit={[1,2;1,1],'\multigate{1}{e^{i\frac{\phi}{2}z z}}','\ghost{e^{i\frac{\phi}{2}z z}}';...
    [1,2;1,1],'\gate{e^{i\frac{\phi}{2}z z}}','\sgate{e^{i\frac{\phi}{2}z z}}{_1_}'};
n=n+1;
%% XX
comp_gates(n).size=2;
comp_gates(n).anc_size=0;
comp_gates(n).step_num=5;
comp_gates(n).steps.index={1,2,[1,2],1,2};
comp_gates(n).steps.gates={'EY','EY','ZZ','EY','EY'};
comp_gates(n).global_phase=exp(-1i*sym('phi')/2);
comp_gates(n).steps.param={[{sym('phi'),sym(pi/2)}],[{sym('phi'),sym(pi/2)}],[{sym('phi'),sym('phi')}],[{sym('phi'),sym(-pi/2)}],[{sym('phi'),sym(-pi/2)}]};
comp_gates(n).circuit_subs={[{'-i\frac{\phi}{2}','-i\frac{\pi}{4}'}],[{'-i\frac{\phi}{2}','-i\frac{\pi}{4}'}],[],[{'-i\frac{\phi}{2}','i\frac{\pi}{4}'}],[{'-i\frac{\phi}{2}','i\frac{\pi}{4}'}]};
comp_gates(n).circuit={[1,2;1,1],'\multigate{1}{e^{-i\frac{\phi}{2}x x}}','\ghost{e^{-i\frac{\phi}{2}x x}}';...
    [1,2;1,1],'\gate{e^{-i\frac{\phi}{2}x x}}','\sgate{e^{-i\frac{\phi}{2}x x}}{_1_}'};
n=n+1;
%% XX_min
comp_gates(n).size=2;
comp_gates(n).anc_size=0;
comp_gates(n).step_num=5;
comp_gates(n).steps.index={1,2,[1,2],1,2};
comp_gates(n).steps.gates={'EY','EY','ZZ_min','EY','EY'};
comp_gates(n).global_phase=exp(-1i*sym('phi')/2);
comp_gates(n).steps.param={[{sym('phi'),pi/2}],[{sym('phi'),pi/2}],[],[{sym('phi'),-pi/2}],[{sym('phi'),-pi/2}]};
comp_gates(n).circuit_subs={[{'-i\frac{\phi}{2}','-i\frac{\pi}{4}'}],[{'-i\frac{\phi}{2}','-i\frac{\pi}{4}'}],[],[{'-i\frac{\phi}{2}','i\frac{\pi}{4}'}],[{'-i\frac{\phi}{2}','i\frac{\pi}{4}'}]};
comp_gates(n).circuit={[1,2;1,1],'\multigate{1}{e^{i\frac{\phi}{2}x x}}','\ghost{e^{i\frac{\phi}{2}x x}}';...
    [1,2;1,1],'\gate{e^{i\frac{\phi}{2}x x}}','\sgate{e^{i\frac{\phi}{2}x x}}{_1_}'};
n=n+1;
%% YY
comp_gates(n).size=2;
comp_gates(n).anc_size=0;
comp_gates(n).step_num=5;
comp_gates(n).steps.index={1,2,[1,2],1,2};
comp_gates(n).steps.gates={'EX','EX','ZZ','EX','EX'};
comp_gates(n).global_phase=exp(-1i*sym('phi')/2);
comp_gates(n).steps.param={[{sym('phi'),sym(-pi/2)}],[{sym('phi'),sym(-pi/2)}],[],[{sym('phi'),sym(pi/2)}],[{sym('phi'),sym(pi/2)}]};
comp_gates(n).circuit_subs={[{'-i\frac{\phi}{2}','i\frac{\pi}{4}'}],[{'-i\frac{\phi}{2}','i\frac{\pi}{4}'}],[],[{'-i\frac{\phi}{2}','-i\frac{\pi}{4}'}],[{'-i\frac{\phi}{2}','-i\frac{\pi}{4}'}]};
comp_gates(n).circuit={[1,2;1,1],'\multigate{1}{e^{-i\frac{\phi}{2}y y}}','\ghost{e^{-i\frac{\phi}{2}y y}}';...
    [1,2;1,1],'\gate{e^{-i\frac{\phi}{2}y y}}','\sgate{e^{-i\frac{\phi}{2}y y}}{_1_}'};
n=n+1;
%% YY_min
comp_gates(n).size=2;
comp_gates(n).anc_size=0;
comp_gates(n).step_num=5;
comp_gates(n).steps.index={1,2,[1,2],1,2};
comp_gates(n).steps.gates={'EX','EX','ZZ_min','EX','EX'};
comp_gates(n).global_phase=exp(-1i*sym('phi')/2);
comp_gates(n).steps.param={[{sym('phi'),-pi/2}],[{sym('phi'),-pi/2}],[],[{sym('phi'),pi/2}],[{sym('phi'),pi/2}]};
comp_gates(n).circuit_subs={[{'-i\frac{\phi}{2}','i\frac{\pi}{4}'}],[{'-i\frac{\phi}{2}','i\frac{\pi}{4}'}],[],[{'-i\frac{\phi}{2}','-i\frac{\pi}{4}'}],[{'-i\frac{\phi}{2}','-i\frac{\pi}{4}'}]};
comp_gates(n).circuit={[1,2;1,1],'\multigate{1}{e^{i\frac{\phi}{2}y y}}','\ghost{e^{i\frac{\phi}{2}y y}}';...
    [1,2;1,1],'\gate{e^{i\frac{\phi}{2}y y}}','\sgate{e^{i\frac{\phi}{2}y y}}{_1_}'};
n=n+1;

%% CX / CNOT1
comp_gates(n).size=2;
comp_gates(n).anc_size=0;
comp_gates(n).step_num=3;
comp_gates(n).steps.index={1,[1,2],1};
comp_gates(n).steps.gates={'H','CZ','H'};
comp_gates(n).steps.param={[],[],[]};
comp_gates(n).circuit={[1,2;1,1],'\targ','\ctrl{_1_}'};
n=n+1;
%% SWAP
comp_gates(n).size=2;
comp_gates(n).anc_size=0;
comp_gates(n).step_num=3;
comp_gates(n).steps.index={[1,2],[2,1],[1,2]};
comp_gates(n).steps.gates={'CX','CX','CX'};
comp_gates(n).steps.param={[],[],[]};
comp_gates(n).circuit={[1,2;1,1],'\qswap','\qswap \qwx'};
n=n+1;
%% CNOT0
comp_gates(n).size=2;
comp_gates(n).anc_size=0;
comp_gates(n).step_num=3;
comp_gates(n).steps.index={2,[1,2],2};
comp_gates(n).steps.gates={'X','CNOT1','X'};
comp_gates(n).steps.param={[],[],[]};
comp_gates(n).circuit={[1,2;1,1],'\targ','\ctrlo{_1_}'};
n=n+1;

%% C0PHASE1
comp_gates(n).size=2;
comp_gates(n).anc_size=0;
comp_gates(n).step_num=3;
comp_gates(n).steps.index={2,[1,2],2};
comp_gates(n).steps.gates={'NOT','CPHASE1','NOT'};
comp_gates(n).steps.param={[],[],[]};
comp_gates(n).circuit={[1,2;1,1],'\gate{\scalebox{0.8}[0.8]{\(\phi_1\)}}','\ctrlo{_1_}'};
n=n+1;

%% C1PHASE0
comp_gates(n).size=2;
comp_gates(n).anc_size=0;
comp_gates(n).step_num=3;
comp_gates(n).steps.index={1,[1,2],1};
comp_gates(n).steps.gates={'NOT','CPHASE1','NOT'};
comp_gates(n).steps.param={[],[],[]};
comp_gates(n).circuit={[1,2;1,1],'\gate{\scalebox{0.8}[0.8]{\(\phi_0\)}}','\ctrl{_1_}'};
n=n+1;

%% C0PHASE0
comp_gates(n).size=2;
comp_gates(n).anc_size=0;
comp_gates(n).step_num=5;
comp_gates(n).steps.index={2,1,[1,2],1,2};
comp_gates(n).steps.gates={'NOT','NOT','CPHASE1','NOT','NOT'};
comp_gates(n).steps.param={[],[],[],[],[]};
comp_gates(n).circuit={[1,2;1,1],'\gate{\scalebox{0.8}[0.8]{\(\phi_0\)}}','\ctrlo{_1_}'};
n=n+1;

[comp_gates(n:n+3).size]=deal(3);
[comp_gates(n:n+3).anc_size]=deal(0);
comp_gates(n).step_num=7;
comp_gates(n+1).step_num=9;
comp_gates(n+2).step_num=11;
comp_gates(n+3).step_num=5;

%% TOFFOLI11
comp_gates(n).steps.index={1,[1,2],[2,3],[1,2],[2,3],[1,3],1};
comp_gates(n).steps.gates={'H','CPHASE1','CNOT1','CPHASE1','CNOT1','CPHASE1','H'};
comp_gates(n).steps.param={[],[sym(pi/2)],[],[sym(-pi/2)],[],[sym(pi/2)],[]};
comp_gates(n).circuit_subs={[],[{'\frac{\phi}{2}','\frac{\pi}{4}'}],[],[{'-i\frac{\phi}{2}','i\frac{\pi}{4}'}],[],[{'\frac{\phi}{2}','\frac{\pi}{4}'}],[]};
comp_gates(n).circuit={[1,2,3;1,1,1],'\targ','\ctrl{_1_}','\ctrl{_1_}'};
n=n+1;
%% TOFFOLI01
comp_gates(n).steps.index={2,1,[1,2],[2,3],[1,2],[2,3],[1,3],1,2};
comp_gates(n).steps.gates={'X','H','CPHASE1','CNOT1','CPHASE1','CNOT1','CPHASE1','H','X'};
comp_gates(n).steps.param={[],[],[sym(pi/2)],[],[-sym(pi/2)],[],[sym(pi/2)],[],[]};
comp_gates(n).circuit_subs={[],[],[{'\frac{\phi}{2}','\frac{\pi}{4}'}],[],[{'-i\frac{\phi}{2}','i\frac{\pi}{4}'}],[],[{'\frac{\phi}{2}','\frac{\pi}{4}'}],[],[]};
comp_gates(n).circuit={[1,2,3;1,1,1],'\targ','\ctrlo{_1_}','\ctrl{_1_}'};
n=n+1;
%% TOFFOLI00
comp_gates(n).steps.index={2,3,1,[1,2],[2,3],[1,2],[2,3],[1,3],1,2,3};
comp_gates(n).steps.gates={'X','X','H','CPHASE1','CNOT1','CPHASE1','CNOT1','CPHASE1','H','X','X'};
comp_gates(n).steps.param={[],[],[],[sym(pi/2)],[],[-sym(pi/2)],[],[sym(pi/2)],[],[],[]};
comp_gates(n).circuit_subs={[],[],[],[{'\frac{\phi}{2}','\frac{\pi}{4}'}],[],[{'-i\frac{\phi}{2}','i\frac{\pi}{4}'}],[],[{'\frac{\phi}{2}','\frac{\pi}{4}'}],[],[],[]};
comp_gates(n).circuit={[1,2,3;1,1,1],'\targ','\ctrlo{_1_}','\ctrlo{_1_}'};
n=n+1;
%% Tophili / CCZ_Phi
comp_gates(n).steps.index={[1,2],[2,3],[1,2],[2,3],[1,3]};
comp_gates(n).steps.gates={,'CPHASE1','CNOT1','CPHASE1','CNOT1','CPHASE1'};
comp_gates(n).steps.param={[sym(phi/2)],[],[-sym(phi/2)],[],[sym(phi/2)]};
comp_gates(n).circuit_subs={[{'\frac{\phi}{2}','\frac{\phi}{4}'}],[],[{'-i\frac{\phi}{2}','i\frac{\phi}{4}'}],[],[{'\frac{\phi}{2}','\frac{\phi}{4}'}]};
comp_gates(n).circuit={[1,2,3;1,1,1],'\gate{e^{-i\frac{\phi}{2}z}}','\ctrl{_1_}','\ctrl{_1_}'};
n=n+1;
%% U
comp_gates(n).size=1;
comp_gates(n).anc_size=0;
comp_gates(n).step_num=8;
comp_gates(n).global_phase=exp(1i*sym('delta'));
comp_gates(n).steps.index={1,1,1};
comp_gates(n).steps.gates={'EZ','EY','EZ'};
comp_gates(n).steps.param={[{sym('phi'),sym('beta','real')}],[{sym('phi'),sym('theta','real')}],[{sym('phi'),sym('alpha','real')}]};
%comp_gates(n).steps.index={1,1,1,1,1,1,1,1};
%comp_gates(n).steps.gates={'EZ','X','EZ','EY','X','EY','EZ','Phi'};
%comp_gates(n).steps.param={[{sym('phi'),1/2*(sym('beta','real')-sym('alpha','real'))}],[],[{sym('phi'),-(sym('beta','real')+sym('alpha','real'))/2}],[{sym('phi'),-1/2*sym('theta','real')}],[],[{sym('phi'),1/2*sym('theta','real')}],[{sym('phi'),sym('alpha','real')}],[{sym('phi'),sym('delta','real')}]};
comp_gates(n).circuit={'\gate{\scalebox{0.8}[0.8]{\(U\)}}'};
n=n+1;
%% CU
comp_gates(n).size=2;
comp_gates(n).anc_size=0;
comp_gates(n).step_num=8;
comp_gates(n).global_phase=exp(1i*1/4*sym('delta'));
comp_gates(n).steps.index={1,[1,2],1,1,[1,2],1,1,2};
comp_gates(n).steps.gates={'EZ','CNOT','EZ','EY','CNOT','EY','EZ','EZ'};
comp_gates(n).steps.param={[{sym('phi'),1/2*(sym('beta','real')-sym('alpha','real'))}],[],[{sym('phi'),-(sym('beta','real')+sym('alpha','real'))/2}],[{sym('phi'),-1/2*sym('theta','real')}],[],[{sym('phi'),1/2*sym('theta','real')}],[{sym('phi'),sym('alpha','real')}],[{sym('phi'),1/2*sym('delta','real')}]};
comp_gates(n).circuit={[1,2;1,1],'\gate{\scalebox{0.8}[0.8]{\(U\)}}','\ctrl{_1_}'};
n=n+1;
%%C0U
comp_gates(n).size=2;
comp_gates(n).anc_size=0;
comp_gates(n).step_num=8;
comp_gates(n).global_phase=exp(1i*1/4*sym('delta'));
comp_gates(n).steps.index={1,[1,2],1,1,[1,2],1,1,2};
comp_gates(n).steps.gates={'EZ','CNOT0','EZ','EY','CNOT0','EY','EZ','EZ'};
comp_gates(n).steps.param={[{sym('phi'),1/2*(sym('beta','real')-sym('alpha','real'))}],[],[{sym('phi'),-(sym('beta','real')+sym('alpha','real'))/2}],[{sym('phi'),-1/2*sym('theta','real')}],[],[{sym('phi'),1/2*sym('theta','real')}],[{sym('phi'),sym('alpha','real')}],[{sym('phi'),-1/2*sym('delta','real')}]};
comp_gates(n).circuit={[1,2;1,1],'\gate{\scalebox{0.8}[0.8]{\(U\)}}','\ctrlo{_1_}'};


%[ matrix anc_matrix ] = Gate2Matrix( elem_gates,comp_gates, comp_gates(9),2)
for i=1:length(comp_gates)
    [ matrix anc_matrix ] = Gate2Matrix( elem_gates,comp_gates, comp_gates(i),2);
    if isa(matrix,'sym')
        comp_gates(i).matrix=simplify(matrix);
    else
        comp_gates(i).matrix=matrix;
    end
    
    curr_name=comp_gates(i).names;
    if isa(curr_name,'char')
        curr_name2=curr_name;
    else
        curr_name2=curr_name{1};
    end
    %Print the resulting Gates Matrix
    fprintf(['\n<strong>Gate:</strong> ' curr_name2 '\n']);
    [M index M_mode]=Sparse2Square(matrix,{'fock',1,comp_gates(i).size});
    fprintf(M_mode)
end

save('comp_gates.mat','comp_gates');