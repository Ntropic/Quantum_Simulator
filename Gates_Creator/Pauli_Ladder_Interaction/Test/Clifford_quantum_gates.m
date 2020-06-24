%Clifford_quantum_gates.m
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

%Load gates
load('../../../Gates_Table/elem_gates.mat','-mat')
load('../../../Gates_Table/comp_gates.mat','-mat')

phi=sym('phi');
%-------------------------------------------------------------------------

name='XXXX';
size=4;
anc_size=0;
gate_steps={'EY','EY','EY','ZZ_min','ZZ_min','ZZ_min','EY','ZZ','ZZ','ZZ','EY','EY','EY'};
index_steps={2,3,4,[1 2],[1 3],[1 4],1,[1 4],[1 3],[1 2],4,3,2};
param_steps={[-pi/2],[-pi/2],[-pi/2],[pi/2],[pi/2],[pi/2],[],[pi/2],[pi/2],[pi/2],[pi/2],[pi/2],[pi/2]};
%gate_string='e^{-i\frac{\phi}{2}(x^{\otimes 4)}}';
gate_string='e^{-i\frac{\phi}{2}xxxx}';

clifford_gates=Create_Comp_Gate(name,size,anc_size,gate_steps,index_steps,param_steps);
clifford_gates=Add_Matrix2Gate(clifford_gates,elem_gates,comp_gates);
fprintf([clifford_gates.names ':\n'])
FockSparsePrint(clifford_gates.matrix)
clifford_gates=Generate_Gate_Circuit(clifford_gates,gate_string,1:size,[],[]);
Gates2Circ2Preview(clifford_gates,elem_gates,comp_gates,name,'gate=expansion')

%-------------------------------------------------------------------------

name='XXXX_min';
size=4;
anc_size=0;
gate_steps={'EY','EY','EY','ZZ_min','ZZ_min','ZZ_min','EY','ZZ','ZZ','ZZ','EY','EY','EY'};
index_steps={2,3,4,[1 2],[1 3],[1 4],1,[1 4],[1 3],[1 2],4,3,2};
param_steps={[pi/2],[pi/2],[pi/2],[pi/2],[pi/2],[pi/2],[],[pi/2],[pi/2],[pi/2],[-pi/2],[-pi/2],[-pi/2]};
%gate_string='e^{i\frac{\phi}{2}(x^{\otimes 4)}}';
gate_string='e^{i\frac{\phi}{2}xxxx}';

comp_gates=Comp_Gate_Merger(comp_gates,clifford_gates);

gate=Create_Comp_Gate(name,size,anc_size,gate_steps,index_steps,param_steps);
gate=Add_Matrix2Gate(gate,elem_gates,comp_gates);
fprintf([gate.names ':\n'])
FockSparsePrint(gate.matrix)
gate=Generate_Gate_Circuit(gate,gate_string,1:size,[],[]); 
Gates2Circ2Preview(gate,elem_gates,comp_gates,name,'gate=expansion')
clifford_gates=Comp_Gate_Merger(clifford_gates,gate);

%% -----------------------------------------------------------------------

name='YYYY';
size=4;
anc_size=0;
gate_steps={'EX','EX','EX','ZZ_min','ZZ_min','ZZ_min','EX','ZZ','ZZ','ZZ','EX','EX','EX'};
index_steps={2,3,4,[1 2],[1 3],[1 4],1,[1 4],[1 3],[1 2],4,3,2};
param_steps={[-pi/2],[-pi/2],[-pi/2],[pi/2],[pi/2],[pi/2],[],[pi/2],[pi/2],[pi/2],[pi/2],[pi/2],[pi/2]};
%gate_string='e^{-i\frac{\phi}{2}(Y^{\otimes 4)}}';
gate_string='e^{-i\frac{\phi}{2}yyyy}';

comp_gates=Comp_Gate_Merger(comp_gates,clifford_gates);

gate=Create_Comp_Gate(name,size,anc_size,gate_steps,index_steps,param_steps);
gate=Add_Matrix2Gate(gate,elem_gates,comp_gates);
fprintf([gate.names ':\n'])
FockSparsePrint(gate.matrix)
gate=Generate_Gate_Circuit(gate,gate_string,1:size,[],[]); 
Gates2Circ2Preview(gate,elem_gates,comp_gates,name,'gate=expansion')
clifford_gates=Comp_Gate_Merger(clifford_gates,gate);

%% -----------------------------------------------------------------------
name='YYYY_min';
size=4;
anc_size=0;
gate_steps={'EX','EX','EX','ZZ_min','ZZ_min','ZZ_min','EX','ZZ','ZZ','ZZ','EX','EX','EX'};
index_steps={2,3,4,[1 2],[1 3],[1 4],1,[1 4],[1 3],[1 2],4,3,2};
param_steps={[pi/2],[pi/2],[pi/2],[pi/2],[pi/2],[pi/2],[],[pi/2],[pi/2],[pi/2],[-pi/2],[-pi/2],[-pi/2]};
%gate_string='e^{i\frac{\phi}{2}(Y^{\otimes 4)}}';
gate_string='e^{i\frac{\phi}{2}yyyy}';

comp_gates=Comp_Gate_Merger(comp_gates,clifford_gates);

gate=Create_Comp_Gate(name,size,anc_size,gate_steps,index_steps,param_steps);
gate=Add_Matrix2Gate(gate,elem_gates,comp_gates);
fprintf([gate.names ':\n'])
FockSparsePrint(gate.matrix)
gate=Generate_Gate_Circuit(gate,gate_string,1:size,[],[]); 
Gates2Circ2Preview(gate,elem_gates,comp_gates,name,'gate=expansion')
clifford_gates=Comp_Gate_Merger(clifford_gates,gate);

%% -----------------------------------------------------------------------

name='XYXY';
size=4;
anc_size=0;
gate_steps={'EY','EX','EY','ZZ_min','ZZ_min','ZZ_min','EX','ZZ','ZZ','ZZ','EY','EX','EY'};
index_steps={2,3,4,[1 2],[1 3],[1 4],1,[1 4],[1 3],[1 2],4,3,2};
param_steps={[-pi/2],[-pi/2],[-pi/2],[pi/2],[pi/2],[pi/2],[],[pi/2],[pi/2],[pi/2],[pi/2],[pi/2],[pi/2]};
gate_string='e^{-i\frac{\phi}{2}xyxy}';

comp_gates=Comp_Gate_Merger(comp_gates,clifford_gates);

gate=Create_Comp_Gate(name,size,anc_size,gate_steps,index_steps,param_steps);
gate=Add_Matrix2Gate(gate,elem_gates,comp_gates);
fprintf([gate.names ':\n'])
FockSparsePrint(gate.matrix)
gate=Generate_Gate_Circuit(gate,gate_string,1:size,[],[]); 
Gates2Circ2Preview(gate,elem_gates,comp_gates,name,'gate=expansion')
clifford_gates=Comp_Gate_Merger(clifford_gates,gate);

%% -----------------------------------------------------------------------
name='XYXY_min';
size=4;
anc_size=0;
gate_steps={'EY','EX','EY','ZZ_min','ZZ_min','ZZ_min','EX','ZZ','ZZ','ZZ','EY','EX','EY'};
index_steps={2,3,4,[1 2],[1 3],[1 4],1,[1 4],[1 3],[1 2],4,3,2};
param_steps={[pi/2],[pi/2],[pi/2],[pi/2],[pi/2],[pi/2],[],[pi/2],[pi/2],[pi/2],[-pi/2],[-pi/2],[-pi/2]};
gate_string='e^{i\frac{\phi}{2}xyxy}';

comp_gates=Comp_Gate_Merger(comp_gates,clifford_gates);

gate=Create_Comp_Gate(name,size,anc_size,gate_steps,index_steps,param_steps);
gate=Add_Matrix2Gate(gate,elem_gates,comp_gates);
fprintf([gate.names ':\n'])
FockSparsePrint(gate.matrix)
gate=Generate_Gate_Circuit(gate,gate_string,1:size,[],[]); 
Gates2Circ2Preview(gate,elem_gates,comp_gates,name,'gate=expansion')
clifford_gates=Comp_Gate_Merger(clifford_gates,gate);

%% -----------------------------------------------------------------------

name='YXYX';
size=4;
anc_size=0;
gate_steps={'EX','EY','EX','ZZ_min','ZZ_min','ZZ_min','EY','ZZ','ZZ','ZZ','EX','EY','EX'};
index_steps={2,3,4,[1 2],[1 3],[1 4],1,[1 4],[1 3],[1 2],4,3,2};
param_steps={[-pi/2],[-pi/2],[-pi/2],[pi/2],[pi/2],[pi/2],[],[pi/2],[pi/2],[pi/2],[pi/2],[pi/2],[pi/2]};
gate_string='e^{-i\frac{\phi}{2}yxyx}';

comp_gates=Comp_Gate_Merger(comp_gates,clifford_gates);

gate=Create_Comp_Gate(name,size,anc_size,gate_steps,index_steps,param_steps);
gate=Add_Matrix2Gate(gate,elem_gates,comp_gates);
fprintf([gate.names ':\n'])
FockSparsePrint(gate.matrix)
gate=Generate_Gate_Circuit(gate,gate_string,1:size,[],[]); 
Gates2Circ2Preview(gate,elem_gates,comp_gates,name,'gate=expansion')
clifford_gates=Comp_Gate_Merger(clifford_gates,gate);

%% -----------------------------------------------------------------------
name='YXYX_min';
size=4;
anc_size=0;
gate_steps={'EX','EY','EX','ZZ_min','ZZ_min','ZZ_min','EY','ZZ','ZZ','ZZ','EX','EY','EX'};
index_steps={2,3,4,[1 2],[1 3],[1 4],1,[1 4],[1 3],[1 2],4,3,2};
param_steps={[pi/2],[pi/2],[pi/2],[pi/2],[pi/2],[pi/2],[],[pi/2],[pi/2],[pi/2],[-pi/2],[-pi/2],[-pi/2]};
gate_string='e^{i\frac{\phi}{2}yxyx}';

comp_gates=Comp_Gate_Merger(comp_gates,clifford_gates);

gate=Create_Comp_Gate(name,size,anc_size,gate_steps,index_steps,param_steps);
gate=Add_Matrix2Gate(gate,elem_gates,comp_gates);
fprintf([gate.names ':\n'])
FockSparsePrint(gate.matrix)
gate=Generate_Gate_Circuit(gate,gate_string,1:size,[],[]); 
Gates2Circ2Preview(gate,elem_gates,comp_gates,name,'gate=expansion')
clifford_gates=Comp_Gate_Merger(clifford_gates,gate);

%% -----------------------------------------------------------------------
name='L_{+}';
size=4;
anc_size=0;
gate_steps={'XXXX','YXYX','YYYY','XYXY'};
index_steps={[1 2 3 4],[1 2 3 4],[1 2 3 4],[1 2 3 4]};
param_steps={[{phi,1/4*phi}],[{phi,1/4*phi}],[{phi,1/4*phi}],[{phi,1/4*phi}]};
gate_string='L^4_{-\phi}';
circ_subs={[{'\frac{\phi}{2}','\frac{\phi}{4}'}],[{'\frac{\phi}{2}','\frac{\phi}{4}'}],[{'\frac{\phi}{2}','\frac{\phi}{4}'}],[{'\frac{\phi}{2}','\frac{\phi}{4}'}]};

comp_gates=Comp_Gate_Merger(comp_gates,clifford_gates);

gate=Create_Comp_Gate(name,size,anc_size,gate_steps,index_steps,param_steps);
gate=Add_Matrix2Gate(gate,elem_gates,comp_gates);
fprintf([gate.names ':\n'])
FockSparsePrint(gate.matrix)
gate=Generate_Gate_Circuit(gate,gate_string,1:size,[],[],circ_subs); 
Gates2Circ2Preview(gate,elem_gates,comp_gates,name,'gate=expansion')
clifford_gates=Comp_Gate_Merger(clifford_gates,gate);

%% -----------------------------------------------------------------------

name='L_{-}';
size=4;
anc_size=0;
gate_steps={'XXXX_min','YXYX_min','YYYY_min','XYXY_min'};
index_steps={[1 2 3 4],[1 2 3 4],[1 2 3 4],[1 2 3 4]};
param_steps={[{phi,1/4*phi}],[{phi,1/4*phi}],[{phi,1/4*phi}],[{phi,1/4*phi}]};
gate_string='L^4_{\phi}';
circ_subs={[{'\frac{\phi}{2}','\frac{\phi}{4}'}],[{'\frac{\phi}{2}','\frac{\phi}{4}'}],[{'\frac{\phi}{2}','\frac{\phi}{4}'}],[{'\frac{\phi}{2}','\frac{\phi}{4}'}]};

comp_gates=Comp_Gate_Merger(comp_gates,clifford_gates);

gate=Create_Comp_Gate(name,size,anc_size,gate_steps,index_steps,param_steps);
gate=Add_Matrix2Gate(gate,elem_gates,comp_gates);
fprintf([gate.names ':\n'])
FockSparsePrint(gate.matrix)
gate=Generate_Gate_Circuit(gate,gate_string,1:size,[],[],circ_subs); 
Gates2Circ2Preview(gate,elem_gates,comp_gates,name,'gate=expansion')
clifford_gates=Comp_Gate_Merger(clifford_gates,gate);

%% -----------------------------------------------------------------------

filename=Save_Gate_list('clifford_gates','clifford_gates',clifford_gates);
fprintf('File has been saved\n');