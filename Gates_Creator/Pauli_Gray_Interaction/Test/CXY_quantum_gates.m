%CXY_quantum_gates.m
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

%-------------------------------------------------------------------------

name='XX_YY';
size=2;
anc_size=0;
gate_steps={'XX','YY'};
index_steps={[1 2],[1 2]};
param_steps={[],[]};
gate_string='e^{-i\frac{\phi}{2}(x x+y y)}';

cxy_gates=Create_Comp_Gate(name,size,anc_size,gate_steps,index_steps,param_steps);
cxy_gates=Add_Matrix2Gate(cxy_gates,elem_gates,comp_gates);
fprintf([cxy_gates.names ':\n'])
FockPrint(cxy_gates.matrix)
cxy_gates=Generate_Gate_Circuit(cxy_gates,gate_string,1:size,[],[]);
Gates2Circ2Preview(cxy_gates,elem_gates,comp_gates,name,'gate=expansion')

%-------------------------------------------------------------------------
name='C_ZZ';
size=3;
anc_size=0;
gate_steps={'X','TOPHILI','X','X','TOPHILI','X'};
index_steps={1,[1 2 3],1,2,[1 2 3],2};
param_steps={[],[],[],[],[],[]};
gate_string='e^{-i\frac{\phi}{2}z z}';


comp_gates=Comp_Gate_Merger(comp_gates,cxy_gates);

gate=Create_Comp_Gate(name,size,anc_size,gate_steps,index_steps,param_steps);
gate=Add_Matrix2Gate(gate,elem_gates,comp_gates);
fprintf([gate.names ':\n'])
FockPrint(gate.matrix)
gate=Generate_Gate_Circuit(gate,gate_string,1:size,[1],[3]); 
Gates2Circ2Preview(gate,elem_gates,comp_gates,name,'gate=expansion')
cxy_gates=Comp_Gate_Merger(cxy_gates,gate);

%-------------------------------------------------------------------------
name='C_XX';
size=3;
anc_size=0;
gate_steps={'EY','EY','C_ZZ','EY','EY','P1'};
index_steps={1,2,[1,2,3],1,2,3};
param_steps={[pi/2],[pi/2],[],[-pi/2],[-pi/2],[-sym('phi')/2]};
circuit_subs={[{'-i\frac{\phi}{2}','-i\frac{\pi}{4}'}],[{'-i\frac{\phi}{2}','-i\frac{\pi}{4}'}],[],[{'-i\frac{\phi}{2}','i\frac{\pi}{4}'}],[{'-i\frac{\phi}{2}','i\frac{\pi}{4}'}],[]};
gate_string='e^{-i\frac{\phi}{2}x x}';


comp_gates=Comp_Gate_Merger(comp_gates,cxy_gates);

gate=Create_Comp_Gate(name,size,anc_size,gate_steps,index_steps,param_steps);
gate=Add_Matrix2Gate(gate,elem_gates,comp_gates);
fprintf([gate.names ':\n'])
FockPrint(gate.matrix)
gate=Generate_Gate_Circuit(gate,gate_string,1:size,[1],[3],circuit_subs); 
Gates2Circ2Preview(gate,elem_gates,comp_gates,name,'gate=expansion')
cxy_gates=Comp_Gate_Merger(cxy_gates,gate);

%-------------------------------------------------------------------------
name='C_YY';
size=3;
anc_size=0;
gate_steps={'EX','EX','C_ZZ','EX','EX','P1'};
index_steps={1,2,[1,2,3],1,2,3};
param_steps={[-pi/2],[-pi/2],[],[pi/2],[pi/2],[-sym('phi')/2]};
circuit_subs={[{'-i\frac{\phi}{2}','i\frac{\pi}{4}'}],[{'-i\frac{\phi}{2}','i\frac{\pi}{4}'}],[],[{'-i\frac{\phi}{2}','-i\frac{\pi}{4}'}],[{'-i\frac{\phi}{2}','-i\frac{\pi}{4}'}],[]};
gate_string='e^{-i\frac{\phi}{2}y y}';


comp_gates=Comp_Gate_Merger(comp_gates,cxy_gates);

gate=Create_Comp_Gate(name,size,anc_size,gate_steps,index_steps,param_steps);
gate=Add_Matrix2Gate(gate,elem_gates,comp_gates);
fprintf([gate.names ':\n'])
FockPrint(gate.matrix)
gate=Generate_Gate_Circuit(gate,gate_string,1:size,[1],[3],circuit_subs); 
Gates2Circ2Preview(gate,elem_gates,comp_gates,name,'gate=expansion')
cxy_gates=Comp_Gate_Merger(cxy_gates,gate);

%-------------------------------------------------------------------------

name='C_XX_YY';
size=3;
anc_size=0;
gate_steps={'C_XX','C_YY'};
index_steps={[1 2 3],[1 2 3]};
param_steps={[],[]};
gate_string='e^{-i\frac{\phi}{2}(x x+y y)}';


comp_gates=Comp_Gate_Merger(comp_gates,cxy_gates);

gate=Create_Comp_Gate(name,size,anc_size,gate_steps,index_steps,param_steps);
gate=Add_Matrix2Gate(gate,elem_gates,comp_gates);
fprintf([gate.names ':\n'])
FockPrint(gate.matrix)
gate=Generate_Gate_Circuit(gate,gate_string,1:size,[1],[3]);
Gates2Circ2Preview(gate,elem_gates,comp_gates,name,'gate=expansion')
cxy_gates=Comp_Gate_Merger(cxy_gates,gate);

%-------------------------------------------------------------------------

name='Cx_XX_YY';
size=3;
anc_size=0;
gate_steps={'X','C_XX','C_YY','X'};
index_steps={1,[1 2 3],[1 2 3],1};
param_steps={[],[],[],[]};
gate_string='e^{-i\frac{\phi}{2}(x x-y y)}';


comp_gates=Comp_Gate_Merger(comp_gates,cxy_gates);

gate=Create_Comp_Gate(name,size,anc_size,gate_steps,index_steps,param_steps);
gate=Add_Matrix2Gate(gate,elem_gates,comp_gates);
fprintf([gate.names ':\n'])
FockPrint(gate.matrix)
gate=Generate_Gate_Circuit(gate,gate_string,1:size,[1],[3]);
Gates2Circ2Preview(gate,elem_gates,comp_gates,name,'gate=expansion')
cxy_gates=Comp_Gate_Merger(cxy_gates,gate);


%% -----------------------------------------------------------------------

filename=Save_Gate_list('cxy_gates','cxy_gates',cxy_gates);
fprintf('File has been saved\n');