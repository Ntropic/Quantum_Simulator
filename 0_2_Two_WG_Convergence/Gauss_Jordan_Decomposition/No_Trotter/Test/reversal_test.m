%reversal_test.m
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
load('../../../../Gates_Table/elem_gates.mat','-mat')
load('../../../../Gates_Table/comp_gates.mat','-mat')

%Photon number
n=1;


%Iteration time
t=pi/4;

%% Approximate decomposition
mode=4;

if mode==1
    gate2=Create_Comp_Gate( 'testgate',2,0,{'CU'},{[1,2]})
    %gate2=Create_Comp_Gate( 'testgate',2,0,{'CNOT','CNOT'},{[2,1],[1,2]})
    other_gates2=Comp_Gate_Merger(comp_gates,gate2);

    U_approx=Gate2Matrix(elem_gates,comp_gates,gate2)

    other_gates2=Gates_Tables_Prepare(other_gates2,elem_gates);   

    ind=Gate_Index_by_Name(gate2.names,elem_gates,other_gates2);
    U_approx2=simplify(full(Expansion2Matrix(other_gates2(ind(2)))),10)


    simplify(U_approx-U_approx2,10)
elseif mode==2
    param=[{sym('alpha'),-(3*pi)/2,sym('beta'),(3*pi)/2,sym('delta'),0,sym('theta'),pi}];
    gate2=Create_Comp_Gate( 'testgate',2,0,{'CU'},{[1,2]},{param})
    
    other_gates2=Comp_Gate_Merger(comp_gates,gate2);

    U_approx=Gate2Matrix(elem_gates,comp_gates,gate2)

    other_gates2=Gates_Tables_Prepare(other_gates2,elem_gates);   

    ind=Gate_Index_by_Name(gate2.names,elem_gates,other_gates2);
    U_approx2=full(Expansion2Matrix(other_gates2(ind(2))))


   	U_approx-U_approx2
elseif mode==3
    gate2=Create_Comp_Gate( 'testgate',1,0,{'U'},{[1]})
    
    other_gates2=Comp_Gate_Merger(comp_gates,gate2);

    U_approx=simplify(Gate2Matrix(elem_gates,comp_gates,gate2),10)

    other_gates2=Gates_Tables_Prepare(other_gates2,elem_gates);   

    ind=Gate_Index_by_Name(gate2.names,elem_gates,other_gates2);
    U_approx2=simplify(full(Expansion2Matrix(other_gates2(ind(2)))),10)


   	simplify(U_approx-U_approx2,10)
elseif mode==4
    param=[{sym('alpha'),-(3*pi)/2,sym('beta'),(3*pi)/2,sym('delta'),0,sym('theta'),pi}];
    gate2=Create_Comp_Gate( 'testgate',1,0,{'U'},{[1]},{param})
    
    other_gates2=Comp_Gate_Merger(comp_gates,gate2);

    U_approx=Gate2Matrix(elem_gates,comp_gates,gate2)

    other_gates2=Gates_Tables_Prepare(other_gates2,elem_gates);   

    ind=Gate_Index_by_Name(gate2.names,elem_gates,other_gates2);
    U_approx2=full(Expansion2Matrix(other_gates2(ind(2))))


   	U_approx-U_approx2
elseif mode==5
    gate2=Create_Comp_Gate( 'testgate',1,0,{'EZ'},{[1]})
    
    other_gates2=Comp_Gate_Merger(comp_gates,gate2);

    U_approx=simplify(Gate2Matrix(elem_gates,comp_gates,gate2),10)

    other_gates2=Gates_Tables_Prepare(other_gates2,elem_gates);   

    ind=Gate_Index_by_Name(gate2.names,elem_gates,other_gates2);
    U_approx2=simplify(full(Expansion2Matrix(other_gates2(ind(2)))),10)


   	simplify(U_approx-U_approx2,10)
elseif mode==6
    param=[{sym('phi'),sym('beta','real')}];
    gate=Create_Comp_Gate( 'test',1,0,{'EZ'},{[1]},{param})
    param2=[{sym('alpha'),-(3*pi)/2,sym('beta'),(3*pi)/2,sym('delta'),0,sym('theta'),pi}];
    gate2=Create_Comp_Gate( 'test2',1,0,{'test'},{[1]},{param2})
    
    other_gates2=Comp_Gate_Merger(comp_gates,gate,gate2);

    U_approx=Gate2Matrix(elem_gates,other_gates2,gate2)

    other_gates2=Gates_Tables_Prepare(other_gates2,elem_gates);   

    ind=Gate_Index_by_Name(gate2.names,elem_gates,other_gates2);
    U_approx2=full(Expansion2Matrix(other_gates2(ind(2))))


   	U_approx-U_approx2
end