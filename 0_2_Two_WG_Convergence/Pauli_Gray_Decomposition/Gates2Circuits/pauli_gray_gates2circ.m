%pauli_gray_gates2circ.m
clear all;
close all;
clc;

mode=2;

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
load('../cxy_gates.mat','-mat')
other_gates=Comp_Gate_Merger(comp_gates,cxy_gates);



if mode==1
    curr=comp_gates(7);
     depth=1;
    file_name=[curr.names{1} 'gate_exp_' num2str(depth)];

    mode='gate=expansion';

    sep_dist=25;

    file_names=Gates2Tex(curr,elem_gates,other_gates,file_name,mode,depth,sep_dist);
    %Plot the results
    for i=1:length(file_names)
        tex2pdf2preview(file_names{i},1,800);
    end
elseif mode==2
    for i=1:5
        depth=1;
        if i==1
            gates_opt=cxy_gates(2);
            file_name=[gates_opt.names 'gate_exp_' num2str(depth)];
        elseif i==2
            gates_opt=cxy_gates(3);
            file_name=[gates_opt.names 'gate_exp_' num2str(depth)];
        elseif i==3
            gates_opt=cxy_gates(4);
            file_name=[gates_opt.names 'gate_exp_' num2str(depth)];
        elseif i==4
            gates_opt=comp_gates(16);
            file_name=[gates_opt.names{1} 'gate_exp_' num2str(depth)];
        elseif i==5
            gates_opt=comp_gates(13);
            file_name=[gates_opt.names 'gate_exp_' num2str(depth)];
        end


        mode='gate=expansion';

        sep_dist=25;

        file_names=Gates2Tex(gates_opt,elem_gates,other_gates,file_name,mode,depth,sep_dist);
        %Plot the results
        for i=1:length(file_names)
            tex2pdf2preview(file_names{i},1,800);
        end  
    end
elseif mode==3
    %gates_opt=clifford_gates(9);

    new_gate=Add_Gate(new_gate,{'L_{+}'},{[2,3,4,5]},{[]});

    depth=1;
    file_name=[new_gate.names 'gate_exp_' num2str(depth)];

    mode='expansion';

    sep_dist=25;

    file_names=Gates2Tex(new_gate,elem_gates,other_gates,file_name,mode,depth,sep_dist);
    %Plot the results
    for i=1:length(file_names)
        tex2pdf2preview(file_names{i},1,800);
    end  
end