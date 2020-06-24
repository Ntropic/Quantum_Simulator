%pauli_gate_tex.m
clear all;
close all;
clc;

mode=5;
n=5;

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

%Create/Load Controlled Unitaries for different sizes
PU_gates=CreateCompleteSet_Cn_U2(n-1,elem_gates,comp_gates,[1,1,0]);

if mode==1
    for j=1:3
        if j==1
            gates_opt=elem_gates(1+j);
            file_name=[gates_opt.names{1}];
        else
            gates_opt=elem_gates(1+j);
            file_name=[gates_opt.names];
        end

        mode='expansion';
        depth=1;
        sep_dist=15;

        file_names=Gates2Tex(gates_opt,elem_gates,comp_gates,file_name,mode,depth,sep_dist);
        %Plot the results
        for i=1:length(file_names)
            tex2pdf2preview(file_names{i},1,800);
        end
    end
elseif mode==2
    gates_opt=elem_gates(12);
    file_name=[gates_opt.names{3}];

    mode='gate';
    depth=1;
    sep_dist=15;
    
    file_names=Gates2Tex(gates_opt,elem_gates,comp_gates,file_name,mode,depth,sep_dist);
    %Plot the results
    for i=1:length(file_names)
        tex2pdf2preview(file_names{i},1,800);
    end    
    %[filename]=SparseFull2Tex(gates_opt.matrix,['C_not_mat'],'tikz_gray_fock_dense_def_subs_var_width');
elseif mode==3
    gates_opt=comp_gates(4);
    file_name=[gates_opt.names 'gate_exp'];

    mode='gate=expansion';
    depth=1;
    sep_dist=15;
    
    file_names=Gates2Tex(gates_opt,elem_gates,comp_gates,file_name,mode,depth,sep_dist);
    %Plot the results
    for i=1:length(file_names)
        tex2pdf2preview(file_names{i},1,800);
    end    
    %[filename]=SparseFull2Tex(gates_opt.matrix,['Toffoli_mat'],'tikz_gray_fock_dense_def_subs_var_width');
elseif mode==4
    gates_opt=comp_gates(8);
    file_name=[gates_opt.names{1} 'gate_exp'];

    mode='gate=expansion';
    depth=1;
    sep_dist=15;
    
    file_names=Gates2Tex(gates_opt,elem_gates,comp_gates,file_name,mode,depth,sep_dist);
    %Plot the results
    for i=1:length(file_names)
        tex2pdf2preview(file_names{i},1,800);
    end    
    %[filename]=SparseFull2Tex(gates_opt.matrix,['CU_mat'],'tikz_gray_fock_dense_def_subs_var_width');
elseif mode==5
    gates_opt=PU_gates(8);
    file_name=[gates_opt.names 'gate_exp'];

    mode='gate=expansion';
    depth=1;
    sep_dist=15;
    
    file_names=Gates2Tex(gates_opt,elem_gates,comp_gates,file_name,mode,depth,sep_dist);
    %Plot the results
    for i=1:length(file_names)
        tex2pdf2preview(file_names{i},1,800);
    end    
    [filename]=SparseFull2Tex(gates_opt.matrix,['PU_mat'],'tikz_gray_fock_dense_def_subs_var_width');

elseif mode==6  % Test Multigates (together)
    gate_name='U';
    gates_opt=Gate_by_Name(gate_name,elem_gates,comp_gates);
    file_name=[gate_name 'gate_exp'];

    mode='gate=expansion';
    depth=1;
    sep_dist=15;
    
    file_names=Gates2Tex(gates_opt,elem_gates,comp_gates,file_name,mode,depth,sep_dist);
    %Plot the results
    for i=1:length(file_names)
        tex2pdf2preview(file_names{i},1,800);
    end    
    %[filename]=SparseFull2Tex(gates_opt.matrix,['PU_mat'],'tikz_gray_fock_dense_def_subs_var_width');
end