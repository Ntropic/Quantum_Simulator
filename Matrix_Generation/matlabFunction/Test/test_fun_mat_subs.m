%test_fun_mat_subs.m
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


%% ------------------------------------------------------------------------
i=10;
mode=2;

if mode==1  %Not cell - for the replacement data (fun_var_subs)
    fun_mat=elem_gates(i).fun_mat
    fun_vars=elem_gates(i).fun_vars
    fun_var_subs=[sym('x/l')];

    [fun_mat fun_vars]=fun_mat_subs(fun_mat,fun_vars,fun_var_subs)


    fun_var_subs=[sym('a'),pi/2];

    [fun_mat_new fun_vars_new]=fun_mat_subs(fun_mat,fun_vars,fun_var_subs)
    
elseif mode==2    %Cell - for the replacement data (fun_var_subs)
    fun_mat=elem_gates(i).fun_mat
    fun_vars=elem_gates(i).fun_vars
    fun_var_subs=[{sym('phi'),sym('beta','real')}];
    %fun_var_subs=[{sym('phi'),sym('beta','real') - sym('alpha','real')'}];
    
    [fun_mat2 fun_vars2]=fun_mat_subs(fun_mat,fun_vars,fun_var_subs)
    
    param=[{sym('alpha'),-(3*pi)/2,sym('beta'),(3*pi)/2,sym('delta'),0,sym('theta'),pi}];

    fun_mat2((3*pi)/2)
    [fun_mat_new fun_vars_new]=fun_mat_subs(fun_mat2,fun_vars2,param)
elseif mode==3
    fun_mat=elem_gates(i).fun_mat
    fun_vars=elem_gates(i).fun_vars
    fun_var_subs=[sym('pi/2')];
    
    [fun_mat fun_vars]=fun_mat_subs(fun_mat,fun_vars,fun_var_subs)
elseif mode==4
    fun_mat=elem_gates(i).fun_mat
    fun_vars=elem_gates(i).fun_vars
    fun_var_subs=[sym('phi/2')];
    
    [fun_mat fun_vars]=fun_mat_subs(fun_mat,fun_vars,fun_var_subs)
end
    