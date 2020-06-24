%test_simple_errors.m
clc;
clear all;
close all;

p=pwd;
if any(strfind(p,'\'));
    elem=strsplit(p,'\');
else
    elem=strsplit(p,'/');
end
shortened=fullfile(elem{1:end-2});
addpath(genpath(shortened));
%Load gates
load('../../Gates_Table/elem_gates.mat','-mat')
load('../../Gates_Table/comp_gates.mat','-mat')

single_error=1;
double_error=10;

[elem_gates, comp_gates]=Simple_Errors(single_error,double_error, elem_gates,comp_gates);

fprintf(['Elementary Gates: \n'])

for i=1:length(elem_gates(:))
    if isa(elem_gates(i).names,'cell')
        fprintf(['  ' elem_gates(i).names{1} '\n'])
    else
        fprintf(['  ' elem_gates(i).names '\n'])
    end
    fprintf(['    ' num2str(elem_gates(i).cost) '\n'])
    e(i)=elem_gates(i).cost;
end

fprintf('-------------------------------------------------------\n')

for i=1:length(comp_gates(:))
    if isa(comp_gates(i).names,'cell')
        fprintf(['  ' comp_gates(i).names{1} '\n'])
    else
        fprintf(['  ' comp_gates(i).names '\n'])
    end
    fprintf(['    ' num2str(comp_gates(i).cost) '\n'])
    c(i)=comp_gates(i).cost;
end

plot(e)

figure()
plot(c)