%test_sparsesym.m
clc;
clear all;
close all;

p=pwd;
elem=strsplit(p,'\');
shortened=strjoin(elem((1:length(elem)-1)),'\');
addpath(genpath(shortened))

%Generate symbolic matrices
prob=0.25;
test=3;
repeat=20;

dim1=[4,4];
dim2=[4,4];

symlist=[sym('x'),sym('y'),sym('c'),sym('d')];

%Create sparsesym


%% TESTS
if any(test==1)
    %Check if reverse works
    
    fprintf('Create sparsesym from sym and reverse:\n')
    for i=1:repeat
        [ A,B ] = testrandomsymgenerator( symlist,dim1,dim2,prob );
        a=sparsesym(A);
        b=sparsesym(B);
        A_rec(i)=max(max(full(a)-A));
        B_rec(i)=max(max(full(b)-B));
    end
    max(A_rec)
    max(B_rec)
elseif any(test==2)
    fprintf('Test addition:\n')
    
    for i=1:repeat
        [ A,B ] = testrandomsymgenerator( symlist,dim1,dim2,prob );
        a=sparsesym(A);
        b=sparsesym(B);
        diff(i)=max(max(full(a+b)-(A+B)));
    end
    max(diff)
    fprintf('Test subtraction:\n')
    
   	for i=1:repeat
        [ A,B ] = testrandomsymgenerator( symlist,dim1,dim2,prob );
        a=sparsesym(A);
        b=sparsesym(B);
        diff(i)=max(max(full(a-b)-(A-B)));
    end
    max(diff)
elseif any(test==3)
    fprintf('Test multiplication:\n')
    
   	for i=1:repeat
        [ A,B ] = testrandomsymgenerator( symlist,dim1,dim2,prob );
        a=sparsesym(A);
        b=sparsesym(B);
        diff(i)=max(max(full(a*b)-(A*B)));
    end
    max(diff)
end