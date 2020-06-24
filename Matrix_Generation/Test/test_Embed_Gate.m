%test_Embed_Gate.m
clc;
clear all;
close all;

p=pwd;
elem=strsplit(p,'\');
shortened=strjoin(elem((1:length(elem)-2)),'\');
addpath(genpath(shortened))


xxx=RandU(8);

N=10;

t1=tic();
for j=1:200
    list=1:N;
    for i=1:3
        a=ceil((N+1-i)*rand(1));
        indexes(i)=list(a);
        list(a)=[];
    end
    A=full(Embed_Gate(xxx,indexes,N));
end
t_old=toc(t1)

t1=tic();
for j=1:200
    list=1:N;
    for i=1:3
        a=ceil((N+1-i)*rand(1));
        indexes(i)=list(a);
        list(a)=[];
    end
    A_new=full(Embed_Gate_New(xxx,indexes,N));
end
t_new=toc(t1)

A_new=full(Embed_Gate_New(xxx,indexes,N));
A=full(Embed_Gate(xxx,indexes,N));
max(max(A-A_new))