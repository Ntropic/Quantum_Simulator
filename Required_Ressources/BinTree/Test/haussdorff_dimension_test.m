%haussdorff_dimension_test.m
clc;
clear all;
close all;

%Formula:
mini=1000000;
maxi=1000000;

dim=@(x,m) -log(sum(2.^(m-x).*sqrt(2.^(2*x-2)+1)))/log(2^(-m));
dim_adv=@(x,m) 1+log(sum(2.^(-x).*sqrt(2.^(2*x-2)+1)))/(m*log(2));
dim_opt=@(x,m) 1+log(sum(sqrt(2^(-2)+2.^(-2*x))))/(m*log(2));
for s=mini:maxi
    dim_opt(1:s,s)
end