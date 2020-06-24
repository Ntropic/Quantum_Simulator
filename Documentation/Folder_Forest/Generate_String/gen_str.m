%gen_str.m
clc;
clear all;
close all;

d=dir;

str='';

for i=3:length(d)
    curr_name=d(i).name;
    str=[str '[\\path{' curr_name '}]\n'];
end
fprintf(str)