%test_higher_trotter_time_steps.m
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

[seq3,d_seq3]=Higher_Trotter_Time_Steps(0);
[seq5,d_seq5]=Higher_Trotter_Time_Steps(1);
plot(0:length(seq3)-1,seq3,'r');
hold on;
plot(0:length(seq3)-1,seq3,'ro');
plot(0:length(seq5)-1,seq5,'b');
plot(0:length(seq5)-1,seq5,'bo');
axis tight;
x_lim=xlim(gca);
y_lim=ylim(gca);
plot(xlim,[0 0],'k:')
plot(xlim,[1 1],'k:')
ylabel('t')
xlabel('steps')
hold off;

%name='first_trotter_it';
%[filename]=Plot2Tikz(name);
%pause();

for i=1:4
    [seq3,d_seq3]=Higher_Trotter_Time_Steps(zeros(1,i));
    [seq5,d_seq5]=Higher_Trotter_Time_Steps(ones(1,i));
    plot(0:length(seq3)-1,seq3,'r');
    hold on;
    plot(0:length(seq5)-1,seq5,'b');
    axis tight;
    x_lim=xlim(gca);
    y_lim=ylim(gca);
    plot(xlim,[0 0],'k:')
    plot(xlim,[1 1],'k:')
    ylabel('t')
    xlabel('steps')
    a=gca;
    a.YTick=[0,1];
    hold off;
    pause();
    name=['trotter_it_' num2str(i)];
    [filename]=Plot2Tikz(name);
end
