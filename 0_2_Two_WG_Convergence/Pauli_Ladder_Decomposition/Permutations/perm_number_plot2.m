%perm_number_plot2.m
clc;
clear all;
close all;

% p=pwd;
% if any(strfind(p,'\'));
%     elem=strsplit(p,'\');
% else
%     elem=strsplit(p,'/');
% end
% clear p;
% shortened=fullfile(elem{1:end-3});
% addpath(genpath(shortened));

%Determine the number of decomposition gates of a Trotter step and the
%number of permutations, then compare the result for the use of fixed
%total photon number among the 2 waveguides participating in the exchange 
%and for the use of the symmetry n_max-i <-> 1+i

%All photon numbers
all_n_num_gates=@(n_max) n_max.*(n_max+1)/2;
num_perms=@(n) factorial(n);

%Single photon number
one_n_num_gates=@(n) n;
eff_num_perms=@(n) (n>1).*factorial(n)/2+(n==1); %using the symmetry of exchange between n_max-i <-> 1+i

%% Calculate
n_max=1:10;

%% Gate numbers
h=figure('Position', [10 10 800 400]); 
plot(n_max,all_n_num_gates(n_max),'b');
hold on;
plot(n_max,one_n_num_gates(n_max),'r');

xlabel('$N_\text{max}$');
ylabel('\# of Operators');
legend({'$n\leq N_\text{max}$','$n=N_\text{max}$'},'Location','northwest');
matlab2tikz('Permutation_Numbers/operators.tex','standalone',true,'parseStrings',false,'width','8cm','height','4cm')

%% Permutation number
h2=figure('Position', [10 10 800 400]); 
plot(n_max,num_perms(all_n_num_gates(n_max)),'b');
hold on;
plot(n_max,(eff_num_perms(one_n_num_gates(n_max))),'r');
plot(n_max,cumsum(eff_num_perms(one_n_num_gates(n_max))),'b--');
set(gca,'YScale', 'log');
ax2ytick=get(gca,'YTick');
set(gca,'YTick',ax2ytick);

xlabel('$N_\text{max}$');
ylabel('\# of Permutations');
legend({'Unoptimized ($n\leq N_\text{max}$)','Optimized ($n=N_\text{max}$)','Optimized ($n\leq N_\text{max}$)'},'Location','northwest');
matlab2tikz('Permutation_Numbers/permutations.tex','standalone',true,'parseStrings',false,'width','8cm','height','4cm')