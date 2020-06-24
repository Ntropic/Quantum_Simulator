%perm_number_plot.m
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

%% All photon numbers (unoptimized)
h=figure('Position', [10 10 800 400]); 
plot(n_max,all_n_num_gates(n_max),'b');
hold on;
ax1=gca;
pos=get(ax1,'Position');
ax1x=get(ax1,'xlim');
ax1y=get(ax1,'ylim');
set(ax1,'box','off')
plot(1,0,'r');
axis([ax1x ax1y]);
legend({'\# of Operators','\# of Permutations'},'Location','northwest');
xlabel('$N_\text{max}$');

ax2=axes('Position',pos,'Color','none');
ax2.YAxisLocation='right';
plot(n_max,num_perms(all_n_num_gates(n_max)),'r');
set(ax2,'YAxisLocation','right');
set(ax2,'XAxisLocation','top');
set(ax2,'XTick',[]);
set(ax2,'XTickLabel',{});
set(ax2, 'YScale', 'log');
ax2ytick=get(ax2,'YTick');
set(ax2,'YTick',ax2ytick);
set(ax2,'Color','none');
set(ax2,'box','off')

ylabel(ax1,'\# of Operators');
ylabel(ax2,'\# of Permutations');

matlab2tikz('Permutation_Numbers/unoptimized_all_num.tex','standalone',true,'parseStrings',false,'width','8cm','height','4cm')

%% Single photon number (symmetry optimized)
h2=figure('Position', [10 10 800 400]); 
plot(n_max,one_n_num_gates(n_max),'b');
hold on;
ax1=gca;
pos=get(ax1,'Position');
ax1x=get(ax1,'xlim');
ax1y=get(ax1,'ylim');
set(ax1,'box','off')
plot(1,0,'r');
axis([ax1x ax1y]);
legend({'\# of Operators','\# of Permutations'},'Location','northwest');
xlabel('$N_\text{max}$');

ax2=axes('Position',pos,'Color','none');
ax2.YAxisLocation='right';
plot(n_max,eff_num_perms(one_n_num_gates(n_max)),'r');
set(ax2,'YAxisLocation','right');
set(ax2,'XAxisLocation','top');
set(ax2,'XTick',[]);
set(ax2,'XTickLabel',{});
set(ax2, 'YScale', 'log');
ax2ytick=get(ax2,'YTick');
set(ax2,'YTick',ax2ytick);
set(ax2,'Color','none');
set(ax2,'box','off')

ylabel(ax1,'\# of Operators');
ylabel(ax2,'\# of Permutations');

matlab2tikz('Permutation_Numbers/sym_optimized_single_num.tex','standalone',true,'parseStrings',false,'width','8cm','height','4cm')

%% Cummulative photon number (symmetry optimized)
h2=figure('Position', [10 10 800 400]); 
plot(n_max,all_n_num_gates(n_max),'b');
hold on;
ax1=gca;
pos=get(ax1,'Position');
ax1x=get(ax1,'xlim');
ax1y=get(ax1,'ylim');
set(ax1,'box','off')
plot(1,0,'r');
axis([ax1x ax1y]);
legend({'\# of Operators','\# of Permutations'},'Location','northwest');
xlabel('$N_\text{max}$');

ax2=axes('Position',pos,'Color','none');
ax2.YAxisLocation='right';
plot(n_max,cumsum(eff_num_perms(one_n_num_gates(n_max))),'r');
set(ax2,'YAxisLocation','right');
set(ax2,'XAxisLocation','top');
set(ax2,'XTick',[]);
set(ax2,'XTickLabel',{});
set(ax2, 'YScale', 'log');
ax2ytick=get(ax2,'YTick');
set(ax2,'YTick',ax2ytick);
set(ax2,'Color','none');
set(ax2,'box','off')

ylabel(ax1,'\# of Operators');
ylabel(ax2,'\# of Permutations');

matlab2tikz('Permutation_Numbers/sym_optimized_all_num.tex','standalone',true,'parseStrings',false,'width','8cm','height','4cm')