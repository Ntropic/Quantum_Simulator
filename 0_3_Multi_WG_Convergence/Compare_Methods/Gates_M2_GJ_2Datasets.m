%Compare_All_GateNumbers2Datasets.m
clear all;
close all;
clc;

p=pwd;
if any(strfind(p,'\'));
    elem=strsplit(p,'\');
else
    elem=strsplit(p,'/');
end
shortened=fullfile(elem{1:end-3});
addpath(genpath(shortened));

    %Photon number
    n_m=1:18;%1:7;

    %Trotter Iteration number
    steps=1;%0;

    N_m=2;%2:4;

    r_h=0:3;

    %Iteration time
    t=pi/4;

    %% Fidelities and gate numbers (gn)
    %% ------------------------------------------------------------------------------------------------------------------------
    exits={'_pauli_ladder','_pauli_gray','_gj'};%,'_gj2'};
    exit_Folders={'Pauli_Ladder_Decomposition','Pauli_Gray_Decomposition','Gauss_Jordan_Decomposition'};%,'Gauss_Jordan_Decomposition'};

    trot_name=['mwg_trot_compare_fid_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_s_' num2str(steps) '_t_' num2str(t) '_N_' num2str(max(N_m))];
    trot_name_gn=['mwg_trot_compare_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_s_' num2str(steps) '_t_' num2str(t) '_N_' num2str(max(N_m))];

    fprintf('-> Check if simulation results are located in the designated folders:\n\n')

    %% -------------------------------------------------------------------------------------------------------------------------
    for i=1:length(exits) %Only search for gate numbers
        fprintf(['- Searching for calculations of type ' exit_Folders{i} '\n'])

%         if exist(['../' exit_Folders{i} '/Simple_Trotter/' trot_name exits{i} '.mat'])
%             a=load(['../' exit_Folders{i} '/Simple_Trotter/' trot_name exits{i} '.mat']);
%             F_trot=a.F_av_ind;
%             fprintf(' - Found Simple Trotter calculations\n')
%         else
%             error('Not found (trott fid)')
%         end
        if exist(['../' exit_Folders{i} '/Simple_Trotter/' trot_name_gn exits{i} '.mat'])
            a=load(['../' exit_Folders{i} '/Simple_Trotter/' trot_name_gn exits{i} '.mat']);
            s_q_trot=a.s_q;
            fprintf(' - Found Simple Trotter (gate number) calculations\n')
        else
            error('Not found (trott gate numbers)')
        end
        
        if i==1
%            F={};
            s_q={};
        end
%        F={F{:},F_trot};
        s_q={s_q{:},s_q_trot};
        if i==1
            len_sq=length(s_q);
        end
    end

%F=reshape(F,len_sq,length(exits));
%F=F';
s_q=reshape(s_q,len_sq,length(exits));
s_q=s_q';


fprintf('\n--> All files have been found, all variables extracted. <--\n')
fprintf('         --> !Proceeding with data analysis! <--\n\n')

%% Analyze Data---------------------------------------------------------------------------------------------------------
x=1:18;
%GJ
s_q_ind=s_q{3,1};
s_q_1=s_q_ind(1,:,1,1);
s_q_2=s_q_ind(1,:,1,2);

h=figure('Position',[100,100,800,600]);    
fill([x x(end:-1:1)],[s_q_2+s_q_1 zeros(size(x))],[1 1 0]);
hold on;
fill([x x(end:-1:1)],[s_q_1 zeros(size(x))],[0 0 1]);

%set(gca, 'YScale', 'log')
xlabel('\# of photons')
ylabel('\# of gates')


s_q_ind=s_q{1,1};
s_q_1l=s_q_ind(1,:,1,1);
s_q_2l=s_q_ind(1,:,1,2);
s_q_2l=s_q_2l+s_q_1l;

plot(s_q_2l,'r')

s_q_ind=s_q{2,1};
s_q_1g=s_q_ind(1,:,1,1);
s_q_2g=s_q_ind(1,:,1,2);
s_q_2g=s_q_2g+s_q_1g;
ax=plot(s_q_2g);
ax.Color=[1 0.5 0];

legend({'Single qubit','Two qubits','Ladder','Gray'},'interpreter','latex','location','northwest');
%matlab2tikz(['gate_number_gj.tex'],'standalone',true,'parsestrings',false,'width','5cm','height','4cm')

%How many Trotter steps are necessary to get to the same amount of gates
a=(s_q_2+s_q_1)./s_q_2l; %Ladder
b=(s_q_2+s_q_1)./s_q_2g; %Gray

%% Fit data
fprintf('Ladder --------------------------------------\n')
h2=figure('Position',[100,100,800,600]);    
[xData,yData]=prepareCurveData(x,a);
ft = fittype( 'poly1' );
[fitresult, gof] = fit( xData, yData, ft );
f= fit( xData, yData, ft );

p1=fitresult.p1
p2=fitresult.p2
conf=confint(fitresult)
c_p1=conf(:,1);
c_p2=conf(:,2);
conf=confint(fitresult,0.68)
c2_p1=conf(:,1);
c2_p2=conf(:,2);
gof

fun=@(y,p_1,p_2) p_1*y+p_2;
% Plot fit with data.
fill([x,x(end:-1:1)],[fun(x,c_p1(1),c_p2(1)),fun(x(end:-1:1),c_p1(2),c_p2(2))],0.9*[1 1 1],'EdgeColor','none')
hold on;

fprintf('Gray ----------------------------------------\n')
[xData,yData]=prepareCurveData(x,b);
ft = fittype( 'poly1' );
[fitresult2, gof2] = fit( xData, yData, ft );
f2= fit( xData, yData, ft );

p12=fitresult2.p1
p22=fitresult2.p2
conf2=confint(fitresult2)
c_p12=conf2(:,1);
c_p22=conf2(:,2);
conf2=confint(fitresult2,0.68)
c2_p12=conf2(:,1);
c2_p22=conf2(:,2);
gof2

% Plot fit with data.
fill([x,x(end:-1:1)],[fun(x,c_p12(1),c_p22(1)),fun(x(end:-1:1),c_p12(2),c_p22(2))],0.9*[1 1 1],'EdgeColor','none')
fill([x,x(end:-1:1)],[fun(x,c2_p1(1),c2_p2(1)),fun(x(end:-1:1),c2_p1(2),c2_p2(2))],0.8*[1 1 1],'EdgeColor','none')
fill([x,x(end:-1:1)],[fun(x,c2_p12(1),c2_p22(1)),fun(x(end:-1:1),c2_p12(2),c2_p22(2))],0.8*[1 1 1],'EdgeColor','none')
ax=plot(x,f(x),'r');%fun(x,p1,p2),'r')
plot(a,'ro')
ax2=plot(x,f2(x));%fun(x,p12,p22));
ax2.Color=[1 0.5 0];
ax3=plot(x,b,'o');
ax3.Color=[1 0.5 0];
xlabel('\# of Trotter steps')
ylabel('ratio of \# of gates')
legend([ax ax2],{'Ladder','Gray'})

y_lim=ylim();
axis([x(1) x(end) y_lim ]);
%matlab2tikz(['ratio_gate_number_gj.tex'],'standalone',true,'parsestrings',false,'width','5cm','height','4cm')