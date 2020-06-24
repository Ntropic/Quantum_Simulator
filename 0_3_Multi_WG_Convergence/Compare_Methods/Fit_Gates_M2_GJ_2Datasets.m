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
x2=linspace(1,18,100);
%GJ
s_q_ind=s_q{3,1};
s_q_1=s_q_ind(1,:,1,1);
s_q_2=s_q_ind(1,:,1,2);

s=s_q_2+s_q_1;


%Fit 
[xData, yData] = prepareCurveData( x, s );
ft = fittype( 'a*(x+1)*(x+2)*(2*x+1)*1/6', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = 0.983317132027202;

% Fit model to data.
[fitresult,gof]=fit(xData,yData,ft,opts);
a=fitresult.a
gof.adjrsquare
da=confint(fitresult)
da2=confint(fitresult,0.68)

fun=@(x,a) a*(x+1).*(x+2).*(2*x+1)*1/6;
h=figure('Position',[100,100,800,600]);    
fill([x2,x2(end:-1:1)],[fun(x2,da(2)),fun(x2(end:-1:1),da(1))],0.9*[1 1 1],'EdgeColor','none')
hold on;
fill([x2,x2(end:-1:1)],[fun(x2,da2(2)),fun(x2(end:-1:1),da2(1))],0.8*[1 1 1],'EdgeColor','none')
plot(x2,fun(x2,a),'k');
plot(1:18,s,'ko');

axis([1 18 ylim]);
%set(gca, 'YScale', 'log')
xlabel('$N_\text{max}$')
ylabel('$n_\text{G}(N_\text{max})$')


matlab2tikz(['fit_gate_number_gj.tex'],'standalone',true,'parsestrings',false,'width','6cm','height','4cm')
