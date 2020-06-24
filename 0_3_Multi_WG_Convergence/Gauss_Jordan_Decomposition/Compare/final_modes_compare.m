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

%Convergence Plot
N=6;
M=4;


%Photon number
n_m=1:7;

%Trotter Iteration number
steps=10;

N_m=2:4;

r_h=0:3;

%Iteration time
t=pi/4;

%% Fidelities and gate numbers (gn)
%% ------------------------------------------------------------------------------------------------------------------------
exits='_gj';
trot_name=['mwg_trot_compare_fid_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_s_' num2str(steps) '_t_' num2str(t) '_N_' num2str(max(N_m))];
trot_name_gn=['mwg_trot_compare_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_s_' num2str(steps) '_t_' num2str(t) '_N_' num2str(max(N_m))];

sym_name=['mwg_sym_compare_fid_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_s_' num2str(steps) '_t_' num2str(t) '_N_' num2str(max(N_m))];
sym_name_gn=['mwg_sym_compare_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_s_' num2str(steps) '_t_' num2str(t) '_N_' num2str(max(N_m))];

fprintf('-> Check if simulation results are located in the designated folders:\n\n')

%% -------------------------------------------------------------------------------------------------------------------------
if exist(['../Simple_Trotter/' trot_name exits '.mat'])
    a=load(['../Simple_Trotter/' trot_name exits '.mat']);
    F_trot=a.F_av_ind;
    fprintf(' - Found Simple Trotter calculations\n')
else
    error('Not found (trott calc)')
end
if exist(['../Simple_Trotter/' trot_name_gn exits '.mat'])
    a=load(['../Simple_Trotter/' trot_name_gn exits '.mat']);
    s_q_trot=a.s_q;
    fprintf(' - Found Simple Trotter (gate number) calculations\n')
else
    error('Not found (trott gate numbers)')
end

%% -------------------------------------------------------------------------------------------------------------------------
if exist(['../Symmetric_Trotter/' sym_name exits '.mat'])
    a=load(['../Symmetric_Trotter/' sym_name exits '.mat']);
    F_sym=a.F_av_ind;
    fprintf(' - Found Symetric Trotter calculations\n')
else
    error('Not found (sym gcalc)');
end
if exist(['../Symmetric_Trotter/' sym_name_gn exits '.mat'])
    a=load(['../Symmetric_Trotter/' sym_name_gn exits '.mat']);
    s_q_sym=a.s_q;
    fprintf(' - Found Symetric Trotter (gate number) calculations\n')
else
    error('Not found (sym gate numbers)');
end

%% -------------------------------------------------------------------------------------------------------------------------
counter=0;
for r=r_h
    counter=counter+1;
    higher_name_r0=['mwg_higher_compare_fid_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_s_' num2str(steps) '_t_' num2str(t) '_r_' num2str(r) '_N_' num2str(max(N_m))];
    higher_name_r0_gn=['mwg_higher_compare_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_s_' num2str(steps) '_t_' num2str(t) '_r_' num2str(r) '_N_' num2str(max(N_m))];

    if exist(['../Higher_Order_Trotter/' higher_name_r0 exits '.mat'])
        a=load(['../Higher_Order_Trotter/' higher_name_r0 exits '.mat']);
        F_high_r0{counter}=a.F_av_ind;
        fprintf(' - Found Higher_Order_Trotter calculations\n')
    else
        error(['Not found (Higher_Order_Trotter calc) [r=' num2str(r) ']'])
    end
    if exist(['../Higher_Order_Trotter/' higher_name_r0_gn exits '.mat'])
        a=load(['../Higher_Order_Trotter/' higher_name_r0_gn exits '.mat']);
        s_q_high_r0{counter}=a.s_q;
        fprintf([' - Found Higher_Order_Trotter [r=' num2str(r) '] (gate number) calculations\n'])
    else
        error(['Not found (Higher_Order_Trotter gate numbers) [r=' num2str(r) ']'])
    end
end

F={};
s_q={};

F={F_trot,F_sym,F_high_r0{:}};
s_q={s_q_trot,s_q_sym,s_q_high_r0{:}};
len_sq=length(s_q);

    
fprintf('\n--> All files have been found, all variables extracted. <--\n')
fprintf('         --> !Proceeding with data analysis! <--\n\n')

%each element of f=F{i} contains a 3 dimensional tensor with the
%dimensions f(a,b,c) are 
%       a=number of modes N_m
%       b=number of particles n_m
%       c=number of trotter steps

%each element of s=s_q{i} contains a 4 dimensional tensor with the
%dimensions s(a,b,c,d) are 
%       a=number of modes N_m
%       b=number of particles n_m
%       c=number of trotter steps
%       d==1 -> single qubit gates
%       d==2 -> two qubit gates
%% Order of convergence
figure()

names={'\hat{T}','\hat{S}_1','\hat{S}_3^{(3)}','\hat{S}_3^{(5)}','\hat{S}_3^{(7)}','\hat{S}_3^{(9)}'};
stri={'-.r','r','--k','--b','--m','--c'};
for i=1:len_sq
    f=F{i};
    s=s_q{i};
    %   M, N,
    t=s(M-1,N ,:,1)+s(1,6,:,2);
    t=t(:);
    g(:)=f(M-1,N,:);
    g=1-g(:);
    h=Order_Of_Convergence(g);
    
    stringer=stri{i};
    plot(1:length(h),h,stringer);
    hold on;
    conv_order(i)=round(h(end));
end
ylabel('Order of Convergence');
xlabel('\# Trotter steps');
legend(names,'interpreter','latex')
matlab2tikz(['gj_convergence_N_' num2str(N) '_M_' num2str(M) '.tex'],'standalone',true,'parseStrings',false)

%% Errors
figure()
N=6;
names={'$\hat{T}$','$\hat{S}_1$','$\hat{S}_3^{(3)}$','$\hat{S}_3^{(5)}$','$\hat{S}_3^{(7)}$','$\hat{S}_3^{(9)}$'};
stri={'-.r','r','--k','--b','--m','--c'};
stri2={':r',':r',':k',':b',':m',':c'};
for i=1:len_sq
    f=F{i};
    s=s_q{i};
    %   2M, N,
    t=s(M-1 ,N ,:,1)+s(1,6,:,2);
    t=t(:);
    g(:)=f(M-1,N,:);
    g=1-g(:);
    
    %Extrapolate
    
    
    stringer=stri{i};
    semilogy(t,g,stringer);
    hold on;
end
x_lim=xlim();
y_lim=ylim();
for i=1:len_sq
    f=F{i};
    s=s_q{i};
    %   2M, N,
    t=s(M-1 ,N ,:,1)+s(1,6,:,2);
    t=t(:);
    g(:)=f(M-1,N,:);
    g=1-g(:);
    %Extrapolate
    if t(end)<x_lim(2)*0.95
        diff_t=diff(t);
        diff_t=diff_t(end);
        how_many=round((x_lim(2)-t(end))/diff_t);
        if how_many<1
            how_many=1;
        end
        n=length(t);
        
        t2=[t(end) t(end)+diff_t*(1:how_many)];
        
        eps=g(end)*n^conv_order(i);
        g2=eps./(n+(1:how_many)).^(conv_order(i));
        g2=[g(end) g2];
    end

    stringer=stri2{i};
    semilogy(t2,g2,stringer);
    hold on;
end
xlim(x_lim);
ylim(y_lim);
ylabel('$1-F$');
xlabel('\# Gates');
legend(names,'interpreter','latex')
matlab2tikz(['gj_errors_N_' num2str(N) '_M_' num2str(M) '.tex'],'standalone',true,'parseStrings',false)
