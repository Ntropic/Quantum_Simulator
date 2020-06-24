%Compare_Photon_Numbers_Circ_Pauli_Gray_All_Interaction.m
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
n_m=1:7;
%n_m=1:3;

%Trotter Iteration number
steps=10;

%Waveguide number
N_m=2:4;
%N_m=2:6;

%Iteration time
t=pi/4;

name=['sym_compare_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_N_' num2str(max(N_m)) '_s_' num2str(steps) '_t_' num2str(t) '_pauli_gray'];
F_av_ind=zeros(length(N_m),length(n_m),steps);
F_av=zeros(length(N_m),length(n_m),steps);
for j=1:length(N_m)
    N=N_m(j);
    fprintf(['' num2str(N) ' Waveguides\n'])
    
    for i=1:length(n_m)
        fprintf(['  ->' num2str(i) ' Photons\n'])
        n=n_m(i);
        timer=tic();
        
        %% Exact decomposition
        s=ceil(log2(n+1));
      	Hij=Gray_Exchange_Hamiltonian_Particles(n);
        Hij_N=zeros(2^(N*s));
        for k=1:N-1
            i1=diag(ones(2^(s*(k-1)),1));
            i2=diag(ones(2^(s*(N-1-k)),1));
            Hij_N=Hij_N+kron(i2,kron(Hij,i1));
        end
        U_exact=expm(-1i*t*Hij_N);        %FockPrint(Hij)
        
        tree=bintree(s,1);
        gray=all_num_algorithm(tree);

        gray_ind_list=[];
        prefactor=[];
        for num1=2:n+1
            for num2=1:n+2-num1
                gray_num=gray([num1,num2]);
                graybin=dec2bin(gray_num-1,s)';
                graybin=graybin(:)';
                graybin=graybin-'0';
                
                gray_num2=gray([num1-1,num2+1]);
                graybin2=dec2bin(gray_num2-1,s)';
                graybin2=graybin2(:)';
                graybin2=graybin2-'0';
                gray_ind_list=[gray_ind_list; bin2dec(num2str(graybin))+1,bin2dec(num2str(graybin2))+1];
                prefactor=[prefactor;sqrt(num2)*sqrt(num1-1)];  %Strength of interaction 
            end
        end
        
        H=zeros(size(Hij));
        for l=1:steps
            fprintf(['        ' num2str(l) ' Trotter steps  - '])
  
            m=size(gray_ind_list,1);
            H=sparse(size(Hij,1),size(Hij,2));
            U_1=sparse(diag(ones(size(Hij,1),1)));
            U_2=sparse(diag(ones(size(Hij,1),1)));
            for o=1:m
                H_step=H;
                H_step(gray_ind_list(o,1),gray_ind_list(o,2))=prefactor(o);
                H_step(gray_ind_list(o,2),gray_ind_list(o,1))=prefactor(o);
                U_step{o}=expm(-1i*H_step*t/l/2);
            end
            for o=1:m
                U_1=U_step{o}*U_1;
                U_2=U_step{m-o+1}*U_2;
            end

            fprintf([' Elements Done  - '])

            %for the different waveguides
            U_step_N=diag(ones(2^(N*s),1));
            for o=1:N-1
                i1=diag(ones(2^(s*(o-1)),1));
                i2=diag(ones(2^(s*(N-1-o)),1));
                U_step_N=kron(i2,kron(U_1,i1))*U_step_N;
            end
            for o=N-1:-1:1
                i1=diag(ones(2^(s*(o-1)),1));
                i2=diag(ones(2^(s*(N-1-o)),1));
                U_step_N=kron(i2,kron(U_2,i1))*U_step_N;
            end
            fprintf([' Multiplication done \n'])

            %Test fidelity
            U_approx=U_step_N^l;
            
            %Find indexes
            indexes=Multi_Gray_Indexes(N,n);
            [F_av(j,i,l) F_av_ind(j,i,l)]=Average_Fidelity(U_approx,U_exact,indexes);
        end
    end
end


save([name '.mat'],'F_av_ind','F_av_ind','N_m');


%% Plot 1-Fidelity
h=figure('Position',[100,100,1400,600]);
m=size(F_av_ind,1);
for i=1:m
    ax{i}=subplot(1,m,i);
    err=reshape(1-F_av_ind(i,:,:),size(F_av_ind,2),size(F_av_ind,3));
    surf(err');
    if i==1
        xlabel('Number of photons')
        zlabel('1-Fidelity')
    end
    if i==m
        ylabel('Trotter steps')
    end
    title([num2str(N_m(i)) ' Waveguides']);
    view(-135,30)
    axis tight;
    drawnow;
end

%% Plot Convergence of 

% For Trotter steps
h_c_t=figure('Position',[100,100,1400,600]);
for i=1:m
    ax_h_t{i}=subplot(1,m,i);
    c=[];
    for j=2:size(F_av_ind,2)
        c(:,j-1)=Order_Of_Convergence(1-F_av_ind(i,j,:),1:steps);
    end
    surf(2:steps,2:max(n_m),c');
    if i==1
        xlabel('Trotter steps')
        zlabel('Convergence Order (# of Trotter steps)')
    end
    if i==m
        ylabel('Number of photons')
    end
    title([num2str(N_m(i)) ' Waveguides']);
    view(-135,30)
    axis tight;
    drawnow;
end

% For Photon numbers
h_c_p=figure('Position',[100,100,1400,600]);
for i=1:m
    ax_h_t{i}=subplot(1,m,i);
    c=[];
    for j=1:size(F_av_ind,3)
        c(j,:)=Order_Of_Convergence(1-F_av_ind(i,2:end,j),2:max(n_m));
    end
    surf(3:max(n_m),1:steps,c);
    if i==1
        xlabel('Number of photons')
        zlabel('Convergence Order (# of Photons)')
    end
    if i==m
        ylabel('Trotter steps')
    end
    title([num2str(N_m(i)) ' Waveguides']);
    view(-135,30)
    axis tight;
    drawnow;
end
