%No_Circ_Higher_Pauli_Gray_All_Interaction.m
%for trippling sheme

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

plotting=0;
%Photon number
n=5;
s=ceil(log2(n+1));

%Waveguide number
N=3;

%Trotter Iteration number
steps=10;
mode={'fock_gray',2^s-1,N,0,16};

%Iteration time
t=pi/4;

[seq,d_seq]=Higher_Trotter_Time_Steps(0);

Hij=Gray_Exchange_Hamiltonian_Particles(n);
Hij_N=zeros(2^(N*s));
for i=1:N-1
    i1=diag(ones(2^(s*(i-1)),1));
    i2=diag(ones(2^(s*(N-1-i)),1));
    Hij_N=Hij_N+kron(i2,kron(Hij,i1));
end
U_exact=expm(-1i*t*Hij_N);

s=ceil(log2(n+1));
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

H=sparse(size(Hij,1),size(Hij,2));

H_step=H;
H_step(gray_ind_list(1,1),gray_ind_list(1,2))=prefactor(1);
H_step(gray_ind_list(1,2),gray_ind_list(1,1))=prefactor(1);
H_s{1}=H_step;

for o=2:size(gray_ind_list,1)
    H_step=H;
    H_step(gray_ind_list(o,1),gray_ind_list(o,2))=prefactor(o);
    H_step(gray_ind_list(o,2),gray_ind_list(o,1))=prefactor(o);
    H_s{o}=H_step;
end
    
u1_cmp=zeros(steps,1);
for l=1:steps
    fprintf([num2str(l) ' Trotter steps  - '])
    
    m=size(gray_ind_list,1);
    %for the different waveguides
    for p=1:ceil(length(d_seq)/2)
        ti=d_seq(p)/l*t;
        U_step=diag(ones(2^(2*s),1));
        U_step2=U_step;
        for o=1:m
            U_step=expm(-1i*H_s{o}*ti/2)*U_step;
            U_step2=expm(-1i*H_s{m+1-o}*ti/2)*U_step2;
        end
        U_s{p}=U_step;
        U_s2{p}=U_step2;    
    end
    
    U_step_N=diag(ones(2^(N*s),1));
    for p=1:ceil(length(d_seq)/2)
        for i=1:N-1
            i1=diag(ones(2^(s*(i-1)),1));
            i2=diag(ones(2^(s*(N-1-i)),1));
            U_step_N=kron(i2,kron(U_s{p},i1))*U_step_N;
        end
        for i=N-1:-1:1
            i1=diag(ones(2^(s*(i-1)),1));
            i2=diag(ones(2^(s*(N-1-i)),1));
            U_step_N=kron(i2,kron(U_s2{p},i1))*U_step_N;
        end
    end
    for p=floor(length(d_seq)/2):-1:1
        for i=1:N-1
            i1=diag(ones(2^(s*(i-1)),1));
            i2=diag(ones(2^(s*(N-1-i)),1));
            U_step_N=kron(i2,kron(U_s{p},i1))*U_step_N;
        end
        for i=N-1:-1:1
            i1=diag(ones(2^(s*(i-1)),1));
            i2=diag(ones(2^(s*(N-1-i)),1));
            U_step_N=kron(i2,kron(U_s2{p},i1))*U_step_N;
        end
    end
    
    U_approx=U_step_N^l;
        
    fprintf([' Multiplication done \n'])

    %Test fidelity
    %Find indexes
    indexes=Multi_Gray_Indexes(N,n);
    [F_av(l) F_av_ind(l)]=Average_Fidelity(U_approx,U_exact,indexes);

    %Plotting
    if plotting==1
        ax1=subplot(1,2,1);
        Add_PColorMat(ax1,U_approx,mode);
        title('Approximation')

        ax2=subplot(1,2,2);
        Add_PColorMat(ax2,U_exact,mode);
        title('Exact Result')
        pause(0.2)
    end
end 

f1_cmp=1-F_av_ind;

figure()
plot(1:steps,f1_cmp,'r')
xlabel('Trotter steps')
ylabel('1-Average Fidelity')
legend({'\mathrm{O}(t^3)'})
axis tight;
%matlab2tikz('Average_Fidelity_Convergence1_50.tex')
drawnow;