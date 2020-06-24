%No_Circ_Pauli_Gray_All_Interaction.m
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
n=4;
s=ceil(log2(n+1));

%Waveguide number
N=2;

%Trotter Iteration number
steps=10;
mode={'fock_gray',2^s-1,N,0,16};

%Iteration time
t=pi/4;


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

u1_cmp=zeros(steps,1);
for l=1:steps
    fprintf([num2str(l) ' Trotter steps  - '])
    
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
    for i=1:N-1
        i1=diag(ones(2^(s*(i-1)),1));
        i2=diag(ones(2^(s*(N-1-i)),1));
        U_step_N=kron(i2,kron(U_1,i1))*U_step_N;
    end
    for i=N-1:-1:1
        i1=diag(ones(2^(s*(i-1)),1));
        i2=diag(ones(2^(s*(N-1-i)),1));
        U_step_N=kron(i2,kron(U_2,i1))*U_step_N;
    end
    fprintf([' Multiplication done \n'])
    
    %Test fidelity
    U_approx=U_step_N^l;
    
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