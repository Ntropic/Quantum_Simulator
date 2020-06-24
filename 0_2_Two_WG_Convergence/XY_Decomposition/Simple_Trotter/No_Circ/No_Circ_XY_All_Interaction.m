%No_Circ_XY_All_Interaction.m
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
%Load gates
load('../../../../Gates_Table/elem_gates.mat','-mat')
load('../../../../Gates_Table/comp_gates.mat','-mat')

n=6;
N=n-1;%Photon number

%Waveguide number
%Nwg=2;

%Trotter Iteration number
len=3;
mode={'fock_bin',2^n-1,1,0,16};

%Iteration time
t=pi/4;

x=[0 1; 1 0];
y=[ 0 -1i; 1i 0];
I=[0 0; 0 1];
O=[1 0; 0 0];
IO={O,I};

Hij=XY_Exchange_Hamiltonian(n);
U_exact=expm(-1i*t*Hij);

%mode={'fock_gray',n,N};
%PColorMat(Hij,mode);
%pause()

%Create Hamilton operators 
xx=kron(x,x);
yy=kron(y,y);
xxi=cell(n-1,1);
yyi=cell(n-1,1);
for i=1:n-1
    xxi{i}=Embed_Gate(xx,i+[0 1],n);
    yyi{i}=Embed_Gate(yy,i+[0 1],n);
end
Ha=sparse(2^n,2^n);
Hb=sparse(2^n,2^n);

index=2.^(0:n-1)+1;
Ni=N:-1:1;
in=1:N;
Jx=sqrt(in).*sqrt(Ni);
for i=1:2:n-1
    Ha=Ha+Jx(i)/2*(xxi{i}+yyi{i});
end
for i=2:2:n-1
    Hb=Hb+Jx(i)/2*(xxi{i}+yyi{i});
end

for l=1:len
    Ua=expm(-1i*t/l*Ha);
    Ub=expm(-1i*t/l*Hb);
        
    Ua2=expm(-1i*t/l*Ha/2);

    
    %Test if gates are correct
    U_approx1=(Ua*Ub)^l;
    U_approx2=Ua2*(Ub*Ua)^(l-1)*Ub*Ua2;

    %average_fidelity
    [F_av(l) F_av_ind(l)]=Average_Fidelity(U_approx1,U_exact,index);
    [F_av2(l) F_av_ind2(l)]=Average_Fidelity(U_approx2,U_exact,index);
    
    %Plotting
    %ax1=subplot(1,2,1);
    %Add_PColorMat(ax1,U_approx,mode);
    %title('Approximation')
    
    %ax2=subplot(1,2,2);
    %Add_PColorMat(ax2,U_exact,mode);
    %title('Exact Result')
end 

f1_cmp=1-F_av_ind;
f2_cmp=1-F_av_ind2;

figure()
plot(1:len,f1_cmp,'r')
hold on;
plot(1:len,f2_cmp,'r')
xlabel('Trotter steps')
ylabel('1-Average Fidelity')
legend({'\mathrm{O}(t^2)','\mathrm{O}(t^3)'})
axis tight;
%matlab2tikz('Average_Fidelity_Convergence1_50.tex')
drawnow;

figure()
n=2:len;
logn=log((n-1)./n);
u1_diff=log(f1_cmp(2:end))-log(f1_cmp(1:end-1));
u2_diff=log(f2_cmp(2:end))-log(f2_cmp(1:end-1));
semilogx(n,u1_diff./logn,'r')
hold on;
semilogx(n,u2_diff./logn,'r')
xlabel('Trotter steps')
ylabel('Convergence Order')
legend({'\mathrm{O}(t^2)','\mathrm{O}(t^3)'})
axis tight;
%matlab2tikz('Average_Fidelity_Convergence_order1_50.tex')
drawnow;