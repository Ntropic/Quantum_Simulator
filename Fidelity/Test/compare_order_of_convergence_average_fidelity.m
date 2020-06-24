%Average Fidelity comparisons of Order
clc;
clear all;
close all;

p=pwd;
if any(strfind(p,'\'));
    elem=strsplit(p,'\');
else
    elem=strsplit(p,'/');
end
shortened=fullfile(elem{1:end-2});
addpath(genpath(shortened));

x=[0 1; 1 0];
y=[0 -1i; 1i 0];
z=[1 0; 0 -1];
I=diag([1 1]);

xx=kron(x,x);
yy=kron(y,y);
II=kron(I,I);

Ixx=kron(I,xx);
Iyy=kron(I,yy);
xxI=kron(xx,I);
yyI=kron(yy,I);

zII=kron(z,II);
IIz=kron(II,z);
IzI=kron(I,kron(z,I));

A=xxI+yyI;
B=Ixx+Iyy;
C=zII+IzI+IIz;

H2=A+B;
H3=A+B+C;


tn=pi/2./(1:10);

%H2
U2_exact=expm(-1i*H2*pi/2);
for i=1:10
    U2_approx=(expm(-1i*A*tn(i))*expm(-1i*B*tn(i)))^i;
    U2_sym_approx=(expm(-1i*A*tn(i)/2)*expm(-1i*B*tn(i))*expm(-1i*A*tn(i)/2))^i;
    
    F_ap(i)=Average_Fidelity(U2_exact,U2_approx);
    F_sym_ap(i)=Average_Fidelity(U2_exact,U2_sym_approx);
end
O2_ap=Order_Of_Convergence(1-F_ap);
O2_sym_ap=Order_Of_Convergence(1-F_sym_ap);
plot(O2_ap);
hold on;
plot(O2_sym_ap);
title('H2')

%H3
figure()
U3_exact=expm(-1i*H3*pi/2);
for i=1:10
    U3_approx=(expm(-1i*A*tn(i))*expm(-1i*B*tn(i))*expm(-1i*C*tn(i)))^i;
    U3_sym_approx=(expm(-1i*A*tn(i)/2)*expm(-1i*B*tn(i)/2)*expm(-1i*C*tn(i))*expm(-1i*B*tn(i)/2)*expm(-1i*A*tn(i)/2))^i;
    
    F3_ap(i)=Average_Fidelity(U3_exact,U3_approx);
    F3_sym_ap(i)=Average_Fidelity(U3_exact,U3_sym_approx);
end
O3_ap=Order_Of_Convergence(1-F3_ap);
O3_sym_ap=Order_Of_Convergence(1-F3_sym_ap);
plot(O3_ap);
hold on;
plot(O3_sym_ap);
title('H3')

%%Random Hamiltonian
figure()
for j=1:100
    Hr=RandH(8);
    Ur_exact=diag(ones(8,1));
    for i=1:10
        Ur_approx=expm(-1i*Hr*tn(i));
        Ur_sym_approx=expm(-1i*Hr*tn(i)^2);

        Fr_ap(i)=Average_Fidelity(Ur_exact,Ur_approx);
        Fr_sym_ap(i)=Average_Fidelity(Ur_exact,Ur_sym_approx);
    end
    Or_ap=Order_Of_Convergence(1-Fr_ap);
    Or_sym_ap=Order_Of_Convergence(1-Fr_sym_ap);
    plot(Or_ap);
    hold on;
    plot(Or_sym_ap);
end
title('H random')

%%Random Trotter
figure()
for j=1:100
    Ha=RandH(8);
    Hb=RandH(8);
    Hr=Ha+Hb;
    Ur_exact=expm(-1i*Hr*pi/2);
    for i=1:10
        Ur_approx=(expm(-1i*Ha*tn(i))*expm(-1i*Hb*tn(i)))^i;
        Ur_sym_approx=(expm(-1i*Ha*tn(i)/2)*expm(-1i*Hb*tn(i))*expm(-1i*Ha*tn(i)/2))^i;

        Fr_ap(i)=Average_Fidelity(Ur_exact,Ur_approx);
        Fr_sym_ap(i)=Average_Fidelity(Ur_exact,Ur_sym_approx);
    end
    Or_ap=Order_Of_Convergence(1-Fr_ap);
    Or_sym_ap=Order_Of_Convergence(1-Fr_sym_ap);
    plot(Or_ap);
    hold on;
    plot(Or_sym_ap);
end
title('H random')