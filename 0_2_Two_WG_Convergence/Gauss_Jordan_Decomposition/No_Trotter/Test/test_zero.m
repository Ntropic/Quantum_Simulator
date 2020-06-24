%test.m
clc;
clear all;
close all;

a=RandU(2);
a=1/sqrt(2)*[1 0;1 0 ]

index_now=[1,2];
u_fun=@(alpha,beta,delta,theta)reshape([cos(theta/2)*exp(-(alpha*1i)/2-(beta*1i)/2+delta*1i) , sin(theta/2)*exp((alpha*1i)/2-(beta*1i)/2+delta*1i) , -sin(theta/2)*exp(-(alpha*1i)/2+(beta*1i)/2+delta*1i) , cos(theta/2)*exp((alpha*1i)/2 + (beta*1i)/2 + delta*1i) ],[2,2]);

%--------------------------------------------------------------
ind2zero=1;
u=sqrt(abs(a(1))^2+abs(a(2))^2);
U1=[a(2)/u,-a(1)/u;a(1)'/u,a(2)'/u];


if ind2zero~=min(index_now)
    U1(1:2,:)=U1([2,1],:);
    U1(2,:)=-U1(2,:);
end

[alpha,beta,delta,theta]=Unitary2Angles(U1);
u=[alpha,beta,delta,theta];

U1=u_fun(u(1),u(2),u(3),u(4))*a

%--------------------------------------------------------------
ind2zero=2;
u=sqrt(abs(a(1))^2+abs(a(2))^2);
U1=[a(2)/u,-a(1)/u;a(1)'/u,a(2)'/u];


if ind2zero~=min(index_now)
    U1(1:2,:)=U1([2,1],:);
    U1(2,:)=-U1(2,:);
end

[alpha,beta,delta,theta]=Unitary2Angles(U1);
u=[alpha,beta,delta,theta];

U2=u_fun(u(1),u(2),u(3),u(4))*a