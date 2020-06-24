%fidelity_statistics.m
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

%Create random unitary matrices and determine their average fidelitieswith
%another random unitary matrix
s=4;
how_many=10^4;
bins=1000;

x=linspace(0,1-1/bins,bins)+1/bins/2;
for j=1:2^s-1
    fprintf([num2str(j+1) '/' num2str(2^s) '\n']);
    n=j+1;
    %Create random Unitary
    U=RandU(n);

    for i=1:how_many
        U2=RandU(n);

        F_av(i)=Average_Fidelity(U,U2);
    end

    A(j,:)=hist(F_av,x);
end
%Fidelity distribution can be described by 
image(x,2:2^s,A)
axis tight;
xlabel('size(U)')
ylabel('F')
zlabel('#F')
rotate3d on;

hold on
plot(1./((2:2^s)+1),2:2^s,'k')

figure()
x0=1/(m+1);
not_found=0;
i=0;
while not_found==0
    i=i+1;
    if x(i)>x0
        not_found=1;
    end
end
f=A(m,i:end);
xi=x(i:end);
plot(xi,f);
hold on;