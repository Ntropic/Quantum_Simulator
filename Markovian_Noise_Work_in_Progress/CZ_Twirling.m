function [ K,p ] = CZ_Twirling( phi,delta,E1 )
%CZ_TWIRLING 

%Pauli matrix
x=[0 1; 1 0];
y=[0 -1i; 1i 0];
z=[1 0; 0 -1];
i=[1 0; 0 1];

I=kron(i,i);
iZ=kron(i,z);
Zi=kron(z,i);
XX=kron(x,x);
YY=kron(y,y);
XY=kron(x,y);
YX=kron(y,x);
ZZ=kron(z,z);

p(1)=abs(1/4*(1+2*sqrt(1-E1)+exp(1i*delta)))^2;
p(2)=abs(1/4*(1-exp(1i*delta)))^2;
p(3)=abs(1/4*(1-exp(1i*delta)))^2;
p(4)=abs(1i/2*sqrt(E1)*sin(phi))^2;
p(5)=abs(1i/2*sqrt(E1)*sin(phi))^2;
p(6)=abs(1i/2*sqrt(E1)*cos(phi))^2;
p(7)=abs(1i/2*sqrt(E1)*cos(phi))^2;
p(8)=abs(1/4*(1-2*sqrt(1-E1)+exp(1i*delta)))^2;

K={iZ,Zi,XX,YY,XY,YX,ZZ};
end

