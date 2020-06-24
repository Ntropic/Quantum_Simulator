function [ alpha,beta,delta,theta ] = Unitary2Angles( U )
%UNITARY2ANGLES creates a separation of a unitary (2x2) matrix into angles 
%alpha,beta, delta, theta with the property
%U=[cos(theta/2)*exp(1i*(delta-alpha/2-beta/2)),-sin(theta/2)*exp(1i*(delta-alpha/2+beta/2))]
%  [sin(theta/2)*exp(1i*(delta+alpha/2-beta/2)), cos(theta/2)*exp(1i*(delta+alpha/2+beta/2))];

%Check if matrix is 2x2
if all(size(U)~=[2,2])
    error('Unitary is of wrong size. Must be 2x2 matrix.')
end
%Check if U is unitary
if abs(U(:,1)'*U(:,2))>5e-14
    error('Matrix is not unitary.')
end


%First separate unitary into:
% U=[a                  ,-b                 ]
%   [e^(i*2*delta)*conj(b),e^(i*2*delta)*conj(a)]
a=U(1,1);
b=-U(1,2);
delta=angle(U(2,1)/conj(b))/2;
if isnan(delta)
    delta=0;
end
%U2=[a ,-b ;exp(1i*2*delta)*conj(b),exp(1i*2*delta)*conj(a)]               %For testing the expansion

%Second separate the unitary into:
% U=e^(i*delta)[e^(i*phi1)*cos(theta) ,-e^(i*phi2)*sin(theta)]
%              [e1(-i*phi2)*sin(theta),e^(-i*phi1)*cos(theta)]
phi1=mod(atan(imag(a)/real(a))+[0 pi]-delta+2*pi,2*pi);
if isnan(phi1)
    phi1=[0 0];
end
phi2=mod(atan(imag(b)/real(b))+[0 pi]-delta+2*pi,2*pi);
if isnan(phi2)
    phi2=[0 0];
end
maxi=5;
ind=[3,3];
for i=1:2
    for j=1:2
        ct=real(a)*cos(delta+phi1(i))+imag(a)*sin(delta+phi1(i));
        st=real(b)*cos(delta+phi2(j))+imag(b)*sin(delta+phi2(j));
        theta2=mod((angle(ct+1i*st)+2*pi)*2,2*pi);
        
        Uij=exp(1i*delta)*[cos(theta2/2)*exp(1i*phi1(i)), -sin(theta2/2)*exp(1i*phi2(j)) ; sin(theta2/2)*exp(-1i*phi2(j)) , cos(theta2/2)*exp(-1i*phi1(i))];
        if max(abs(U(:)-Uij(:)))<maxi
            maxi=max(abs(U(:)-Uij(:)));
            ind=[i,j];
        end
    end
end

ct=real(a)*cos(delta+phi1(ind(1)))+imag(a)*sin(delta+phi1(ind(1)));
st=real(b)*cos(delta+phi2(ind(2)))+imag(b)*sin(delta+phi2(ind(2)));
phi1=phi1(ind(1));
phi2=phi2(ind(2));
theta=mod((angle(ct+1i*st)+2*pi)*2,2*pi);

%Third represent via:
% U=[cos(theta/2)*exp(1i*(delta-alpha/2-beta/2)),-sin(theta/2)*exp(1i*(delta-alpha/2+beta/2))]
%   [sin(theta/2)*exp(1i*(delta+alpha/2-beta/2)), cos(theta/2)*exp(1i*(delta+alpha/2+beta/2))];
%For which we substitute:
% phi1=-alpha/2-beta/2
% phi2=-alpha/2+beta/2

alpha=-(phi1+phi2);
beta=(phi2-phi1);
end