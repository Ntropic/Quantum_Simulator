function [ H ] = Single_Particle_Exchange_Hamiltonian(n,p)
%Exchange_Hamiltonian Create n qubits exchange matrix for N modes with
%exchange between modes ij and ij+1
%   - log2(n) quibts per mode for storing particle number (via gray code)
%   - 2*log2(n) qubits in total for interaction
%   - p=[n0,n1_start,n+-] total number of particles n0 and starting
%   particles in mode1 and particles are added to or reduced from mode 1
%   (p(3)=+-1)

if length(p)~=3
    error('Length of p must be 3.')
end
if max(p(2))>p(1)
    error('Total number of particles p(1) must be larger than mode 1s.')
end
n1=p(2);
n0=p(1);
n_sign=sign(p(3));
if n1+n_sign>n0 | n1+n_sign<0
    error('Mode particles and its addition has to lie within 0 to n0')
end

H=zeros(2^(2*n),2^(2*n));


%Find correct elements within gray permutation
num1=n1+1;
num2=n0-n1+1;

bin_num=[num1,num2];
bin=dec2bin(bin_num-1,n)';
bin=bin(:)';
index=bin2dec(bin)+1;

gray_num2=[num1+n_sign,num2-n_sign];
bin2=dec2bin(gray_num2-1,n)';
bin2=bin2(:)';
index2=bin2dec(bin2)+1;
H(index,index2)=sqrt(num2)*sqrt(num1-1);

end