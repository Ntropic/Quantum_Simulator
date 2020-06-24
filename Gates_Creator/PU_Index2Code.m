function [ code ] = PU_Index2Code( n , index )
%PU_INDEX2CODE converts a PU Gates Index in a CompleteSet of PU Gates into
%the Code of it's conditions (eg. 010x would be conditions 0 and 1 on the
%bits 2 to 4, that control bit 1) here x will be replaced by NaN
%n is the number of qbits

%The complete PU_gates list consists of 2^(n-1)xn entries, single index
%leads to 2^(n-1) gates controlling the same qbit 
index=index(:);
red_index=mod(index-1,2^(n-1))+1;
con_bit=floor((index-1)/2^(n-1));

red_index_bin=de2bi(red_index-1,n-1)';
num_list=repmat([n:-1:1],length(index),1)';
for i=1:length(index)
    num_list(:,i)=num_list(:,i)+n*(i-1);
end
%Remove con_bit elements from matrix
lin_con_bit=(con_bit(:)+1)+n*((1:length(index))-1)';
num_list(lin_con_bit)=[];
%num_list=reshape(num_list,n-1,length(index))'
code=NaN(n,length(index));
code(num_list)=red_index_bin(:);
code=reshape(code,n,length(index))';
end

