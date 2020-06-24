function [ mat ] = Plot_Gauss_Jordan_Decomposition( matrix,plotting )
%PLOT_GAUSS_JORDAN_DECOMPOSITION decomposes a a unitary matrix into a set 
%of controlled U2 rotations -> in order to implement on a quantum computer
% shows the evolution of the state in plot form
% plotting    - 0=matrix steps, 1=cummulative matrix steps
mini=10^-12;
n=log2(size(matrix,1));

if n==1
    mat{1}=matrix;
    mat{2}=diag(ones(1,2^n));
else
    %Prepare variables for gates

    %Create path for making elements zero
    index_mat=NumberPaths(n);
    index_comp_mat=zeros(size(index_mat,1)-1,size(index_mat,2)-1,2);
    index_comp_diff=zeros(size(index_mat)-[1,1]);

    index_comp_mat(:,:,1)=index_mat(2:end,1:end-1);
    index_comp_mat(:,:,2)=index_mat(1:end-1,1:end-1);
    index_comp_mat(:,:,1)=index_comp_mat(:,:,1)-diag((2:2^n-1),1);

    index_comp_diff(:,:)=log2(abs(index_comp_mat(:,:,2)-index_comp_mat(:,:,1)))+1;

    index_comp_mat2=reshape(index_comp_mat(:,:,2),(2^n-1)^2,1);
    index_i=dec2bin(max(index_comp_mat2(:)-1,0));
    for i=1:n
        index_i(index_comp_diff(:)==i,n+1-i)='x';
    end
    index_i(index_comp_diff(:)==-Inf,1)='x';
    index_i=reshape(index_i',1,(2^n-1)^2*n);
    index_i(index_i=='x')=[];
    index_i=reshape(index_i,n-1,length(index_i)/(n-1))';
    index_i=reshape(bin2dec(index_i)+1,2^n-1,2^n-1);

    pu_index_mat=index_i+(index_comp_diff-1)*2^(n-1);   


    %Which CN_rotation to use 
    pu_index_list=[];
    for i=1:2^n-1
        pu_index_list=[pu_index_list;pu_index_mat(end:-1:i,i)]; %From bottom up the list
    end
    how_long=zeros(length(pu_index_list),n);
    c_bits=PU_Index2Code(n,pu_index_list(:));

    if n>2
        gates.anc_size=n-2;
        o=[n:-1:1 n+1:2*n-2];
    else
        gates.anc_size=0;
        o=[n:-1:1];
    end

    %Determine angles first, then reverse order and create gates
    angles_list=zeros(length(pu_index_list),4);

    counter=1;
    pu_gate=diag(ones(1,2^n));
    for i=1:2^n-1           %Columns
        for j=2^n:-1:i+1    %Rows
            %Determine order 
            indexes=[index_mat(j,i),index_mat(j-1,i)];

            %Create a zero on the element
            ind2zero=indexes(1);
            a=matrix(sort(indexes),i);

            if max(abs(a))>mini
                u=sqrt(abs(a(1))^2+abs(a(2))^2);
                U1=[a(2)/u,-a(1)/u;a(1)'/u,a(2)'/u];

                if i==2^n-1
                    a2=matrix([2^n-1,2^n],[2^n-1,2^n]);
                    u=sqrt(abs(a(1))^2+abs(a(2))^2);
                    U1=[a(2)/u,-a(1)/u;a(1)'/u,a(2)'/u];
                    U1(1,:)=U1(1,:)*u/(a2(2,2)*a2(1,1)-a2(1,2)*a2(2,1));
                end

                if ind2zero~=min(indexes)
                    U1(1:2,:)=U1([2,1],:);
                    U1(2,:)=-U1(2,:);
                end
                
                [alpha,beta,delta,theta]=Unitary2Angles(U1);

                %Check if angles are trivial
                if max([alpha,beta,delta,theta])>mini
                    pu_matrix=pu_gate;
                    pu_matrix(sort(indexes),sort(indexes))=U1;
                    mat{counter}=pu_matrix';
                    matrix=pu_matrix*matrix;
                    counter=counter+1;
                else
                    pu_index_list(counter)=[];
                    angles_list(counter,:)=[];
                    c_bits(counter,:)=[];
                    how_long(counter,:)=[];
                end
            else
                pu_index_list(counter)=[];
                angles_list(counter,:)=[];
                c_bits(counter,:)=[];
                how_long(counter,:)=[];
            end
        end
    end
    mat{counter}=pu_gate;
end

%% Output 
mat=mat(end:-1:1);
if plotting==1    %Return cummulative matrix steps
    for i=2:length(mat)
        mat{i}=mat{i}*mat{i-1};
    end
end
end