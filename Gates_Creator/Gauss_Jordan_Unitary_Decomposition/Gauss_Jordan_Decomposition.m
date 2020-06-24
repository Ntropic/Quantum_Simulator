function [ gates,check_matrix ] = Gauss_Jordan_Decomposition( matrix,PU_gates,name,circ_text,plotting )
%GAUSS_JORDAN_DECOMPOSITION decomposes a a unitary matrix into a set of
%controlled U2 rotations -> in order to implement on a quantum computer
mini=10^-12;
n=log2(size(matrix,1));

%Gates List
if nargin==2
    gates=struct('names',['U_n_' num2str(n) '_mat']);
    mode='opti';
elseif nargin>=3
    gates=struct('names',name);
end
if nargin<4
    circ_text='U_n';
    plotting=0;
elseif nargin==4
    plotting=0;
end 
if length(plotting)==1
    plotting(2)=0;
end

gates.matrix=matrix;
gates.size=n;
gates=Generate_Gate_Circuit(gates,circ_text,[1:n]);
gates.steps.index={};
gates.steps.gates={};
gates.steps.param={};
    
if n==1
    gates.anc_size=0;
    gates.steps.index{[1]};
    gates.steps.gates{['U']};
    gates.steps.param{Unitary2Angles(matrix)};
    check_matrix=[1 0; 0 1];
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
                    pu_matrix=double(subs(PU_gates(pu_index_list(counter)).matrix,[sym('alpha'),sym('beta'),sym('delta'),sym('theta')],[alpha,beta,delta,theta]));
                    matrix=pu_matrix*matrix;
                    angles_list(counter,:)=[-beta,alpha,-delta,-theta]; %Conjugate of unitary angles

                    if plotting(1)>=1
                        fprintf(['\n<strong>Gate: </strong> ' PU_gates(pu_index_list(counter)).names '\n']);
                        FockPrint(matrix);
                        if plotting(2)==1
                            FockPrint(PU_gates(pu_index_list(counter)).matrix);
                        end
                    end

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
    check_matrix=matrix-diag(ones(2^n,1));

    %% Create the gates in an optimized fashion

    %How long do the CN code elements remain -> for ordering of conditions
    c_bits=c_bits(end:-1:1,:);

    for j=1:n %Row for row - qbit for qbit
        %Determine repetitions
        k=[1 find(diff(c_bits(:,j))')+1];
        p=[k numel(c_bits(:,j))]-[0 k];
        c=arrayfun(@(x) x:-1:1,p,'un',0);
        how_long(:,j)=cat(2,c{:})-1;

        k2=[find(diff(c_bits(:,j))')];
        p2=[k2 numel(c_bits(:,j))]-[0 k2];
        c2=arrayfun(@(x) x:-1:1,p2,'un',0);
        how_long2(:,j)=cat(2,c2{:})-1;
        how_long2(isnan(c_bits(:,j)),j)=-1;
    end
    [~,order]=sort(how_long2,2,'descend');
    if plotting(2)>=1
        [c_bits how_long how_long2 order]
    end

    %Adding Gates
    curr_order=zeros(1,n);
    for i=1:size(how_long,1)
        if max(curr_order)==0 
            %Create all new
            curr_order=order(i,:);
            
            index=[n+1,curr_order(1:2)];
            c_bit=c_bits(i,curr_order(1:2));
            if n>2
                gates=Add_Toffoli(gates,o(index),c_bit);
            end

            for j=1:n-3
                index=[n+1+j,curr_order(j+2),n+j];
                c_bit=[c_bits(i,curr_order(j+2)) 1];
                gates=Add_Toffoli(gates,o(index),c_bit);
            end

            gates.steps.gates={gates.steps.gates{:},'CU'};
            if n>2
                gates.steps.index={gates.steps.index{:},[o(curr_order(n)),n*2-2]};
            else
                if o(curr_order(n))==1
                    gates.steps.index={gates.steps.index{:},[o(curr_order(n)),2]};
                else
                    gates.steps.index={gates.steps.index{:},[o(curr_order(n)),1]};
                end
            end
            gates.steps.param={gates.steps.param{:},angles_list(counter-i,:)};

            %gates.steps.gates'
            %gates.steps.index{:}
        else
            %Which controls have to change? (New ordering?)
            if i<size(how_long,1)
                index_change=find(how_long(i,:)==0);
            else
                [k l]=find(diff(c_bits(i-1:i,:))); %l are the changing indexes
                index_change=l;
            end
            prev_order=curr_order;
            min_change=min(index_change);         %Check how deep the list has to be changed
            prev_ind=sort(prev_order(min_change:end));  
            [~,new_order]=sort(how_long2(i,prev_ind),2,'descend');
            curr_order=[curr_order(1:min_change-1),prev_ind(new_order)];
            
            if n>2
                %Undo the Toffoli Tree sufficiently and rebuild it
                if min_change>2
                    for j=n-1:-1:min_change
                        index=[n-1+j,prev_order(end+j-n),n-2+j];
                        c_bit=[c_bits(i-1,prev_order(j)) 1];
                        gates=Add_Toffoli(gates,o(index),c_bit);
                    end
                    for j=min_change:n-1
                        index=[n-1+j,curr_order(end+j-n),n-2+j];
                        c_bit=[c_bits(i,curr_order(j)) 1];
                        gates=Add_Toffoli(gates,o(index),c_bit);
                    end
                else
                    for j=n-1:-1:3
                        index=[n-1+j,prev_order(end+j-n),n-2+j];
                        c_bit=[c_bits(i-1,prev_order(end+j-n)) 1];
                        gates=Add_Toffoli(gates,o(index),c_bit);
                    end
                    index=[n+1,prev_order(1:2)];
                    c_bit=c_bits(i-1,prev_order(1:2));
                    gates=Add_Toffoli(gates,o(index),c_bit);

                    index=[n+1,curr_order(1:2)];
                    c_bit=c_bits(i,curr_order(1:2));
                    gates=Add_Toffoli(gates,o(index),c_bit);
                    for j=3:n-1
                        index=[n-1+j,curr_order(end+j-n),n-2+j];
                        c_bit=[c_bits(i,curr_order(j)) 1];
                        gates=Add_Toffoli(gates,o(index),c_bit);
                    end
                end
            end
            
            gates.steps.gates={gates.steps.gates{:},'CU'};
         	if n>2
                gates.steps.index={gates.steps.index{:},[o(curr_order(n)),n*2-2]};
            else
                if o(curr_order(n))==1
                    gates.steps.index={gates.steps.index{:},[o(curr_order(n)),2]};
                else
                    gates.steps.index={gates.steps.index{:},[o(curr_order(n)),1]};
                end
            end
            gates.steps.param={gates.steps.param{:},angles_list(counter-i,:)};

            %For plotting gates and corresponding indexes
            %for k=1:length(gates.steps.gates)
            %    fprintf([gates.steps.gates{k} ': ' num2str(gates.steps.index{k}) '\n'])
            %end
        end
    end
    if n>2
        for j=n-1:-1:3
            index=[n-1+j,curr_order(end+j-n),n-2+j];
            c_bit=[c_bits(i,curr_order(end+j-n)) 1];
            gates=Add_Toffoli(gates,o(index),c_bit);
        end
        index=[n+1,curr_order(1:2)];
        c_bit=c_bits(i,curr_order(1:2));
        gates=Add_Toffoli(gates,o(index),c_bit);
    end
end
end

