function [ gates ] = Fast_Gauss_Jordan_Decomposition( matrix,plotting,name,circ_text )
%GAUSS_JORDAN_DECOMPOSITION decomposes a a unitary matrix into a set of
%controlled U2 rotations -> in order to implement on a quantum computer
mini=10^-12;
n=log2(size(matrix,1));

%Gates List
if nargin==1
    plotting=0;
    gates=struct('names',['U_n_' num2str(n) '_mat']);
    mode='opti';
elseif nargin==2
    gates=struct('names',['U_n_' num2str(n) '_mat']);
    mode='opti';
elseif nargin>2
    gates=struct('names',name);
end
if nargin<4
    circ_text='U_n';
end

matrix2=diag(ones(1,length(matrix)));

gates.matrix=matrix;
gates.size=n;
gates=Generate_Gate_Circuit(gates,circ_text,[1:n]);
gates.steps.index={};
gates.steps.gates={};
gates.steps.param={};
    
if n==1         
    gates.anc_size=0;
    gates.steps.index={[1]};
    gates.steps.gates={['U']};
    [a,b,d,t]=Unitary2Angles(matrix);
    gates.steps.param={[{sym('alpha'),a,sym('beta'),b,sym('delta'),d,sym('theta'),t}]};
else
    
    if plotting>=1
        FockPrint(matrix);
    end
    if n>2
        gates.anc_size=n-2;
        o=[n:-1:1 n+1:2*n-2];
    else
        gates.anc_size=0;
        o=[n:-1:1];
    end
    par1=sym('alpha');
    par2=sym('beta');
    par3=sym('delta');
    par4=sym('theta');
    par5=sym('phi');
    u_fun=@(alpha,beta,delta,theta)reshape([cos(theta/2)*exp(-(alpha*1i)/2-(beta*1i)/2+delta*1i) , sin(theta/2)*exp((alpha*1i)/2-(beta*1i)/2+delta*1i) , -sin(theta/2)*exp(-(alpha*1i)/2+(beta*1i)/2+delta*1i) , cos(theta/2)*exp((alpha*1i)/2 + (beta*1i)/2 + delta*1i) ],[2,2]);
    one=sparse(diag(ones(1,2^n)));
    c_bits=[];
    %Determine angles first, then reverse order and create gates
    angles_list=[];
    
    indexes2=[];
    bi=[];
    for i=1:2^n           %Columns
        if plotting>=1
          	fprintf([num2str(i) '------------------------------------------------------------------\n'])
        end
        
        %Find non zero elements in column i
        mat_i=matrix(:,i);
        ind=find(abs(mat_i)>mini); %Non trivial elements in column i
%         to_do=unique([i;ind(ind>=i)]);
%         bits=dec2bin(to_do-1,n)-'0';
%         intermed=bits(1,:);
%         last=bits(1,:);
%         for k=2:size(bits,1)
%             differ=bits(k,:)-bits(k-1,:);
%             for l=1:length(differ)
%                 if differ(l)~=0                
%                     last(l)=last(l)+differ(l);
%                     intermed=[intermed;last];
%                 end
%             end
%         end
    %% New faster approach    
        to_do=unique([i;ind(ind>=i)]);
        bits=dec2bin(to_do-1,n)-'0';  %Bits representation of non trivial bits
        bits2=dec2bin(i-1,n)-'0';     %Bit representation of diagonal element
        intermeds={};
        %Find the distance between bits 
        len_to_do=length(to_do);
        dist_mat=NaN(len_to_do); %Bit distance
        dist_diag=[];
        for k=1:len_to_do
            a=bits(k,:);
            for j=k+1:len_to_do
                b=sum(abs(a-bits(j,:)));
                dist_mat(k,j)=b;
                dist_mat(j,k)=b;
            end
            dist_diag(k)=sum(abs(a-bits2));
        end
        % Start with furthest (max(dist diag))
        % Find minimum distance of 
        [~,order]=sort(dist_diag,'descend');
        
        curr_intermed=0;
        for k=1:len_to_do-1 %Go through the elements from furthest to closest
            curr_bits=bits(order(k),:);
            %Find min distance from curr_bits
            [~,min_elem]=min(dist_mat(:,order(k)));
            goal_bits=bits(min_elem,:);
            differ=goal_bits-curr_bits;
            if k==1
                curr_intermed=curr_intermed+1;
                intermeds{curr_intermed}=[curr_bits];
            elseif isequal(intermeds(end,:),curr_bits)==0
                curr_intermed=curr_intermed+1;
                intermeds{curr_intermed}=[curr_bits];
            end
            
            if bi2de(curr_bits)>bi2de(goal_bits)
                index_order=1:length(differ);
            else
                index_order=length(differ):-1:1;
            end
            for l=index_order
                if differ(l)~=0                
                    curr_bits(l)=curr_bits(l)+differ(l);
                    intermeds{curr_intermed}=[intermeds{curr_intermed};curr_bits];
                end
            end
            dist_mat(order(k),:)=NaN;
            dist_mat(:,order(k))=NaN;
        end
        for k=1:length(intermeds)
            inter=intermeds{k};
            intermeds{k}=inter(end:-1:1,:);
        end
        
%Old approach
%         intermed=bits(1,:);
%         last=bits(1,:);
%         for k=2:size(bits,1)
%             differ=bits(k,:)-bits(k-1,:);
%             for l=1:length(differ)
%                 if differ(l)~=0                
%                     last(l)=last(l)+differ(l);
%                     intermed=[intermed;last];
%                 end
%             end
%         end
%         pause()
%         intermed
%         bits(order,:)
        if length(intermeds)==0
            intermeds{1}=bits;
        end
        
        for k=1:length(intermeds)
            intermed=intermeds{k};
            indexes=bin2dec(num2str(intermed,'%i'))+1;

            how_many=0;
            if length(indexes)>1 && i<=2^n-1
                indexes2=[indexes2;indexes];
                for j=length(indexes):-1:2   %Rows
                    index_now=[indexes(j),indexes(j-1)];
                    %find pu index
                    e=find(abs(intermed(j,:)-intermed(j-1,:)));
                    indi=intermed(j,:);
                    indi(e)=[];
                    b=bin2dec(num2str(indi,'%i'))+1;
                    pu_index_list=[b,e];
                    %Create a zero on the element
                    ind2zero=index_now(1);
                    a=matrix(sort(index_now),i);

                    if max(abs(a))>mini
                        u=sqrt(abs(a(1))^2+abs(a(2))^2);
                        U1=[a(2)/u,-a(1)/u;a(1)'/u,a(2)'/u];

                        if ind2zero~=min(index_now)
                            U1(1:2,:)=U1([2,1],:);
                            U1(2,:)=-U1(2,:);
                        end

                        [alpha,beta,delta,theta]=Unitary2Angles(U1);
                        [alpha2,beta2,delta2,theta2]=Unitary2Angles(U1');

                        %Check if angles are trivial
                        if max([alpha,beta,delta,theta])>mini
                            bi=[bi;ind2zero];
                            how_many=how_many+1;
                            d=dec2bin(indexes(j)-1,n)-'0';
                            c=d;
                            c(e)=NaN(1);   
                            c_bits=[c_bits;c];

                            u=[alpha,beta,delta,theta];
                            u2=[alpha2,beta2,delta2,theta2];
                            pu_matrix=one;
                            pu_matrix(sort(index_now),sort(index_now))=u_fun(u(1),u(2),u(3),u(4));
                            matrix=pu_matrix*matrix;
                            angles_list=[angles_list;u2(1),u2(2),u2(3),u2(4)]; %Conjugate of unitary angles

                            if plotting>=1
                                fprintf('cu ------------------------------------------------------------------\n')
                                %FockPrint(pu_matrix);
                                FockPrint(matrix);
                            end
                        end
                    end
                end
            end
            if how_many==0
                s=angle(matrix(i,i));
                if abs(s)>mini
                    %Rotate via controlled rotation
                    pu_matrix=one;
                    pu_matrix(i,i)=exp(-1i*s);
                    matrix=pu_matrix*matrix;
                    angles_list=[angles_list;NaN,NaN,NaN,s]; %Conjugate of unitary angles
                    bi=[bi;i];

                    d=dec2bin(i-1,n)-'0';
                    c_bits=[c_bits;d];

                    if plotting>=1
                        fprintf('cz ------------------------------------------------------------------\n')
                        %FockPrint(pu_matrix);
                        FockPrint(matrix);
                    end
                end
            end
        end
    end
    %check_matrix=matrix-diag(ones(2^n,1));

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
    %if plotting>=1
    %    [c_bits how_long how_long2 order]
    %end
    
    bi=dec2bin(bi-1,n)-'0';
    bi=bi(:,end:-1:1);
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

            if n>2            
                if isnan(angles_list(end-i+1,1))
                    %Controlled Z
                    if bi(end-i+1,o(curr_order(n)))==1
                        gates.steps.gates={gates.steps.gates{:},'C1PHASE1'};
                    else
                        gates.steps.gates={gates.steps.gates{:},'C1PHASE0'};
                    end

                    params=[{par5,angles_list(end-i+1,4)}];
                    gates.steps.param={gates.steps.param{:},params};
                else
                    gates.steps.gates={gates.steps.gates{:},'CU'};

                    params=[{par1,angles_list(end-i+1,1),par2,angles_list(end-i+1,2),par3,angles_list(end-i+1,3),par4,angles_list(end-i+1,4)}];
                    gates.steps.param={gates.steps.param{:},params};
                end
                gates.steps.index={gates.steps.index{:},[o(curr_order(n)),n*2-2]};
            else
                if isnan(angles_list(end-i+1,1))
                    %Controlled Z
                    if bi(end-i+1,o(curr_order(1)))==1
                        if bi(end-i+1,o(curr_order(2)))==1
                            gates.steps.gates={gates.steps.gates{:},'C1PHASE1'};
                        else
                            gates.steps.gates={gates.steps.gates{:},'C1PHASE0'};
                        end
                    else
                        if bi(end-i+1,o(curr_order(2)))==1
                            gates.steps.gates={gates.steps.gates{:},'C0PHASE1'};
                        else
                            gates.steps.gates={gates.steps.gates{:},'C0PHASE0'};
                        end
                    end

                    params=[{par5,angles_list(end-i+1,4)}];
                    gates.steps.param={gates.steps.param{:},params};
                else
                    if bi(end-i+1,o(curr_order(1)))==1
                        gates.steps.gates={gates.steps.gates{:},'CU'};
                    else
                        gates.steps.gates={gates.steps.gates{:},'C0U'};
                    end

                    params=[{par1,angles_list(end-i+1,1),par2,angles_list(end-i+1,2),par3,angles_list(end-i+1,3),par4,angles_list(end-i+1,4)}];
                    gates.steps.param={gates.steps.param{:},params};
                end
                
                if o(curr_order(n))==1
                    gates.steps.index={gates.steps.index{:},[o(curr_order(n)),2]};
                else
                    gates.steps.index={gates.steps.index{:},[o(curr_order(n)),1]};
                end
            end

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
            
           %Add Unitary rotations
           if n>2            
                if isnan(angles_list(end-i+1,1))
                    %Controlled Z
                    if bi(end-i+1,o(curr_order(n)))==1
                        gates.steps.gates={gates.steps.gates{:},'C1PHASE1'};
                    else
                        gates.steps.gates={gates.steps.gates{:},'C1PHASE0'};
                    end

                    params=[{par5,angles_list(end-i+1,4)}];
                    gates.steps.param={gates.steps.param{:},params};
                else
                    gates.steps.gates={gates.steps.gates{:},'CU'};

                    params=[{par1,angles_list(end-i+1,1),par2,angles_list(end-i+1,2),par3,angles_list(end-i+1,3),par4,angles_list(end-i+1,4)}];
                    gates.steps.param={gates.steps.param{:},params};
                end
                gates.steps.index={gates.steps.index{:},[o(curr_order(n)),n*2-2]};
            else
                if isnan(angles_list(end-i+1,1))
                    %Controlled Z
                    if bi(end-i+1,o(curr_order(1)))==1
                        if bi(end-i+1,o(curr_order(2)))==1
                            gates.steps.gates={gates.steps.gates{:},'C1PHASE1'};
                        else
                            gates.steps.gates={gates.steps.gates{:},'C1PHASE0'};
                        end
                    else
                        if bi(end-i+1,o(curr_order(2)))==1
                            gates.steps.gates={gates.steps.gates{:},'C0PHASE1'};
                        else
                            gates.steps.gates={gates.steps.gates{:},'C0PHASE0'};
                        end
                    end

                    params=[{par5,angles_list(end-i+1,4)}];
                    gates.steps.param={gates.steps.param{:},params};
                else
                    if bi(end-i+1,o(curr_order(1)))==1
                        gates.steps.gates={gates.steps.gates{:},'CU'};
                    else
                        gates.steps.gates={gates.steps.gates{:},'C0U'};
                    end

                    params=[{par1,angles_list(end-i+1,1),par2,angles_list(end-i+1,2),par3,angles_list(end-i+1,3),par4,angles_list(end-i+1,4)}];
                    gates.steps.param={gates.steps.param{:},params};
                end
                
                if o(curr_order(n))==1
                    gates.steps.index={gates.steps.index{:},[o(curr_order(n)),2]};
                else
                    gates.steps.index={gates.steps.index{:},[o(curr_order(n)),1]};
                end
            end
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
gates.step_num=length(gates.steps.index);
end

