function [ init main connector outit gate_names ] = Split_Strang_Pauli_Gray_Gates( n,name,circ_string )
%SPLIT_STRANG_PAULI_GRAY_GATES (Strang Splitting - Symmetric) Creates the 
%gates for Pauli Gray encoding and 2 mode interactions for n photons but 
%splits the gate into a init gate that creates the conditions on the 
%ancilla qubits with Strang Splitting the main gate and an outit, 
%that gets rid of the conditions the gate can be made symmetric 

if nargin==3
    circ_string_in= [circ_string '^{(1),(ie)}(\phi)'];
    circ_string_out=[circ_string '^{(1),(fe)}(\phi)'];
    circ_string_con=[circ_string '^{(1),(e)}(\phi)'];
    circ_string=    [circ_string '^{(1),(u)}(\phi)'];
    name_in=[name '_in'];
    name_out=[name '_out'];
    name_con=[name '_con'];
    name=[name '_cen'];
elseif nargin==2
    circ_string_in= ['S_{N\leq ' num2str(n) '}^{G,(ie)}(\phi)'];
    circ_string_out=['S_{N\leq ' num2str(n) '}^{G,(fe)}(\phi)'];
    circ_string_con=['S_{N\leq ' num2str(n) '}^{G,(e)}(\phi)'];
    circ_string=    ['S_{N\leq ' num2str(n) '}^{G,(u)}(\phi)'];
    name_in=[name '_in'];
    name_out=[name '_out'];
    name_con=[name '_con'];
    name=[name '_cen'];
elseif nargin==1
    circ_string_in= ['S_{N\leq ' num2str(n) '}^{G,(ie)}(\phi)'];
    circ_string_out=['S_{N\leq ' num2str(n) '}^{G,(fe)}(\phi)'];
    circ_string_con=['S_{N\leq ' num2str(n) '}^{G,(e)}(\phi)'];
    circ_string=    ['S_{N\leq ' num2str(n) '}^{G,(u)}(\phi)'];
    name_in=['Gray_n_' num2str(n) '_sym_in'];
    name_out=['Gray_n_' num2str(n) '_sym_out'];
    name_con=['Gray_n_' num2str(n) '_sym_con'];
    name=['Gray_n_' num2str(n) '_sym_cen'];
end
gate_names={name_in,name_out,name,name_con};
%% Prepare the lists of exchanges -> see comments for explanation of the meanings of each variable 
s=ceil(log2(n+1));
tree=bintree(s,1);
gray=all_num_algorithm(tree);

%ind_list=[];
gray_ind_list1=[];
gray_ind_list2=[];
prefactor=[];
subsi=[];
code_list=[];
in_out=[];
control_bit_change=[ones(1,2*s)];

%     for i=[1:2:n 2:2:n]
%         for l=1:n
%             j=l-i+1;
%             if j>=1
%                 %prefactor=[prefactor;1+0*(i*j)];  %Strength of interaction 
%                 prefactor=[prefactor;sqrt(i*j)];
%                 subsi=[subsi,i*j];
% 
%                 gray_num=gray([i+1,j]);
%                 graybin=dec2bin(gray_num-1,s)';
%                 graybin=graybin(:)';
%                 graybin=graybin-'0';
% 
%                 gray_num2=gray([i,j+1]);
%                 graybin2=dec2bin(gray_num2-1,s)';
%                 graybin2=graybin2(:)';
%                 graybin2=graybin2-'0';
% 
%                 gray_ind_list1=[gray_ind_list1; graybin ];      %qubits before
%                 gray_ind_list2=[gray_ind_list2; graybin2];      %qubits after operation
%                 code=(graybin~=graybin2);
%                 code_list=[code_list; code];                    %Which qubits are being changed in the interaction? 1=changed , 0 =stays
%                 c_1=find(code);
%                 in_out=[in_out;mod(sum(graybin(c_1)),2)];       %0=> 00<->11  |  1=> 01<->10 which 2 qubit subspace is transformed into which
%                 if size(gray_ind_list1,1)>1
%                     control_bit_change=[control_bit_change; mod(gray_ind_list1(end-1,:)+gray_ind_list1(end,:),2)];
%                 end
%             end
%         end
%     end   
    for i=[1:2:n 2:2:n]
        for l=1:n
            j=l-i+1;
            if j>=1
                %prefactor=[prefactor;1+0*(i*j)];  %Strength of interaction 
                prefactor=[prefactor;sqrt(i*j)];
                subsi=[subsi,i*j];

                gray_num=gray([j,i+1]);
                graybin=dec2bin(gray_num-1,s)';
                graybin=graybin(:)';
                graybin=graybin-'0';
                
                gray_num2=gray([j+1,i]);
                graybin2=dec2bin(gray_num2-1,s)';
                graybin2=graybin2(:)';
                graybin2=graybin2-'0';

                gray_ind_list1=[gray_ind_list1; graybin];%([s+1:2*s 1:s])];     %qubits before
                gray_ind_list2=[gray_ind_list2; graybin2];%([s+1:2*s 1:s])];      %qubits after operation
                code=(graybin~=graybin2);
                code_list=[code_list; code];                   %Which qubits are being changed in the interaction? 1=changed , 0 =stays
                c_1=find(code);
                in_out=[in_out;mod(sum(graybin(c_1)),2)];       %0=> 00<->11  |  1=> 01<->10 which 2 qubit subspace is transformed into which
                if size(gray_ind_list1,1)>1
                    control_bit_change=[control_bit_change; mod(gray_ind_list1(end-1,:)+gray_ind_list1(end,:),2)];
                end
            end
        end
        if mod(i,2)==1 %Uneven
            how_many_even=length(prefactor);
        end
    end   
    
code_list=[code_list;code_list(1,:)];
control_bit_change=[control_bit_change; mod(gray_ind_list1(end,:)+gray_ind_list1(1,:),2)];
gray_ind_list1=[gray_ind_list1; gray_ind_list1(1,:)];
gray_ind_list2=[gray_ind_list2; gray_ind_list2(1,:)];
prefactor=[prefactor;prefactor(1)];
in_out=[in_out;in_out(1)];
control_bit_change=((control_bit_change+code_list)>0);

%% How long
for j=1:2*s %Row for row - qbit for qbit
    %Determine repetitions
    k=find(control_bit_change(:,j)');               %At which step do the (control and flip) qubits change
    p=diff([k length(control_bit_change(:,j))+1]);    %how long do we keep the same control bit?
    c=arrayfun(@(x) x:-1:1,p,'un',0);
    how_long(:,j)=cat(2,c{:});
end

how_long(find(code_list))=0;            %Reduce how_long of qubits that are interacting as they need to be last in C^n stack
[~,order]=sort(how_long,2,'descend');   %Generate ordering for control qubits %2*s,2*s-1 are the interaction qubits, start stack at 1,2...
for i=2:size(code_list,1) %Find equal indexes that are switched without reason
     how_curr=how_long(i,:);
     un=unique(how_curr);
     
     for j=1:length(un)
         fi=find(how_curr==un(j));
         if length(fi)>1    
             curr_loc=[];
             prev_loc=[];
             for k=1:length(fi);
                curr_loc(k)=find(order(i,:)==fi(k));
                prev_loc(k)=find(order(i-1,:)==fi(k));
             end
             [p_fi b]=sort(prev_loc);
             order(i,curr_loc)=fi(b);
         end
     end
end
order(end,:)=order(1,:);
%% Add Gates
curr_order=zeros(1,2*s); %The current order of control gates
re=[2*s:-1:1 2*s+1:2*2*s-2];            %Index reverser ->1 on the right, 2*s on the left 
num_qubits=[1:2*s-4 2*s-3 2*s-3 ];   %how many Toffoli's you need for getting x (being the index here) control qubits
gray_ind_list1=gray_ind_list1(:,end:-1:1);

if s>1
    init=Create_Empty_Comp_Gate( name_in,2*s,2*s-3 );
    init=Generate_Gate_Circuit(init,circ_string_in,[1:2*s]);
    init.circuit_subs={};
    outit=Create_Empty_Comp_Gate( name_out,2*s,2*s-3 );
    outit=Generate_Gate_Circuit(outit,circ_string_out,[1:2*s]);
    outit.circuit_subs={};
    connector=Create_Empty_Comp_Gate( name_con,2*s,2*s-3 );
    connector=Generate_Gate_Circuit(connector,circ_string_con,[1:2*s]);
    connector.circuit_subs={};
    main=Create_Empty_Comp_Gate( name,2*s,2*s-3 );
    main=Generate_Gate_Circuit(main,circ_string,[1:2*s]);
    main.circuit_subs={};
    h_l=size(how_long,1); %Last step just to get back the connector settings
    for i=1:h_l
        if max(curr_order)==0 
            %% First interaction
            %Create all conditions via toffoli gates
            curr_order=re(order(i,:));

            index=[2*s+1,curr_order(1:2)];
            c_bit=gray_ind_list1(i,curr_order(1:2));
            init=Add_Toffoli(init,index,c_bit);
            init.circuit_subs{end+1}=[];

            for j=1:2*s-4
                index=[2*s+1+j,curr_order(j+2),2*s+j];
                c_bit=[gray_ind_list1(i,curr_order(j+2)) 1];
                init=Add_Toffoli(init,index,c_bit);
                init.circuit_subs{end+1}=[];
            end
            if prefactor(i)~=1
                subs={[{sym('phi'),sym('phi')*prefactor(i)/2}]};
                halving=1;
                subs2={[{sym('phi'),sym('phi')*prefactor(i)}]};
                halving2=0;
            else
                subs={[{sym('phi'),sym('phi')/2}]};
                halving=1;
                subs2={[{sym('phi'),sym('phi')}]};
                halving2=0;
            end
            
            if in_out(i)==1
                init=Add_Gate(init,{'C_XX_YY'},{[curr_order(end-1),curr_order(end),4*s-3]},subs);
                init=gray_gate_substitute(init,halving,subsi,i);
                connector=Add_Gate(connector,{'C_XX_YY'},{[curr_order(end-1),curr_order(end),4*s-3]},subs2);
                connector=gray_gate_substitute(connector,halving2,subsi,i);
                outit=Add_Gate(outit,{'C_XX_YY'},{[curr_order(end-1),curr_order(end),4*s-3]},subs);
                outit=gray_gate_substitute(outit,halving,subsi,i);
            else
                init=Add_Gate(init,{'Cx_XX_YY'},{[curr_order(end-1),curr_order(end),4*s-3]},subs);
                init=gray_gate_substitute(init,halving,subsi,i);
                connector=Add_Gate(connector,{'Cx_XX_YY'},{[curr_order(end-1),curr_order(end),4*s-3]},subs2);
                connector=gray_gate_substitute(connector,halving2,subsi,i);
                outit=Add_Gate(outit,{'Cx_XX_YY'},{[curr_order(end-1),curr_order(end),4*s-3]},subs);
                outit=gray_gate_substitute(outit,halving,subsi,i);
            end
        else
            %% Subsequent interactions
            %Which control bits have to change? (New ordering?)
            l=[];
            l2=(find(control_bit_change(i,:)));
            for k=1:length(l2)
                l=[l,find(order(i,:)==l2(k)),find(order(i-1,:)==l2(k))];
            end

            %Order isn't all - also need to check the control qubits
            index_change=min(l):2*s;       %These qubits in the order list have to be changed the earliest one, and oll that are built on top of it  
            cq_rem=length(index_change)-2;  %In total we have 2*s-3 toffoli's, we have to change the
                                            %control of length(index_change)-2 (two for the xx-yy operation
            cq_left=2*s-length(index_change);   %How many qubits are left
            if length(index_change)<=2  %Only change interaction, but all control qubits stay the same
                num_toff_back=0;
            else
                num_toff_back=num_qubits(cq_rem);   %Go back this many of the toffoli's 
            end
            
            % Removing control qubits
            counter=0;
            for j=2*s-4:-1:1
                counter=counter+1;
                if counter<=num_toff_back
                    index=[2*s+1+j,curr_order(j+2),2*s+j];
                    c_bit=[gray_ind_list1(i-1,curr_order(j+2)) 1];
                    if i<=how_many_even
                        init=Add_Toffoli(init,index,c_bit);
                        init.circuit_subs{end+1}=[];         
                        connector=Add_Toffoli(connector,index,c_bit);
                        connector.circuit_subs{end+1}=[];   
                        outit=Add_Toffoli(outit,index,c_bit);
                        outit.circuit_subs{end+1}=[];     
                    elseif i==how_many_even+1
                        init=Add_Toffoli(init,index,c_bit);
                        init.circuit_subs{end+1}=[];
                        connector=Add_Toffoli(connector,index,c_bit);
                        connector.circuit_subs{end+1}=[];   
                    else
                        main=Add_Toffoli(main,index,c_bit);
                        main.circuit_subs{end+1}=[];
                    end

                else
                    break; %quit for loop
                end
            end
            counter=counter+1;
            if counter<=num_toff_back
            	index=[2*s+1,curr_order(1:2)];
                c_bit=gray_ind_list1(i-1,curr_order(1:2));
                if i<=how_many_even
                    init=Add_Toffoli(init,index,c_bit);
                    init.circuit_subs{end+1}=[];      
                    connector=Add_Toffoli(connector,index,c_bit);
                    connector.circuit_subs{end+1}=[];      
                    outit=Add_Toffoli(outit,index,c_bit);
                    outit.circuit_subs{end+1}=[];      
                elseif i==how_many_even+1
                    init=Add_Toffoli(init,index,c_bit);
                    init.circuit_subs{end+1}=[];
                    connector=Add_Toffoli(connector,index,c_bit);
                    connector.circuit_subs{end+1}=[];   
                else
                	main=Add_Toffoli(main,index,c_bit);
                    main.circuit_subs{end+1}=[];                   
                end
                % Add new control qubits
                curr_order=re(order(i,:));
                index=[2*s+1,curr_order(1:2)];
                c_bit=gray_ind_list1(i,curr_order(1:2));
                if i<=how_many_even
                    init=Add_Toffoli(init,index,c_bit);
                    init.circuit_subs{end+1}=[];    
                    connector=Add_Toffoli(connector,index,c_bit);
                    connector.circuit_subs{end+1}=[];    
                    outit=Add_Toffoli(outit,index,c_bit);
                    outit.circuit_subs{end+1}=[];    
                elseif i==h_l
                    connector=Add_Toffoli_Left(connector,index,c_bit);
                    connector.circuit_subs={[],connector.circuit_subs{1:end}};    
                    outit=Add_Toffoli_Left(outit,index,c_bit);
                    outit.circuit_subs={[],outit.circuit_subs{1:end}};    
                else
                    main=Add_Toffoli(main,index,c_bit);
                    main.circuit_subs{end+1}=[];                    
                end
                counter=1;          %Reset counter for other direction
            else
                curr_order=re(order(i,:));
                counter=0;          %Reset counter for other direction
            end
            
            s_p=max([1,cq_left-1]); %Starting point for new control qubits depends on how far we went back  
            for j=s_p:2*s-4   %starting point depends on how far we went back 
                index=[2*s+1+j,curr_order(j+2),2*s+j];
                c_bit=[gray_ind_list1(i,curr_order(j+2)) 1];
                if i<=how_many_even
                    init=Add_Toffoli(init,index,c_bit);
                    init.circuit_subs{end+1}=[];
                    connector=Add_Toffoli(connector,index,c_bit);
                    connector.circuit_subs{end+1}=[];
                    outit=Add_Toffoli(outit,index,c_bit);
                    outit.circuit_subs{end+1}=[];
                elseif i==h_l
                    connector=Add_Toffoli_Left(connector,index,c_bit,counter+1);
                    connector.circuit_subs={connector.circuit_subs{1:counter-1},[],connector.circuit_subs{counter:end}};
                    outit=Add_Toffoli_Left(outit,index,c_bit,counter+1);
                    outit.circuit_subs={outit.circuit_subs{1:counter-1},[],outit.circuit_subs{counter:end}};
                else
                    main=Add_Toffoli(main,index,c_bit);
                    main.circuit_subs{end+1}=[];
                end
                counter=counter+1;
            end
            % Add new xx-yy interaction
            if prefactor(i)~=1
                subs={[{sym('phi'),sym('phi')*prefactor(i)/2}]};
                halving=1;
                halving2=0;
                subs2={[{sym('phi'),sym('phi')*prefactor(i)}]}; %Center piece needs to be a full rotation
            else
                subs={[{sym('phi'),sym('phi')/2}]};
                halving=1;
                subs2={[{sym('phi'),sym('phi')}]};
                halving2=0;
            end
            if in_out(i)==1
                if i<=how_many_even
                    init=Add_Gate(init,{'C_XX_YY'},{[curr_order(end-1),curr_order(end),4*s-3]},subs);
                    init=gray_gate_substitute(init,halving,subsi,i);
                    connector=Add_Gate(connector,{'C_XX_YY'},{[curr_order(end-1),curr_order(end),4*s-3]},subs2);
                    connector=gray_gate_substitute(connector,halving2,subsi,i);
                    outit=Add_Gate(outit,{'C_XX_YY'},{[curr_order(end-1),curr_order(end),4*s-3]},subs);
                    outit=gray_gate_substitute(outit,halving,subsi,i);
                elseif i<h_l
                    main=Add_Gate(main,{'C_XX_YY'},{[curr_order(end-1),curr_order(end),4*s-3]},subs2);
                    main=gray_gate_substitute(main,halving2,subsi,i);
                end
            else
                if i<=how_many_even
                    init=Add_Gate(init,{'Cx_XX_YY'},{[curr_order(end-1),curr_order(end),4*s-3]},subs);
                    init=gray_gate_substitute(init,halving,subsi,i);
                    connector=Add_Gate(connector,{'Cx_XX_YY'},{[curr_order(end-1),curr_order(end),4*s-3]},subs2);
                    connector=gray_gate_substitute(connector,halving2,subsi,i);
                    outit=Add_Gate(outit,{'Cx_XX_YY'},{[curr_order(end-1),curr_order(end),4*s-3]},subs);
                    outit=gray_gate_substitute(outit,halving,subsi,i);
                elseif i<h_l
                    main=Add_Gate(main,{'Cx_XX_YY'},{[curr_order(end-1),curr_order(end),4*s-3]},subs2);
                    main=gray_gate_substitute(main,halving2,subsi,i);
                end
            end
            if i==how_many_even
                %% Remove conditions on outit
                for j=2*s-4:-1:1
                    index=[2*s+1+j,curr_order(j+2),2*s+j];
                    c_bit=[gray_ind_list1(i,curr_order(j+2)) 1];
                    outit=Add_Toffoli(outit,index,c_bit);
                    outit.circuit_subs{end+1}=[];
                end
                index=[2*s+1,curr_order(1:2)];
                c_bit=gray_ind_list1(i,curr_order(1:2));
                outit=Add_Toffoli(outit,index,c_bit);
                outit.circuit_subs{end+1}=[];
            end
        end
    end
        
    %% Main 2 Connect 
        
elseif s==1
    %Only controlbits
    init=Create_Empty_Comp_Gate( name_in,2,0);
    outit=Create_Empty_Comp_Gate( name_out,2,0);
    connector=Create_Empty_Comp_Gate( name_con,2,0);
    main=Create_Empty_Comp_Gate( name,2,0);
    main=Generate_Gate_Circuit(main,circ_string,[1,2]);
    main=Add_Gate(main,{'XX_YY'},{[1,2]},{[]});
elseif s==0
    error('Zero qubits are too few.');
end
end