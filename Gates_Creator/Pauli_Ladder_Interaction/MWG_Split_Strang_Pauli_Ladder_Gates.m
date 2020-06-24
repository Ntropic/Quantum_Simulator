function [ init outit init2 outit2 symit con  gate_names ] = MWG_Split_Strang_Pauli_Ladder_Gates( n,name,circ_string )
%MWG_SPLIT_STRANG_PAULI_LADDER_GATES (Strang Splitting - Symmetric) Creates the 
%gates for Pauli Gray encoding and multi mode interactions for n photons but 
%splits the gate into a init gate that creates the conditions on the 
%ancilla qubits with Strang Splitting the main gate and an outit, 
%that gets rid of the conditions the gate can be made symmetric 
%via mode==1 -> Strang Splitting
if nargin==3
    circ_string_in= [circ_string '^{(1)}'];
    circ_string_out=[circ_string '^{(1)\dagger}'];
    circ_string_in2= [circ_string '^{(1),(r)}'];
    circ_string_out2=[circ_string '^{(1),(r)\dagger}'];
    circ_string_sym=[circ_string '^{(1),(s)}'];
    circ_string_con=[circ_string '^{(1),(c)}'];
    name_in=[name '_in'];
    name_out=[name '_out'];
    name_in2=[name '_in2'];
    name_out2=[name '_out2'];
    name_sym=[name '_sym'];
    name_conn=[name '_con'];
elseif nargin==2
    circ_string_in= ['S_{L,n={' num2str(n) '}}^{(1)}'];
    circ_string_out=['S_{L,n={' num2str(n) '}}^{(1)\dagger}'];
    circ_string_in2= ['S_{L,n={' num2str(n) '}}^{(1),(r)}'];
    circ_string_out2=['S_{L,n={' num2str(n) '}}^{(1),(r)\dagger}'];
    circ_string_sym=['S_{L,n={' num2str(n) '}}^{(1),(s)}'];
    circ_string_con=['S_{L,n={' num2str(n) '}}^{(1),(c)}'];
    name_in=[name '_in'];
    name_out=[name '_out'];
    name_in2=[name '_in2'];
    name_out2=[name '_out2'];
    name_sym=[name '_sym'];
    name_conn=[name '_con'];
elseif nargin==1
    circ_string_in= ['S_{L,n={' num2str(n) '}}^{(1)}'];
    circ_string_out=['S_{L,n={' num2str(n) '}}^{(1)\dagger}'];
    circ_string_in2= ['S_{L,n={' num2str(n) '}}^{(1),(r)}'];
    circ_string_out2=['S_{L,n={' num2str(n) '}}^{(1),(r)\dagger}'];
    circ_string_sym=['S_{L,n={' num2str(n) '}}^{(1),(s)}'];
    circ_string_con=['S_{L,n={' num2str(n) '}}^{(1),(c)}'];
    name_in=['Ladder_n_' num2str(n) '_in'];
    name_out=['Ladder_n_' num2str(n) '_out'];
    name_in2=['Ladder_n_' num2str(n) '_in2'];
    name_out2=['Ladder_n_' num2str(n) '_out2'];
    name_sym=['Ladder_n_' num2str(n) '_sym'];
    name_con=['Ladder_n_' num2str(n) '_con'];
end
gate_names={name_in,name_in2,name_out,name_out2,name_sym,name_con};

s=(n+1)*2; %n+1 qubits for n photons *2 modes

init=Create_Empty_Comp_Gate( name_in,s,0);
init=Generate_Gate_Circuit(init,circ_string_in,[1:s]);
outit=Create_Empty_Comp_Gate( name_out,s,0 );
outit=Generate_Gate_Circuit(outit,circ_string_out,[1:s]);
init2=Create_Empty_Comp_Gate( name_in2,s,0);
init2=Generate_Gate_Circuit(init2,circ_string_in2,[1:s]);
outit2=Create_Empty_Comp_Gate( name_out2,s,0 );
outit2=Generate_Gate_Circuit(outit2,circ_string_out2,[1:s]);
symit=Create_Empty_Comp_Gate( name_sym,s,0 );
symit=Generate_Gate_Circuit(symit,circ_string_sym,[1:s]);
con=Create_Empty_Comp_Gate( name_con,s,0 );
con=Generate_Gate_Circuit(con,circ_string_con,[1:s]);

phi=sym('phi');

%Difference (d) based model for the ladder
block_nums=[];
prefactor=[];
subsi=[];
% for d=0:n
%     i_max=floor((n+d+1)/2);
%     for i=d+1:i_max
%         j=i-d;
%         block_nums=[block_nums;i j];
%         prefactor=[prefactor;2*sqrt(i*j)];  %Strength of interaction 
%         subsi=[subsi;i*j];
%     end
%     if d~=0
%         for i=d+1:i_max
%             j=i-d;
%             block_nums=[block_nums;j i];
%             prefactor=[prefactor;2*sqrt(i*j)];  %Strength of interaction 
%             subsi=[subsi;i*j];
%         end
%     end
% end

for n_now=1:n
    for i=1:2:n_now
        j=n_now+1-i;
        block_nums=[block_nums;i j];
        prefactor=[prefactor;2*sqrt(i*j)];  %Strength of interaction 
        subsi=[subsi;i*j];
    end
end
how_many_even=length(prefactor);
for n_now=1:n
    for i=2:2:n_now
        j=n_now+1-i;
        block_nums=[block_nums;i j];
        prefactor=[prefactor;2*sqrt(i*j)];  %Strength of interaction 
        subsi=[subsi;i*j];
    end
end

for i=1:size(block_nums,1)
    if prefactor(i)~=1
        subs={[{phi,phi*prefactor(i)/2}]};
        subs2={[{phi,phi*prefactor(i)}]};
    else
        subs={[{phi,phi/2}]};
        subs2={[]};
    end
    i2=block_nums(i,1);
    j2=block_nums(i,2);
    
    if i<=how_many_even
        init=Add_Gate(init,{'L_{+}'},{[i2 i2+1 n+1+j2 n+2+j2]},subs);
        if subsi(i)~=1
            if abs(sqrt(subsi(i))/2-1)>10^-14
                if abs(sqrt(subsi(i))-round(sqrt(subsi(i))))<10^-14
                    if mod(sqrt(subsi(i))/2,1)<10^-14
                        init.circuit_subs{i}=[{'\phi',[ num2str(round(sqrt(subsi(i))/2)) '\phi' ]}];
                    else
                        init.circuit_subs{i}=[{'\phi',['\frac{' num2str(round(sqrt(subsi(i)))) '}{2}\phi' ]}];
                    end
                else
                    if mod(sqrt(subsi(i))/2,1)<10^-14
                        init.circuit_subs{i}=[{'\phi',['\sqrt{' num2str(round(sqrt(subsi(i))/2)) '}\phi']}];
                    else    
                        init.circuit_subs{i}=[{'\phi',['\frac{\sqrt{' num2str(subsi(i)) '}}{2}\phi']}];
                    end
                end
            else
                init.circuit_subs{i}=[];
            end
        else
            init.circuit_subs{i}=[{'\phi','\frac{\phi}{2}'}];
        end
    end
    
    
    if i>how_many_even
        init2=Add_Gate(init2,{'L_{+}'},{[i2 i2+1 n+1+j2 n+2+j2]},subs);
       	if subsi(i)~=1
            if abs(sqrt(subsi(i))/2-1)>10^-14
                if abs(sqrt(subsi(i))-round(sqrt(subsi(i))))<10^-14
                    if mod(sqrt(subsi(i))/2,1)<10^-14
                        init2.circuit_subs{end+1}=[{'\phi',[ num2str(round(sqrt(subsi(i))/2)) '\phi' ]}];
                    else
                        init2.circuit_subs{end+1}=[{'\phi',['\frac{' num2str(round(sqrt(subsi(i)))) '}{2}\phi' ]}];
                    end
                else
                    if mod(sqrt(subsi(i))/2,1)<10^-14
                        init2.circuit_subs{end+1}=[{'\phi',['\sqrt{' num2str(round(sqrt(subsi(i))/2)) '}\phi']}];
                    else    
                        init2.circuit_subs{end+1}=[{'\phi',['\frac{\sqrt{' num2str(subsi(i)) '}}{2}\phi']}];
                    end
                end
            else
                init2.circuit_subs{end+1}=[];
            end
        else
            init2.circuit_subs{i-1}=[{'\phi','\frac{\phi}{2}'}];
        end
    end
    
    if i<=how_many_even
        con=Add_Gate(con,{'L_{+}'},{[i2 i2+1 n+1+j2 n+2+j2]},subs2);
    	if subsi(i)~=1
            if abs(sqrt(subsi(i))-round(sqrt(subsi(i))))<10^-14
                con.circuit_subs{end+1}=[{'\phi',[num2str(round(sqrt(subsi(i)))) '\phi' ]}];
            else
                con.circuit_subs{end+1}=[{'\phi',['\sqrt{' num2str(subsi(i)) '}\phi']}];
            end
        else
            con.circuit_subs{end+1}=[];
        end
    end
    
    if i>how_many_even
        symit=Add_Gate(symit,{'L_{+}'},{[i2 i2+1 n+1+j2 n+2+j2]},subs2);
    	if subsi(i)~=1
            if abs(sqrt(subsi(i))-round(sqrt(subsi(i))))<10^-14
                symit.circuit_subs{end+1}=[{'\phi',[num2str(round(sqrt(subsi(i)))) '\phi' ]}];
            else
                symit.circuit_subs{end+1}=[{'\phi',['\sqrt{' num2str(subsi(i)) '}\phi']}];
            end
        else
            symit.circuit_subs{end+1}=[];
        end
    end
end
outit.steps=init.steps;
outit.step_num=init.step_num;
outit.circuit_subs=init.circuit_subs;

outit2.steps=init2.steps;
outit2.step_num=init2.step_num;
outit2.circuit_subs=init2.circuit_subs;
%outit2.circuit_subs=init2.circuit_subs(end:-1:1);

end