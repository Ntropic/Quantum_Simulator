function [ init main connector outit gate_names ] = Split_Strang_Pauli_Ladder_Gates( n,name,circ_string )
%SPLIT_STRANG_PAULI_LADDER_GATES (Strang Splitting - Symmetric) Creates the 
%gates for Pauli Gray encoding and 2 mode interactions for n photons but 
%splits the gate into a init gate that creates the conditions on the 
%ancilla qubits with Strang Splitting the main gate and an outit, 
%that gets rid of the conditions the gate can be made symmetric 
%via mode==1 -> Strang Splitting
if nargin==3
    circ_string_in= [circ_string '^{(1),(i)}'];
    circ_string_out=[circ_string '^{(1),(f)}'];
    circ_string_con=[circ_string '^{(1),(c)}'];
    circ_string=    [circ_string '^{(1),(m)}'];
    name_in=[name '_in'];
    name_out=[name '_out'];
    name_con=[name '_con'];
    name=[name '_cen'];
elseif nargin==2
    circ_string_in= ['S_{L,n={' num2str(n) '}}^{(1),(i)}'];
    circ_string_out=['S_{L,n={' num2str(n) '}}^{(1),(f)}'];
    circ_string_con=['S_{L,n={' num2str(n) '}}^{(1),(c)}'];
    circ_string=    ['S_{L,n={' num2str(n) '}}^{(1),(m)}'];
    name_in=[name '_in'];
    name_out=[name '_out'];
    name_con=[name '_con'];
    name=[name '_cen'];
elseif nargin==1
    circ_string_in= ['S_{L,n={' num2str(n) '}}^{(1),(i)}'];
    circ_string_out=['S_{L,n={' num2str(n) '}}^{(1),(f)}'];
    circ_string_con=['S_{L,n={' num2str(n) '}}^{(1),(c)}'];
    circ_string=    ['S_{L,n={' num2str(n) '}}^{(1),(m)}'];
    name_in=['Ladder_n_' num2str(n) '_sym_in'];
    name_out=['Ladder_n_' num2str(n) '_sym_out'];
    name_con=['Ladder_n_' num2str(n) '_sym_con'];
    name=['Ladder_n_' num2str(n) '_sym_cen'];
end
gate_names={name_in,name_out,name,name_con};

s=(n+1)*2; %n+1 qubits for n photons *2 modes

init=Create_Empty_Comp_Gate( name_in,s,0);
init=Generate_Gate_Circuit(init,circ_string_in,[1:s]);
outit=Create_Empty_Comp_Gate( name_out,s,0 );
outit=Generate_Gate_Circuit(outit,circ_string_out,[1:s]);
connector=Create_Empty_Comp_Gate( name_con,s,0 );
connector=Generate_Gate_Circuit(connector,circ_string_con,[1:s]);
main=Create_Empty_Comp_Gate( name,s,0 );
main=Generate_Gate_Circuit(main,circ_string,[1:s]);

phi=sym('phi');

%Difference (d) based model for the ladder
block_nums=[];
prefactor=[];
subsi=[];

%Even uneven ordering model
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

subsi2=[];
for i=1:size(block_nums,1)
    subs={[{phi,phi*prefactor(i)}]};
    subs2={[{phi,phi*prefactor(i)/2}]};

    i2=block_nums(i,1);
    j2=block_nums(i,2);
    if i<=how_many_even   
        %Half
        init=Add_Gate(init,{'L_{+}'},{[i2 i2+1 n+1+j2 n+2+j2]},subs2);
        if abs(sqrt(subsi(i))-round(sqrt(subsi(i))))<10^-14
            if mod(sqrt(subsi(i))/2,1)<10^-14
                if abs(sqrt(subsi(i))/2-1)<10^-14
                    init.circuit_subs{i}=[];
                else
                    init.circuit_subs{i}=[{'\phi',[ num2str(round(sqrt(subsi(i))/2)) '\phi' ]}];
                end 
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
        
        %Full
        connector=Add_Gate(connector,{'L_{+}'},{[i2 i2+1 n+1+j2 n+2+j2]},subs);	
        if subsi(i)~=1 
            if abs(sqrt(subsi(i))-round(sqrt(subsi(i))))<10^-14
                connector.circuit_subs{i}=[{'\phi',[num2str(round(sqrt(subsi(i)))) '\phi' ]}];
            else
                connector.circuit_subs{i}=[{'\phi',['\sqrt{' num2str(subsi(i)) '}\phi']}];
            end
        else
            connector.circuit_subs{i}=[];
        end
    else
        %Full
        main=Add_Gate(main,{'L_{+}'},{[i2 i2+1 n+1+j2 n+2+j2]},subs);	
        if subsi(i)~=1 
            if abs(sqrt(subsi(i))-round(sqrt(subsi(i))))<10^-14
                main.circuit_subs{i-how_many_even}=[{'\phi',[num2str(round(sqrt(subsi(i)))) '\phi' ]}];
            else
                main.circuit_subs{i-how_many_even}=[{'\phi',['\sqrt{' num2str(subsi(i)) '}\phi']}];
            end
        else
            main.circuit_subs{i-how_many_even}=[];
        end
    end
end
k=0;
for i=how_many_even:-1:1
    k=k+1;
    subs2={[{phi,phi*prefactor(i)/2}]};

    i2=block_nums(i,1);
    j2=block_nums(i,2);
    if i<=how_many_even   
        %Half
        outit=Add_Gate(outit,{'L_{+}'},{[i2 i2+1 n+1+j2 n+2+j2]},subs2);
        if abs(sqrt(subsi(i))-round(sqrt(subsi(i))))<10^-14
            if mod(sqrt(subsi(i))/2,1)<10^-14
                if abs(sqrt(subsi(i))/2-1)<10^-14
                    outit.circuit_subs{k}=[];
                else
                    outit.circuit_subs{k}=[{'\phi',[ num2str(round(sqrt(subsi(i))/2)) '\phi' ]}];
                end 
            else
                outit.circuit_subs{k}=[{'\phi',['\frac{' num2str(round(sqrt(subsi(i)))) '}{2}\phi' ]}];
            end
        else
            if mod(sqrt(subsi(i))/2,1)<10^-14
                outit.circuit_subs{k}=[{'\phi',['\sqrt{' num2str(round(sqrt(subsi(i))/2)) '}\phi']}];
            else    
                outit.circuit_subs{k}=[{'\phi',['\frac{\sqrt{' num2str(subsi(i)) '}}{2}\phi']}];
            end
        end   
    end
end    
end