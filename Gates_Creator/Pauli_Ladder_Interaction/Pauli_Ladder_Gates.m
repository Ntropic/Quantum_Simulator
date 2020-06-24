function [ gate ] = Pauli_Ladder_Gates( n,mode,name,circ_string )
%PAULI_LADDER_GATES Creates the gates for Pauli ladder encoding and 2 mode
%interactions for n photons
% the gate can be made symmetric via mode==1 -> Strang Splitting
if nargin==3
    if mode==0
        circ_string=['L_{n=' num2str(n) '}'];
    else
        circ_string=['S^{(' num2str(n) ')}_{L,1}'];
    end
elseif nargin==2
    if mode==0
        circ_string=['L_{n=' num2str(n) '}'];
        name=['Ladder_n_' num2str(n) '_step'];
    else
        circ_string=['S^{(' num2str(n) ')}_{L,1}'];
        name=['Ladder_n_' num2str(n) '_sym'];
    end
elseif nargin==1
    mode=0; %Non symmetric process#
    name=['Ladder_n_' num2str(n) '_step'];
    circ_string=['L_{n=' num2str(n) '}'];
end
phi=sym('phi');

s=(n+1)*2; %n+1 qubits for n photons *2 modes

block_nums=[];
prefactor=[];
subsi=[];
%% Difference (d) based model for the ladder
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
%% Even uneven model for the ladder (generalized xy approach)
% for n_now=1:n
%     for i=1:2:n_now
%         j=n_now+1-i;
%         block_nums=[block_nums;i j];
%         prefactor=[prefactor;2*sqrt(i*j)];  %Strength of interaction 
%         subsi=[subsi;i*j];
%     end
%     for i=2:2:n_now
%         j=n_now+1-i;
%         block_nums=[block_nums;i j];
%         prefactor=[prefactor;2*sqrt(i*j)];  %Strength of interaction 
%         subsi=[subsi;i*j];
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
for n_now=1:n
    for i=2:2:n_now
        j=n_now+1-i;
        block_nums=[block_nums;i j];
        prefactor=[prefactor;2*sqrt(i*j)];  %Strength of interaction 
        subsi=[subsi;i*j];
    end
end

gate=Create_Empty_Comp_Gate(name,s,0);
gate=Generate_Gate_Circuit(gate,circ_string,[1:s]);
for i=1:size(block_nums,1)
    if prefactor(i)~=1
        if mode==0
            subs={[{phi,phi*prefactor(i)}]};
        else
            if i~=length(prefactor)
                subs={[{phi,phi*prefactor(i)/2}]};
            else
                subs={[{phi,phi*prefactor(i)}]};
            end
        end
    else
        if mode==0
            subs={[]};
        else
            if i~=length(prefactor)
                subs={[{phi,phi/2}]};
            else
                subs={[{phi,phi*2}]};
            end
        end
    end
    i2=block_nums(i,1);
    j2=block_nums(i,2);
    gate=Add_Gate(gate,{'L_{+}'},{[i2 i2+1 n+1+j2 n+2+j2]},subs);
    if subsi(i)~=1
        if abs(sqrt(subsi(i))-round(sqrt(subsi(i))))<10^-14
            if mode==0  || i==size(block_nums,1)
                gate.circuit_subs{i}=[{'\phi',[num2str(round(sqrt(subsi(i)))) '\phi' ]}];
            else
                gate.circuit_subs{i}=[{'\phi',['\frac{' num2str(round(sqrt(subsi(i)))) '\phi}{2}' ]}];
            end
        else
            if mode==0 || i==size(block_nums,1)
            	gate.circuit_subs{i}=[{'\phi',['\sqrt{' num2str(subsi(i)) '}\phi']}];
            else
                gate.circuit_subs{i}=[{'\phi',['\frac{\sqrt{' num2str(subsi(i)) '}\phi}{2}' ]}];
            end
        end
    else
        if mode==0 || i==size(block_nums,1)
            gate.circuit_subs{i}=[];
        else
            gate.circuit_subs{i}=[{'\phi','\frac{\phi}{2}'}];
        end
    end
end
if mode~=0 %Add symmetrisation
    for i=size(block_nums,1)-1:-1:1
        if prefactor(i)*2~=1
            subs={[{phi,phi*prefactor(i)/2}]};
        else
            subs={[{phi,phi/2}]};
        end
        i2=block_nums(i,1);
        j2=block_nums(i,2);
        gate=Add_Gate(gate,{'L_{+}'},{[i2 i2+1 n+1+j2 n+2+j2]},subs);
    if subsi(i)~=1
        if abs(sqrt(subsi(i))-round(sqrt(subsi(i))))<10^-14
            gate.circuit_subs{end+1}=[{'\phi',['\frac{' num2str(round(sqrt(subsi(i)))) '\phi}{2}' ]}];
        else
            gate.circuit_subs{end+1}=[{'\phi',['\frac{\sqrt{' num2str(subsi(i)) '}\phi}{2}' ]}];
        end
    else
        gate.circuit_subs{end+1}=[{'\phi','\frac{\phi}{2}'}];
    end
end

end