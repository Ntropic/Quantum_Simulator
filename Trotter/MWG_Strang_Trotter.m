function [ gate ] = MWG_Strang_Trotter( gate_names,connections,step_num,name,circ_string,sizle,anc_sizle,t )
%MWG_STRANG TROTTER_INOUT does the Lie-Suzuki-Trotter decomposition for
%gates_names={...} -> size(gates_names)=(weights,2) ->2=(1/2 step, full
%step)
%Create sequence (seq) and step width d(seq) für multiple waveguides
%weights defines the interaction strengths of N modes   
%size(connections)=[2,interactions];

siz=sizle/2;
if siz~=round(siz)
    error('Size needs to be even qubits big.')
end

a=length(connections);
if a==1
    WG=connections;
    connections=[(1:WG-1)',(2:WG)'];
else
    WG=max(connections(:));
end
b=size(connections,1);

if size(gate_names,1)==1
    g={gate_names{1,:}};
    gate_names={g{:}};
    for i=2:size(connections,1)
        gate_names(end+1,:)=g;
    end
end

s=siz*WG;
anc_s=(1:anc_sizle)+s;

gate_steps={};
index_steps={};
param_steps={};

if b>1
    j=1;
    in=[(connections(j,1)-1)*siz+(1:siz),(connections(j,2)-1)*siz+(1:siz),anc_s];
    gate_steps={gate_steps{:},gate_names{j,1}};
    index_steps={index_steps{:},in};
    param_steps={param_steps{:},[{}]};
end
for i=1:step_num
    for j=2:b-1
        in=[(connections(j,1)-1)*siz+(1:siz),(connections(j,2)-1)*siz+(1:siz),anc_s];
        gate_steps={gate_steps{:},gate_names{j,1}};
        index_steps={index_steps{:},in};
        param_steps={param_steps{:},[{}]};
    end
    for j=b:-1:1
        if j==b
            in=[(connections(j,1)-1)*siz+(1:siz),(connections(j,2)-1)*siz+(1:siz),anc_s];
            gate_steps={gate_steps{:},gate_names{j,2}};
            index_steps={index_steps{:},in};
            param_steps={param_steps{:},[{}]};
        elseif j~=1
            in=[(connections(j,1)-1)*siz+(1:siz),(connections(j,2)-1)*siz+(1:siz),anc_s];
            gate_steps={gate_steps{:},gate_names{j,1}};
            index_steps={index_steps{:},in};
            param_steps={param_steps{:},[{}]};
        else %j==1 && 1~=length(weights)
            if i~=step_num
                in=[(connections(j,1)-1)*siz+(1:siz),(connections(j,2)-1)*siz+(1:siz),anc_s];
                gate_steps={gate_steps{:},gate_names{j,2}};
                index_steps={index_steps{:},in};
                param_steps={param_steps{:},[{}]};            
            end
        end
    end
end
if b>1
    j=1;
    in=[(connections(j,1)-1)*siz+(1:siz),(connections(j,2)-1)*siz+(1:siz),anc_s];
    gate_steps={gate_steps{:},gate_names{j,1}};
    index_steps={index_steps{:},in};
    param_steps={param_steps{:},[{}]};
end



if length(name)==0
    name=['T_{N=' num2str(WG) ',n=' num2str(sizle) '}^{(' num2str(length(step_num)) ')}(' latex(sym(t)) ')'];
end
if length(circ_string)==0
    circ_string=['T_{N=' num2str(WG) ',n=' num2str(sizle) '}^{(' num2str(length(step_num)) ')}(' latex(sym(t)) ')'];
end

gate=Create_Comp_Gate(name,s,anc_sizle,gate_steps,index_steps,param_steps);
gate=Generate_Gate_Circuit(gate,circ_string,1:sizle,[],[]);
end