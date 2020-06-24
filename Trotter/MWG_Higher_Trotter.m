function [ gate ] = MWG_Higher_Trotter( names2,names,names_sum,connections,step_num,name,circ_string,sizle,anc_sizle,t )
%MWG_Higher_TROTTER does the Lie-Suzuki-Trotter decomposition for
%names2 -> full steps, names -> half steps, names_sum -> sum of two adjacent
%half steps
%step)
%Create sequence (seq) and step width d(seq) für multiple waveguides
%weights defines the interaction strengths of N modes   
%size(connections)=[2,interactions];

siz=sizle/2;
if siz~=round(siz)
    error('Size needs to be even qubits big.')
end

a=length(connections);
if length(connections)==1
    WG=connections;
    connections=[(1:WG-1)',(2:WG)'];
else
    WG=max(connections(:));
end
b=size(connections,1);

if size(names2,1)==1
    g={names2{1,:}};
    names2={g{:}};
    for i=2:size(connections,1)
        names2(end+1,:)=g;
    end
end
if size(names,1)==1
    g={names{1,:}};
    names={g{:}};
    for i=2:size(connections,1)
        names(end+1,:)=g;
    end
end
if size(names_sum,1)==1
    g={names_sum{1,:}};
    names_sum={g{:}};
    for i=2:size(connections,1)
        names_sum(end+1,:)=g;
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
    gate_steps={gate_steps{:},names{j}};
    index_steps={index_steps{:},in};
    param_steps={param_steps{:},[{}]};
end
for i=1:step_num %Normal Trotter steps
    for k=1:size(names2,2) %steps of recursive sheme
        for j=2:b-1
            in=[(connections(j,1)-1)*siz+(1:siz),(connections(j,2)-1)*siz+(1:siz),anc_s];
            gate_steps={gate_steps{:},names{j,k}};
            index_steps={index_steps{:},in};
            param_steps={param_steps{:},[{}]};
        end
        for j=b:-1:1
            if j==b
                in=[(connections(j,1)-1)*siz+(1:siz),(connections(j,2)-1)*siz+(1:siz),anc_s];
                gate_steps={gate_steps{:},names2{j,k}};
                index_steps={index_steps{:},in};
                param_steps={param_steps{:},[{}]};
            elseif j~=1
                in=[(connections(j,1)-1)*siz+(1:siz),(connections(j,2)-1)*siz+(1:siz),anc_s];
                gate_steps={gate_steps{:},names{j,k}};
                index_steps={index_steps{:},in};
                param_steps={param_steps{:},[{}]};
            else %j==1 && 1~=length(weights)
                if ((i==step_num) && (k==size(names2,2)))==0
                    in=[(connections(j,1)-1)*siz+(1:siz),(connections(j,2)-1)*siz+(1:siz),anc_s];
                    gate_steps={gate_steps{:},names_sum{j,k}};
                    index_steps={index_steps{:},in};
                    param_steps={param_steps{:},[{}]};            
                end
            end
        end
    end
end
if b>1
    j=1;
    k=size(names2,2);
    in=[(connections(j,1)-1)*siz+(1:siz),(connections(j,2)-1)*siz+(1:siz),anc_s];
    gate_steps={gate_steps{:},names{j,k}};
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