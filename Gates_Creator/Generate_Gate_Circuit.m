function [ gate ] = Generate_Gate_Circuit( gate,gate_string,gate_indexes,conditions,conditions_indexes,circuit_subs )
%GENERATE_GATE_CIRCUIT creates a string for creating multiqubit gates and
%conditions 
%   -> the gate to which we want to add the circuit representation/label
%   -> gate_string          -   gives string for the box of the multigate
%   -> gate_indexes         -   which indexes are in the box 
%   -> conditions           -   1 or 0 conditional circals (like in CNot)
%   -> conditions_indexes   -   which indexes get the conditional circles
if nargin==2
    error('Gate indexes need to be specified.')
end


if exist('conditions_indexes')
    %Remove overlap of gates- and conditions-indexes:
    for i=1:length(gate_indexes)
        for j=1:length(conditions_indexes)
            if gate_indexes(i)==conditions_indexes(j)
                gate_indexes(i)=[];
            end
        end
    end
else
    conditions_indexes=[];
end

%% Generate String
index_list=[sort(gate_indexes(:),1,'descend');sort(conditions_indexes(:),1,'descend')];
index_list=[index_list';ones(1,size(index_list,1))];

if length(gate_indexes)==0
    circuit={[],[]};
elseif length(gate_indexes)>1   %multigate and sgate
    multi=['\multigate{' num2str(length(gate_indexes)-1) '}{' gate_string '}'];
    ghost=['\ghost{' gate_string '}'];
    sgate=['\sgate{' gate_string '}{_' num2str(length(gate_indexes)-1) '_}'];
    nosgate=['\gate{' gate_string '}'];
    ctrl='\ctrl{_1_}';
    ctrl0='\ctrl0{_1_}';
    
    circuit={index_list,multi};
    circuit2={index_list,nosgate};
    for i=2:length(gate_indexes)
        circuit={circuit{:},ghost};
        sgate=['\sgate{' gate_string '}'];
        circuit2={circuit2{:},sgate};
    end
    for i=1:length(conditions_indexes)
        if conditions(i)==1
            circuit={circuit{:},ctrl};
            circuit2={circuit2{:},ctrl};
        else
            circuit={circuit{:},ctrl0};
            circuit2={circuit2{:},ctrl0};
        end
    end
    circuit={circuit{:};circuit2{:}};
else                        %gate
    nosgate=['\gate{' gate_string '}'];
    ctrl='\ctrl{_1_}';
    ctrl0='\ctrl0{_1_}';
    
    circuit={index_list,nosgate};
    for i=1:length(conditions_indexes)
        if conditions(i)==1
            circuit={circuit{:},ctrl};
        else
            circuit={circuit{:},ctrl0};
        end
    end
end
gate.circuit=circuit;
if exist('circuit_subs')
    gate.circuit_subs=circuit_subs;
end
end

