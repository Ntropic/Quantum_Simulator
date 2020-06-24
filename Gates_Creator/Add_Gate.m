function [ gate ] = Add_Gate( gate,gate_steps,index_steps,param_steps )
%ADD_GATE Adds subgate to a composite gate (comp_gate)

if length(gate_steps)==length(index_steps)  
    gate.steps.gates={gate.steps.gates{:},gate_steps{:}};
    gate.steps.index={gate.steps.index{:},index_steps{:}};
    if nargin==4
        if length(gate_steps)==length(param_steps)
            gate.steps.param={gate.steps.param{:},param_steps{:}};
        else
            error('Gate steps and Parameter Steps must be of same length! No output was created.');
        end
    else
        for i=1:length(index_steps)
            gate.steps.param={gate.steps.param{:},[]};
        end
    end
else
    error('Gate steps and Index Steps must be of same length! No output was created.');
end

gate.step_num=length(gate.steps.param);   %Change step_num to encompass new subgates

size=gate.size;
anc_size=gate.anc_size;

max_ind=max([index_steps{:}]);
if max_ind>anc_size+size
    fprintf('Gate uses more ancilla qubits, corrected description of gate (Add_Gate)');
    gate.anc_size=max_ind-size;
end
end

