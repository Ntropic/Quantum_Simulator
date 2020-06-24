function [ circuit ] = Add_Toffoli( circuit,index,c_bits )
%ADD_TOFFOLI adds a Toffoli gate to a circuit 
%   circuit is the previously established circuit
%   index are the indexes for the control bits and the not-bit
%   [not_ind,ind_1,ind_2] or [not_ind,ind_1]
%   c_bits are the control bits (00,01,10,11)
    
if length(index)==3
    if length(c_bits)==2
        %Add Toffoli
        if c_bits(1:2)==[1 1]
            circuit.steps.gates={circuit.steps.gates{:},'toffoli11'};
            circuit.steps.index={circuit.steps.index{:},index};
            circuit.steps.param={circuit.steps.param{:},[]};
        elseif c_bits(1:2)==[0 1]
            circuit.steps.gates={circuit.steps.gates{:},'toffoli01'};
            circuit.steps.index={circuit.steps.index{:},index};
            circuit.steps.param={circuit.steps.param{:},[]};
        elseif c_bits(1:2)==[1 0]
            circuit.steps.gates={circuit.steps.gates{:},'toffoli01'};
            circuit.steps.index={circuit.steps.index{:},index([1,3,2])};
            circuit.steps.param={circuit.steps.param{:},[]};
        elseif c_bits(1:2)==[0 0]
            circuit.steps.gates={circuit.steps.gates{:},'toffoli00'};
            circuit.steps.index={circuit.steps.index{:},index};
            circuit.steps.param={circuit.steps.param{:},[]};
        else
            error('Sequence has to consist of 0s and 1s.')
        end
    else
        error('length(index) must be 2.');
    end
elseif length(index)==2
    if length(c_bits)==1
        %Add CNOT
        if c_bits(1)==1
            circuit.steps.gates={circuit.steps.gates{:},'CNOT1'};
            circuit.steps.index={circuit.steps.index{:},index};
            circuit.steps.param={circuit.steps.param{:},[]};
        elseif c_bits(1)==0
            circuit.steps.gates={circuit.steps.gates{:},'CNOT0'};
            circuit.steps.index={circuit.steps.index{:},index};
            circuit.steps.param={circuit.steps.param{:},[]};
        else
            error('Sequence has to consist of 0s and 1s.')
        end
    else
        error('length(index) must be 1.');
    end
else
    error('length(index) must be 3 or 2 (for CNOTs).');
end
if any(strcmp(fields(circuit),'step_num'))
    circuit.step_num=circuit.step_num+1;
else
    circuit.step_num=length(circuit.steps.gates);
end

end

