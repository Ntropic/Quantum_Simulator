function [ circuit ] = Add_Toffoli_Left( circuit,index,c_bits,how_many )
%ADD_TOFFOLI_LEFT adds a Toffoli gate to a circuit from the left 
%(Add_Toffoli adds it to the right) 
%   circuit is the previously established circuit
%   index are the indexes for the control bits and the not-bit
%   [not_ind,ind_1,ind_2] or [not_ind,ind_1]
%   c_bits are the control bits (00,01,10,11)
%   how_many gates are to the left of the new Toffoli +1
if nargin==3
    how_many=1;
end
if length(index)==3
    if length(c_bits)==2
        %Add Toffoli
        if c_bits(1:2)==[1 1]
            circuit.steps.gates={circuit.steps.gates{1:how_many-1},'toffoli11',circuit.steps.gates{how_many:end}};
            circuit.steps.index={circuit.steps.index{1:how_many-1},index,circuit.steps.index{how_many:end}};
            circuit.steps.param={circuit.steps.param{1:how_many-1},[],circuit.steps.param{how_many:end}};
        elseif c_bits(1:2)==[0 1]
            circuit.steps.gates={circuit.steps.gates{1:how_many-1},'toffoli01',circuit.steps.gates{how_many:end}};
            circuit.steps.index={circuit.steps.index{1:how_many-1},index,circuit.steps.index{how_many:end}};
            circuit.steps.param={circuit.steps.param{1:how_many-1},[],circuit.steps.param{how_many:end}};
        elseif c_bits(1:2)==[1 0]
            circuit.steps.gates={circuit.steps.gates{1:how_many-1},'toffoli01',circuit.steps.gates{how_many:end}};
            circuit.steps.index={circuit.steps.index{1:how_many-1},index([1,3,2]),circuit.steps.index{how_many:end}};
            circuit.steps.param={circuit.steps.param{1:how_many-1},[],circuit.steps.param{how_many:end}};
        elseif c_bits(1:2)==[0 0]
            circuit.steps.gates={circuit.steps.gates{1:how_many-1},'toffoli00',circuit.steps.gates{how_many:end}};
            circuit.steps.index={circuit.steps.index{1:how_many-1},index,circuit.steps.index{how_many:end}};
            circuit.steps.param={circuit.steps.param{1:how_many-1},[],circuit.steps.param{how_many:end}};
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
            circuit.steps.gates={circuit.steps.gates{1:how_many-1},'CNOT1',circuit.steps.gates{how_many:end}};
            circuit.steps.index={circuit.steps.index{1:how_many-1},index,circuit.steps.index{how_many:end}};
            circuit.steps.param={circuit.steps.param{1:how_many-1},[],circuit.steps.param{how_many:end}};
        elseif c_bits(1)==0
            circuit.steps.gates={circuit.steps.gates{1:how_many-1},'CNOT0',circuit.steps.gates{how_many:end}};
            circuit.steps.index={circuit.steps.index{1:how_many-1},index,circuit.steps.index{how_many:end}};
            circuit.steps.param={circuit.steps.param{1:how_many-1},[],circuit.steps.param{how_many:end}};
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

