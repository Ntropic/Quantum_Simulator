function [ anc_size ] = determine_anc_size( gate,depth,elem_gates,comp_gates,varargin )
%DETERMINE_ANC_SIZE recursively determines the number of ancilla qbits within a
%quantum gate circuit 
if any(strcmp(fields(gate),'anc_size'))
    primary_anc_size=gate.anc_size;
else
    primary_anc_size=0;
end

if nargin>4
    for i=1:length(varargin)
        comp_gates=Comp_Gate_Merger(comp_gates,varargin{i});
    end
end

secondary_anc_size=[0];

if depth~=0
    %Go through the composite gates
    if any(strcmp(fields(gate),'steps'))
        if length(gate.steps.index)>0
            %There exist subfields! go through the
            for i=1:length(gate.steps.gates)
                sub_gate=Gate_by_Name(gate.steps.gates(i),elem_gates,comp_gates);
                [anc_sizle]=determine_anc_size(sub_gate,depth-1,elem_gates,comp_gates);
                secondary_anc_size=[secondary_anc_size,anc_sizle];
            end
        else
            secondary_anc_size=0;
        end
    else
        secondary_anc_size=0;    
    end
end
anc_size=primary_anc_size+max(secondary_anc_size);
end

