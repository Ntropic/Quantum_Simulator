function [ gate ] = Create_Empty_Comp_Gate( name,size,anc_size,circ  )
%Create_Comp_Gate creates a new gate

gate=struct('names',name);
gate.size=size;
gate.anc_size=anc_size;
gate.step_num=0;

gate.steps.index={};
gate.steps.gates={};
gate.steps.param={};

if nargin>3
    if length(circ)>0
        gate.circuit=circ;
    end
end
gate.circuit_subs={};
end

