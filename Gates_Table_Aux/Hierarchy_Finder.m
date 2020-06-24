function [ other_gates ] = Hierarchy_Finder( k,elem_gates,other_gates,counter )
%HIERARCHY FINDER determines the hierarchy order of gates elem_gates being
%1 and gates consisting only of elem_gates being a 2 a gate always having
%order of the highest gate they are composed of +1. Hope you didn't create
%loops in your gate definitions
%Exits when the order of all elements that gate is composed of have been 
%determined and therefor gates order is also known
if counter>length(other_gates)+length(elem_gates)+2;
    error('Gates must have a loop! This is outrageous!');
end

gates=other_gates(k);
if length(gates.hierarchy)==0
    h_i=zeros(1,length(gates.steps.gate_ind));
    for i=1:length(gates.steps.gate_ind)
        if gates.steps.gate_ind{i}(1)~=0
            h_i(i)=1;
        else
            curr_gate=other_gates(gates.steps.gate_ind{i}(2));
            if length(curr_gate.hierarchy)>0
                h_i(i)=curr_gate.hierarchy;
            else
                other_gates=Hierarchy_Finder(gates.steps.gate_ind{i}(2),elem_gates,other_gates,counter+1);
                h_i(i)=other_gates(gates.steps.gate_ind{i}(2)).hierarchy;
            end
        end
    end
    if length(gates.steps.gate_ind)==0
        h_i=1;
    end
    other_gates(k).hierarchy=max(h_i)+1;
end
end

