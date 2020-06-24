function [ filename ] = Save_Gate_list( filename,gates_name,gates )
%SAVE_GATE_LIST Saves a list of gates
eval([gates_name '=gates;'])
filename=[filename '.mat'];
save(filename,gates_name);
end

