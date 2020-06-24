function [ gate_data ] = Gate_by_Name( gate_name,elem_gates,comp_gates,varargin )
%GATE_BY_NAME finds a gate in lists of gates (Case insensitive)

if nargin>3
    for i=1:length(varargin)
        comp_gates=Comp_Gate_Merger(comp_gates,varargin{i});
    end
end

found=0;
for j=1:length(elem_gates)
    obj=strcmpi(elem_gates(j).names,gate_name);
    if iscell(obj)
        obj=cell2mat(obj);
    end
    if any(obj)
        %Found the gate
        found=1;
        gate_data=elem_gates(j);
        break;
    end
end
if found==0
    if iscell(comp_gates)==0
        for j=1:length(comp_gates)
            obj=strcmpi(comp_gates(j).names,gate_name);
            if iscell(obj)
                obj=cell2mat(obj);
            end
            if any(obj)
                %Found the gate
                found=1;
                gate_data=comp_gates(j);
                break;
            end
        end
    else
        for i=1:length(comp_gates)
            gates=comp_gates{i};
            for j=1:length(gates)
                obj=strcmpi(gates(j).names,gate_name);
                if iscell(obj)
                    obj=cell2mat(obj);
                end
                if any(obj)
                    %Found the gate
                    found=1;
                    gate_data=gates(j);
                    break;
                end
            end
        end
end

if found==0
    if isa(gate_name,'char')
        gate_names=gate_name;
    else
        gate_names=gate_name{1};
    end
    error_txt=strcat('Gate "',gate_names,'" was not found.');
    error(error_txt)
end


end

