function [ index ] = Gate_Index_by_Name( gate_name,elem_gates,comp_gates,varargin )
%GATE_INDEX_BY_NAME finds a gates index in lists of gates (Case insensitive)

index=zeros(1,nargin-1);

found=0;
for j=1:length(elem_gates)
    obj=strcmpi(elem_gates(j).names,gate_name);
    if iscell(obj)
        obj=cell2mat(obj);
    end
    if any(obj)
        %Found the gate
        found=1;
        index(1)=j;
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
                index(2)=j;
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
                    index(2)=i;
                    break;
                end
            end
        end
    end
end

if nargin>3
    k=3;
    while k-2<length(varargin) && found==0
        new_comp_gates=varargin{k-2};
        if iscell(new_comp_gates)==0
            for j=1:length(new_comp_gates)
                obj=strcmpi(new_comp_gates(j).names,gate_name);
                if iscell(obj)
                    obj=cell2mat(obj);
                end
                if any(obj)
                    %Found the gate
                    found=1;
                    index(k)=j;
                    break;
                end
            end
        else
            for i=1:length(new_comp_gates)
                gates=new_comp_gates{i};
                for j=1:length(gates)
                    obj=strcmpi(gates(j).names,gate_name);
                    if iscell(obj)
                        obj=cell2mat(obj);
                    end
                    if any(obj)
                        %Found the gate
                        found=1;
                        index(k)=i;
                        break;
                    end
                end
            end
        end
        k=k+1;
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

