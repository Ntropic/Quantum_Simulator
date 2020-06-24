function [ block ] = CircuitBlocks2Tex( gate,index_list,sub_anc_length,sizler,anc_sizler,elem_gates,comp_gates,modes,depth,circuit_subs)
%CircuitBlocks2Tex creates blocks for a latex circuit diagram of a gates circuit 
%expansion -> usually used via Gates2Tex
%
%up to depth
if nargin==8
    depth=0; %no expansion of the gate
    circuit_subs=[];
elseif nargin==9
    circuit_subs=[];
end
end_reached=0;

%No expansion of this gate
block=cell(0);
block2=cell(0);


if depth~=0 %Check if we have reached the end of the circuit recursion tree
    %% Add layers
    if any(strcmp(fields(gate),'steps'))
        if any(strcmp(fields(gate),'circuit_subs'))
            circuit_subs=gate.circuit_subs;
        end
        if length(gate.steps.index)>0
            index_lists=gate.steps.index;
            gates_list=gate.steps.gates;
            for i=1:length(gates_list)
                %Get the gate
                found_gate=Gate_by_Name(gates_list(i),elem_gates,comp_gates);
                
                if any(strcmp(fields(found_gate),'anc_size'))
                    curr_anc_sizler=found_gate.anc_size;
                else
                    curr_anc_sizler=0;
                end
                if depth~=1
                    sub_anc_len=sub_anc_length+curr_anc_sizler;
                else
                    sub_anc_len=sub_anc_length;
                end
                if sub_anc_len>0
                    index_list=[index_list sizler+anc_sizler+(1:sub_anc_len)];
                end
                found_index_list=index_list(index_lists{i});
                %Recursively get the gates inserted
                if length(circuit_subs)==0        
                    [found_block]=CircuitBlocks2Tex(found_gate,found_index_list,sub_anc_len,sizler,anc_sizler,elem_gates,comp_gates,modes,depth-1);
                else
                    if length(circuit_subs{i})==0
                        [found_block]=CircuitBlocks2Tex(found_gate,found_index_list,sub_anc_len,sizler,anc_sizler,elem_gates,comp_gates,modes,depth-1);
                    else
                        [found_block]=CircuitBlocks2Tex(found_gate,found_index_list,sub_anc_len,sizler,anc_sizler,elem_gates,comp_gates,modes,depth-1,circuit_subs{i});
                    end
                end
                

                %Add found block to the current list of blocks
                old_block=block;
                size_block=max([size(old_block,1),size(found_block,1)]);
                size_block(2)=size(old_block,2)+size(found_block,2);
                block=cell(size_block);
                block(:)={'\\qw'};
                if size_block(1)>0 && size_block(2)>0
                    if size(old_block,1)>0 && size(old_block,2)>0
                        block(size_block(1)+1-size(old_block,1):size_block(1),1:size(old_block,2))=old_block;
                        if size(found_block,1)>0 && size(found_block,2)>0
                            block(size_block(1)+1-size(found_block,1):size_block(1),size(old_block,2)+(1:size(found_block,2)))=found_block;
                        end 
                    else
                        block(size_block(1)+1-size(found_block,1):size_block(1),1:size(found_block,2))=found_block;
                    end
                end

            end
        else
            end_reached=0;
            %Don't add to list
        end
    else
        end_reached=1;
        %Add to list -> this is an elementary gate
    end
end
if depth==0 || end_reached==1
    %% Create Gate representation
    %The ancilliary qubit elements have to go to the next ancilla qubit
    %in the list
    [block indexes block2 indexes2]=Block_Creator(gate,index_list,sizler,anc_sizler,sub_anc_length);
    if length(circuit_subs)>0
        [block]=Block_Corrector( block,indexes,block2,indexes2,circuit_subs);
    else
        block=Block_Corrector( block,indexes,block2,indexes2,[]);
    end
    block(1:end,:)=block(end:-1:1,:);
end
end

