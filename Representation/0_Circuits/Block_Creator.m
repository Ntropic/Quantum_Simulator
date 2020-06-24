function [ block indexes block2 indexes2 ] = Block_Creator( gate,index_list,sizler,anc_sizler,sub_anc_length )
%BLOCK_CREATOR creates blocks from gates

    circs=gate.circuit;
    if isa(gate.circuit,'cell')
        if isnumeric(gate.circuit{1})
            elements=circs(1,2:length(circs));
            indexes=circs{1,1};
        else %Assume the indexes are 1:n
            elements=circs;
            indexes=[1:size(circs,2);ones(1,size(circs,2))];
        end
        if size(circs,1)>1
            %Alternative representation
            if isnumeric(gate.circuit{2,1})
                elements2=circs(2,2:size(circs,2));
                indexes2=circs{2,1};
            else %Assume the indexes are 1:n
                error(['Gate [' gate.names{1} '] has multiple circuit representations. Indexing required, and info needs to be in cell.']);
            end
        end
        if length(circs{1})==0
            dont_add=1;
        else
            dont_add=0;
        end
    else
        elements={circs};
        indexes=[1;1];
        if length(circs)==0
            dont_add=1;
        else
            dont_add=0;
        end
    end
    
    if dont_add==0
        if isa(elements,'cell')
            elements=cellfun(@(x) strrep(x,'\','\\'),elements,'UniformOutput',false);
        else
            elements={strrep(elements,'\','\\')};
        end    
        if exist('elements2')
             if isa(elements2,'cell')
                elements2=cellfun(@(x) strrep(x,'\','\\'),elements2,'UniformOutput',false);
            else
                elements2={strrep(elements2,'\','\\')};
            end
        end
        indexes(1,:)=index_list(indexes(1,:));
    
        rows=sizler+anc_sizler+sub_anc_length;
        columns=max(indexes(2,:));
        block=cell(rows,columns);
        [block{:,:}]=deal('\\qw');

        for i=1:length(elements)
            block(indexes(1,i),indexes(2,i))=elements(i);
            if exist('elements2')
                block2(indexes(1,i),indexes(2,i))=elements2(i);
            end
        end
        if exist('block2')==0
            block2=[];
            indexes2=[];
        end
    else
        rows=sizler+anc_sizler+sub_anc_length;
        columns=0;
        block=cell(rows,columns);
        indexes=[];
        
        block2=[];
        indexes2=[];
    end
end

