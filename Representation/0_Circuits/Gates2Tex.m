function [ file_names ] = Gates2Tex( gate,elem_gates,comp_gates,file_name,modes,depth,sep_dist )
%GATES2TEX creates a latex circuit diagram of gate and in case it is a
%composite gate of it's construction -> can be used in conjunction with
%tex2pdf2png
%mode specifies the kind of gate that will be constructed
%       - 'gate' creates the gate
%       - 'expansion' creates the expansion of the gate
%       - 'no_sep' the structure does not get separated into multiple plots
%       - 'no_head' delete the header of the tex file
%       - gate=expansion' creates both next to each other (not yet implemented)
%               Can be done via ghost qubits for missing ancilla bits
%depth is how deep the hierarchies get unraveled
%sep_dist is how many gate columns are shown per file
if nargin==6
    sep_dist=85;
elseif nargin==5
    depth=0; %no expansion of the gate
    sep_dist=85;
elseif nargin==4
    modes='gate+expansion';
    depth=0;
    sep_dist=85;
end

curr_folder=pwd;
save_folder=fullfile(curr_folder,'Circuits');
if exist(save_folder)==0
    mkdir('Circuits');
end

file_names={};
%Top and Bottom elements
tex_start='\\documentclass[border={25pt 10pt 25pt 10pt}]{standalone} \n\\usepackage[braket, qm]{qcircuit} \n\\usepackage{amsmath} \n\\usepackage{fixltx2e,graphicx,mathpazo}\n\n\\begin{document}\n';
tex_end='\\end{document}';

%The Elements circuit replacement
tex_code='\\Qcircuit @C=.7em @R=.4em @!R {\n';

if length(strfind(modes,'gate'))>0 
    %% Create Gate representation
    sizler=gate.size;
    if any(strcmp(fields(gate),'anc_size'))
        anc_sizler=gate.anc_size;
    else
        anc_sizler=0;
    end
    
        
    [block indexes block2 indexes2]=Block_Creator(gate,1:sizler+anc_sizler,sizler,anc_sizler,0);
    block=Block_Corrector( block,indexes,block2,indexes2,[]);
    block=Block_Parallelizer(block);
    
    if length(strfind(modes,'gate='))==0
        blocker=cell(size(block,1),1);
        for i=1:size(block,1)
            blocker{i,1}=['\\lstick{\\ket{\\phi}_{' num2str(i) '}}'];
        end
        block=[blocker block];
    end
    
    block(1:end,:)=block(end:-1:1,:);
    if length(strfind(modes,'gate='))~=0
        block_unexpanded=block(1:end,1:end);
    end
    
    block_rows=cell(size(block,1),1);
    for i=1:size(block,1)
        block_rows{i}=strjoin(block(i,:),' & ');
        if i<size(block,1)
            block_rows{i}=[block_rows{i},' & \\qw \\\\ \n'];
        else
            block_rows{i}=[block_rows{i}, ' & \\qw  \n'];
        end
    end
    
    %Separate gate and expansion?
    if length(strfind(modes,'gate='))==0
        block_rows=strjoin(block_rows,' ');
        if length(strfind(modes,'no_head'))==0
            tex_current=[tex_start tex_code block_rows '}\n' tex_end];
        else
            tex_current=[tex_code block_rows '}\n'];
        end
        
        save_folder2=fullfile(save_folder,file_name);
        if exist(save_folder2)==0
            mkdir(save_folder2);
        end
        fileID=fopen([fullfile(save_folder2,file_name),'.tex'],'w');
        fprintf(fileID,tex_current);
        fclose(fileID);
        file_names={file_names{:} fullfile(save_folder2,file_name)};
    end
end

%% Create second graphic with the whole circuit of the previously substituted one
if length(strfind(modes,'expansion'))>0
    %We trust this size to be enough

    sizler=gate.size;
    if any(strcmp(fields(gate),'anc_size'))
        anc_sizler=gate.anc_size;
    else
        anc_sizler=0;
    end
    
    %[anc_sizler]=determine_anc_size(gate,depth,elem_gates,comp_gates)
    % Determines sub structure anc_bits
    block=CircuitBlocks2Tex(gate,1:(sizler+anc_sizler),0,sizler,anc_sizler,elem_gates,comp_gates,modes,depth);
    block=Block_Parallelizer(block);
    
    if length(strfind(modes,'gate='))~=0
        %Add unexpanded block and block together
        unex_size=size(block_unexpanded);
        expa_size=size(block);
        if unex_size(1)>expa_size(1)
            %Add qubits for unex
            cell_row=cell(unex_size(1)-expa_size(1),expa_size(2));
            [cell_row{1:end,1:end}]=deal('\\qw');
            block=[cell_row;block];
        elseif unex_size(1)<expa_size(1)
            cell_row=cell(expa_size(1)-unex_size(1),unex_size(2));
            [cell_row{1:end,1:end}]=deal('\\qw');
            block_unexpanded=[cell_row;block_unexpanded];
        end
        cell_col=cell(size(block_unexpanded,1),1);
        [cell_col{1:end,1:end}]=deal('\\qw');
        block_unexpanded=[block_unexpanded,cell_col];
        %Add equality in between
        i=ceil(size(block_unexpanded,1)/2);
        if mod(size(block,1),2)==1
            block_unexpanded{i,end}=['\\rstick{\\raisebox{-',num2str(0.5),'em}{=}} \\qw'];
        else            
            block_unexpanded{i,end}=['\\rstick{\\raisebox{-',num2str(2.5),'em}{=}} \\qw'];
        end
        [cell_col{1:end,1:end}]=deal('');
        block_unexpanded=[block_unexpanded,cell_col,cell_col];
        %Fuse both together
        block=[block_unexpanded,block];
    end
    
    if length(strfind(modes,'no_sep'))>0 %Do not separate the image
        block_init=cell(size(block,1),1);
        for i=1:sizler
            block_init{i,1}=['\\lstick{\\ket{\\phi}_{' num2str(i) '}}'];
        end
        for i=sizler+1:size(block,1)
            block_init{i,1}=['\\lstick{\\ket{0}}'];
        end
        block_init(1:end)=block_init(end:-1:1);
        %block(1:end)=block(end:-1:1);
        block=[block_init,block];
        block_rows=cell(size(block,1),1);
        for i=1:size(block,1)
            block_rows{i}=strjoin(block(i,:),' & ');
            if i<size(block,1)
                block_rows{i}=[block_rows{i},' & \\qw \\\\ \n'];
            else
                block_rows{i}=[block_rows{i}, ' & \\qw  \n'];
            end
        end
        block_rows=strjoin(block_rows,' ');
        
        if length(strfind(modes,'gate='))==0
            tex_current=[tex_start tex_code block_rows '}\n' tex_end];
            new_name=[file_name '_expanded'];
        else
            tex_current=[tex_code block_rows '}\n'];
            new_name=[file_name '_gate_expanded'];
        end
        fileID=fopen([new_name,'.tex'],'w');
        fprintf(fileID,tex_current);
        fclose(fileID);
        file_names={file_names{:} new_name};
            
    else %Default - separate the file into manageable chunks
        block_len=size(block,2);
        block2=block;
        for num=1:ceil(block_len/sep_dist)
            if num==1
                block_init=cell(size(block,1),1);
                for i=1:sizler
                    block_init{i,1}=['\\lstick{\\ket{\\phi}_{' num2str(i) '}}'];
                end
                for i=sizler+1:size(block,1)
                    block_init{i,1}=['\\lstick{\\ket{0}}'];
                end
                block_init(1:end)=block_init(end:-1:1);
            else %Dots
                block_init=cell(size(block,1),1);
                for i=1:size(block,1)
                    block_init{i,1}=['\\lstick{}'];
                end
                if mod(size(block,1),2)==0
                    block_init_end=['\\inputgroup{',num2str(size(block,1)/2),'}{',num2str(size(block,1)/2+1),'}{' num2str(0.9) 'em}{\\cdots}'];
                else
                    block_init_end=['\\inputgroup{',num2str(ceil(size(block,1)/2)),'}{',num2str(ceil(size(block,1)/2)),'}{.0em}{\\cdots}'];
                end
            end
            %block(1:end)=block(end:-1:1);
            if num<ceil(block_len/sep_dist)
                block=block2(:,sep_dist*(num-1)+(1:sep_dist));
            else
                block=block2(:,((sep_dist*(num-1)+1):end));
            end
            
            block=[block_init,block];
            block_rows=cell(size(block,1),1);
            for i=1:size(block,1)
                block_rows{i}=strjoin(block(i,:),' & ');
                if mod(size(block,1),2)==1
                    if i==ceil(size(block,1)/2) && num<ceil(block_len/sep_dist)
                        if i<size(block,1)
                            block_rows{i}=[block_rows{i},' & \\rstick{\\cdots} \\qw  \\\\ \n'];
                        else
                            block_rows{i}=[block_rows{i}, ' & \\rstick{\\cdots} \\qw   \n'];
                        end
                    else
                        if i<size(block,1)
                            block_rows{i}=[block_rows{i},' & \\qw \\\\ \n'];
                        else
                            block_rows{i}=[block_rows{i}, ' & \\qw  \n'];
                        end
                    end
                else
                    if num~=ceil(block_len/sep_dist)
                        if i==size(block,1)/2  && num<ceil(block_len/sep_dist)
                            if i<size(block,1)
                                block_rows{i}=[block_rows{i},' & \\rstick{\\raisebox{-',num2str(2),'em}{$\\cdots$}} \\qw  \\\\ \n'];
                            else
                                block_rows{i}=[block_rows{i}, ' & \\rstick{\\raisebox{-',num2str(2),'em}{$\\cdots$}} \\qw   \n'];
                            end
                        else
                            if i<size(block,1)
                                block_rows{i}=[block_rows{i},' & \\qw \\\\ \n'];
                            else
                                block_rows{i}=[block_rows{i}, ' & \\qw  \n'];
                            end
                        end
                    else
                        if i<size(block,1)
                            block_rows{i}=[block_rows{i},' & \\qw \\\\ \n'];
                        else
                            block_rows{i}=[block_rows{i}, ' & \\qw  \n'];
                        end
                    end
                end
            end
            block_rows=strjoin(block_rows,' ');
            
            if length(strfind(modes,'no_head'))==0
                if num~=1
                    tex_current=[tex_start tex_code block_rows block_init_end '}\n' tex_end];
                else
                    tex_current=[tex_start tex_code block_rows '}\n' tex_end];
                end
            else
                if num~=1
                    tex_current=[tex_code block_rows block_init_end '}\n'];
                else
                    tex_current=[tex_code block_rows '}\n'];
                end
            end
                
            %Check i
            save_folder2=fullfile(save_folder,file_name);
            if exist(save_folder2)==0
                mkdir(save_folder2);
            end
            if length(strfind(modes,'gate='))==0
                new_name=[fullfile(save_folder2,file_name) '_expanded_' num2str(depth) '_part_' num2str(num)];
            else
                new_name=[fullfile(save_folder2,file_name) '_gate_expanded_' num2str(depth) '_part_' num2str(num)];
            end
            
            fileID=fopen([new_name,'.tex'],'w');
            fprintf(fileID,tex_current);
            fclose(fileID);
            file_names={file_names{:} new_name};
        end
    end
end

