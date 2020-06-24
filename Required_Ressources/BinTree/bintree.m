classdef bintree    
    %Depth specifies the number of layers
    %Status specifies the (bit) state of the nodes within the tree
    %Pointer specifies the current element (serves as an index within the status list)
    %Path specifies the path through the tree to the pointer
    %Number specifies the index of the pointer within its respective layer
    %Layer specifies the Layer of the path of the Pointer at which we look
    
    properties
        Depth;
        Status;
        Pointer;
        Path;
        Number;
        Layer;
        
    end
    
    methods
        %% Constructor -----------------------------------------------------------------------------
        function obj=bintree(depth,n)
            if nargin==1
                n=1;
            end
            %Tree is filled with ones up to point n
            %And is filled with zeros from n to 2^m
            status=[ones(n-1,1);zeros(2^depth-n+1,1)];
            log2n=(n-mod(n,2))/2;
            if depth>1
                for i=depth-1:-1:0
                    status2=zeros(2^i,1);  %length 2^i
                    status2(1:round(log2n))=ones(round(log2n),1);
                    log2n=(n+2^(depth-i)-1-mod(n+2^(depth-i)-1,2^(depth-i+1)))/2^(depth-i+1);
                    status=[status2;status];
                end
            end
            obj.Status=status;
            obj.Depth=depth;
            
            %Create Path
            obj.Pointer=2^depth+n-1;
            obj.Path=dec2bin(n-1,depth);
            obj.Number=n;
            obj.Layer=depth;
        end
        
        %Find Element position in tree
        function path=curr_path(obj)
            path=obj.Path;
        end
        
        %Find Element pointer
        function pointer=curr_pointer(obj)
            pointer=obj.Pointer;
        end
        
        %Find Element Number (in layer)
        function number=curr_number(obj)
            number=obj.Number;
        end
        
        %Find Element Number (in layer)
        function layer=curr_layer(obj)
            layer=obj.Layer;
        end
        
        function [obj lay_state]=layer_up(obj)
            layer=obj.Layer-1;
            obj.Layer=layer;
            lay_state=obj.Status(2^layer+bin2dec(obj.Path(1:layer)));
        end
        
        function obj=change_path(obj,path)
            obj.Path=path;
            [pos layer number]=path2pointers(obj,path);
            obj.Pointer=pos;
            obj.Layer=layer;
            obj.Number=number;
        end
        
        function [perm]=all_num_algorithm(obj)
            %Go through all numbers of the objectand create a permutation
            %of the numbers n to 2^depth
            
            %First add primary number to list
            perm=obj.Number;
            obj.Status(obj.Pointer)=1;
            while obj.Status(1)==0 || obj.Layer~=0
                if obj.Status(path2pointers(obj,obj.Path(1:obj.Layer)))==1
                    obj=layer_up(obj);
                elseif obj.Status(path2pointers(obj,obj.Path(1:obj.Layer)))==0
                    obj.Status(path2pointers(obj,obj.Path(1:obj.Layer)))=1;
                    new_path=obj.Path;
                    if new_path(obj.Layer+1)=='0'
                        new_path(obj.Layer+1)='1';
                    else
                        new_path(obj.Layer+1)='0';
                    end
                    obj=change_path(obj,new_path);
                    perm=[perm obj.Number];
                    obj.Status(path2pointers(obj,obj.Path(1:obj.Layer)))=1;
                    obj=layer_up(obj);
                end 
            end
        end
        
        function [perm]=all_num_algorithm_gif(obj,gif_name)
            %Go through all numbers of the objectand create a permutation
            %of the numbers n to 2^depth
            
            %First add primary number to list
            perm=obj.Number;
            
            %Make first frame
            [A,h]=make_image_from_tree(obj,1);
            [imind,cm] = rgb2ind(A,256); 
            imwrite(imind,cm,gif_name,'gif','Loopcount',inf); 
            
            obj.Status(obj.Pointer)=1;
            
            %Add frame
            [A,h]=make_image_from_tree(obj,1,h);
            [imind,cm] = rgb2ind(A,256); 
            imwrite(imind,cm,gif_name,'gif','WriteMode','append'); 
            
            while obj.Status(1)==0 || obj.Layer~=0
                if obj.Status(path2pointers(obj,obj.Path(1:obj.Layer)))==1
                    obj=layer_up(obj);
                    
                    [A,h]=make_image_from_tree(obj,1,h);
                    [imind,cm] = rgb2ind(A,256); 
                    imwrite(imind,cm,gif_name,'gif','WriteMode','append'); 
                elseif obj.Status(path2pointers(obj,obj.Path(1:obj.Layer)))==0
                    obj.Status(path2pointers(obj,obj.Path(1:obj.Layer)))=1;
                    
                    A=make_image_from_tree(obj,1,h);
                    [imind,cm] = rgb2ind(A,256); 
                    imwrite(imind,cm,gif_name,'gif','WriteMode','append');
                    
                    new_path=obj.Path;
                    if new_path(obj.Layer+1)=='0'
                        new_path(obj.Layer+1)='1';
                    else
                        new_path(obj.Layer+1)='0';
                    end
                    obj=change_path(obj,new_path);
                    perm=[perm obj.Number];
                        
                    A=make_image_from_tree(obj,1,h);
                    [imind,cm] = rgb2ind(A,256); 
                    imwrite(imind,cm,gif_name,'gif','WriteMode','append'); 
                    
                    obj.Status(path2pointers(obj,obj.Path(1:obj.Layer)))=1;
                    
                    A=make_image_from_tree(obj,1,h);
                    [imind,cm] = rgb2ind(A,256); 
                    imwrite(imind,cm,gif_name,'gif','WriteMode','append'); 
                    
                    obj=layer_up(obj);
                    
                    A=make_image_from_tree(obj,1,h);
                    [imind,cm] = rgb2ind(A,256); 
                    imwrite(imind,cm,gif_name,'gif','WriteMode','append'); 
                end 
            end
        end
        
        function [perm]=all_num_algorithm_tex(obj)
            %Go through all numbers of the objectand create a permutation
            %of the numbers n to 2^depth
            num=1;
            %First add primary number to list
            perm=obj.Number;
            
            %Make first frame
            make_tex_from_tree(obj,num);
            num=num+1;
            
            obj.Status(obj.Pointer)=1;
            
            %Add frame
            make_tex_from_tree(obj,num);
            num=num+1;
            
            while obj.Status(1)==0 || obj.Layer~=0
                if obj.Status(path2pointers(obj,obj.Path(1:obj.Layer)))==1
                    obj=layer_up(obj);
                    
                    make_tex_from_tree(obj,num);
                    num=num+1;
                elseif obj.Status(path2pointers(obj,obj.Path(1:obj.Layer)))==0
                    obj.Status(path2pointers(obj,obj.Path(1:obj.Layer)))=1;
                    
                    make_tex_from_tree(obj,num);
                    num=num+1;
                    
                    new_path=obj.Path;
                    if new_path(obj.Layer+1)=='0'
                        new_path(obj.Layer+1)='1';
                    else
                        new_path(obj.Layer+1)='0';
                    end
                    obj=change_path(obj,new_path);
                    perm=[perm obj.Number];
                        
                    make_tex_from_tree(obj,num);
                    num=num+1;
                    
                    obj.Status(path2pointers(obj,obj.Path(1:obj.Layer)))=1;
                    
                    make_tex_from_tree(obj,num);
                    num=num+1;
                    
                    obj=layer_up(obj);
                    
                    make_tex_from_tree(obj,num);
                    num=num+1;
                end 
            end
        end
        
        
        %Find arbitrary position in tree pointer2path
        function [path layer number]=pointer2path(obj,pos)
            %Element is in which layer
            layer=floor(log2(pos));
            %i_n-th element in layer
            number=pos-2.^layer+1;
            path=dec2bin(number-1,layer(1));
            if length(unique(layer))>1
                error('Inputs have to be of same layer');
            end
        end
        
        %Path to Pointer
        function [pos layer number]=path2pointers(obj,path)
            %Element is in which layer
            layer=length(path);
            %i_n-th element in layer
            if length(path)>0
                if isa(path,'numeric')
                    number=bin2dec(num2str(path))+1;
                else
                    number=bin2dec(path)+1;
                end
            else
                number=1;
            end
            pos=2^layer-1+number;
        end
        
        function []=make_tex_from_tree(obj,num)
            file_name=['tree_' num2str(num)];
            plot_obj=plot_tree(obj);
            fileID=fopen([file_name,'.tex'],'w');
            fprintf(fileID,plot_obj);
            fclose(fileID);
        end
        
        function [A,h]=make_image_from_tree(obj,plotting,h)
            if nargin==1
                plotting=1;
                h=[];
            elseif nargin==2
                h=[];
            end
            file_name='tree';
            plot_obj=plot_tree(obj);
            fileID=fopen([file_name,'.tex'],'w');
            fprintf(fileID,plot_obj);
            fclose(fileID);
            if length(h)==0
                [A,h]=tex2pdf2preview(file_name,plotting,1000);
            else
                [A,h]=tex2pdf2preview(file_name,plotting,1000,h);
            end
        end
        
        function plot_obj=plot_tree(obj)
            %start text for latex document
            text='\\documentclass{standalone}\n\\usepackage[utf8]{inputenc}\n\\usepackage{forest}\n\\usepackage{tikz}\n\\usetikzlibrary{shapes.geometric}\n\n';
            text=[text '\\tikzset{\n\traute/.style={regular polygon,regular polygon sides=4,inner sep=-0.5pt,shape border rotate=45,minimum size=14pt},\n\tsquare/.style={regular polygon,regular polygon sides=4,inner sep=.4pt,minimum size=19.798pt}\n}\n'];
            text=[text '\n\\begin{document}\n\\begin{forest}\nanchors/.style={anchor=#1,child anchor=#1,parent anchor=#1},\n\tfor tree={\n\t\tchild anchor=north,\n\t\tparent anchor=south,\n\t\tedge path={\\noexpand\\path[\\forestoption{edge}]\n'];
            text=[text '\t\t\t(\\forestOve{\\forestove{@parent}}{name}.parent anchor)\n\t\t\t-- ([shift={(0,-10pt)}]\n\t\t\t\\forestOve{\\forestove{@parent}}{name} -| \\forestove{name})\n'];
            text=[text '\t\t\t-- (\\forestove{name}.child anchor)\n\t\t\t\\forestoption{edge label};\n\t\t}\n\t},\n%%The Tree\n'];
        
            %Create forest/tree
            depth=obj.Depth;
            pointer=obj.Pointer;
            status=obj.Status;
            layer=obj.Layer;
            been_there=zeros(size(status));
            
            %What path?!
            [path layer2 number]=pointer2path(obj,pointer);
            
            curr_layer=depth;
            layer_list=status(2^depth:end)';
            tabs=repmat('\t',1,depth);
            tree_text=cellfun(@(x) horzcat(tabs,'[',x,',raute,draw=white'),strsplit(num2str(layer_list),' '),'UniformOutput',false);
            tree_text(1:2:end)=cellfun(@(x) horzcat(x,',edge label={node[pos=0,left=-0.1cm]{\\scriptsize \\color{white}{',path(end),'}}}]\n'),tree_text(1:2:end),'UniformOutput',false);
            tree_text(2:2:end)=cellfun(@(x) horzcat(x,',edge label={node[pos=0,right=-0.1cm]{\\scriptsize \\color{white}{',path(end),'}}}]\n'),tree_text(2:2:end),'UniformOutput',false);
            if layer2==depth
                if layer~=depth
                    if path(end)=='0'
                        tree_text{number}=horzcat(tabs,'[',num2str(layer_list(number)),',raute,draw=white,edge=red,edge label={node[pos=0,left=-0.1cm]{\\scriptsize ',path(end),'}}]\n');
                    else
                        tree_text{number}=horzcat(tabs,'[',num2str(layer_list(number)),',raute,draw=white,edge=red,edge label={node[pos=0,right=-0.1cm]{\\scriptsize ',path(end),'}}]\n');
                    end
                else   
                    if path(end)=='0'
                        tree_text{number}=horzcat(tabs,'[',num2str(layer_list(number)),',raute,draw,edge=red,edge label={node[pos=0,left=-0.1cm]{\\scriptsize ',path(end),'}}]\n');
                    else
                        tree_text{number}=horzcat(tabs,'[',num2str(layer_list(number)),',raute,draw,edge=red,edge label={node[pos=0,right=-0.1cm]{\\scriptsize ',path(end),'}}]\n');
                    end
                end
            end
            
            for i=depth-1:-1:0
                tabs=repmat('\t',1,i);
                tree_text_old=reshape(tree_text,[2,length(tree_text)/2])';
                tree_text_comb=strcat(tree_text_old(:,1),tree_text_old(:,2));
                layer_list=status(2^i:2^(i+1)-1)';
                tree_text_new=cellfun(@(x) horzcat(tabs,'[',x,',raute,draw=white,'),strsplit(num2str(layer_list),' '),'UniformOutput',false);
                if layer2>=i && i>0
                    path2=path(1:i);
                    number2=bin2dec(path2)+1;
                    if layer~=i
                        if path2(end)=='0'
                            tree_text_new{number2}=horzcat(tabs,'[',num2str(layer_list(number2)),',raute,draw=white,edge=red,edge label={node[pos=0,left=-0.1cm]{\\scriptsize ',path2(end),'}}');
                        else
                            tree_text_new{number2}=horzcat(tabs,'[',num2str(layer_list(number2)),',raute,draw=white,edge=red,edge label={node[pos=0,right=-0.1cm]{\\scriptsize ',path2(end),'}}');
                        end
                    else
                        if path2(end)=='0'
                            tree_text_new{number2}=horzcat(tabs,'[',num2str(layer_list(number2)),',raute,draw,edge=red,edge label={node[pos=0,left=-0.1cm]{\\scriptsize ',path2(end),'}}');
                        else
                            tree_text_new{number2}=horzcat(tabs,'[',num2str(layer_list(number2)),',raute,draw,edge=red,edge label={node[pos=0,right=-0.1cm]{\\scriptsize ',path2(end),'}}');
                        end
                    end
                elseif layer==i && i==0
                    tree_text_new{1}=horzcat(tabs,'[',num2str(layer_list(1)),',raute,draw');
                end
                tree_text=cellfun(@(x,x2) horzcat(x2,'\n',x,tabs,']\n'),tree_text_comb,tree_text_new','UniformOutput',false);
            end
            plot_obj=[text,tree_text{1}];
            plot_obj(end-1:end)=[];
            %End text
            plot_obj=[plot_obj '\n\\end{forest}\n\\end{document}'];
        end
    end
    
end

