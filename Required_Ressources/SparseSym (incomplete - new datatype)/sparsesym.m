classdef sparsesym < matlab.mixin.CustomDisplay
    %SPARSESYM creates sparse symbolic matrices.
    %S=sparse(M) transforms a symbolic matrix M to sparse form by squeezing
    %out any zero elements.
    %
    %S=sparse(i,j) creates an empty sparse matrix of size i*j
    
    
    properties
        Index       %indexes of the sparse elements (for nth Element A_ij of matrix A Index(n,:)=[i,j,Symtypeindex]
        Symtypes    %The different possible elements (association is saved within the third column of Index)
        SymtypesChar
        Symsubs     %Possible Substiutions for Representations (string)
        Size        %Size of the matrix 
    end
    
    methods
%% Constructor -----------------------------------------------------------------------------

        function obj=sparsesym(varargin) 
            if nargin==1
                M=varargin{1};
                if isa(M,'numeric')
                    M=sym(M);
                end
                if isa(M,'sparsesym')
                    obj=M; 
                elseif isa(M,'sym') && isa(M,'sparsesym')==0
                    dim=size(M);
                    obj.Size=dim;
                    index=[];
                    symtypes=[];
                    symtypeschar=[];
                    symsubs=[];

                    M=M(:);
                    [M_sort,ind]=sort(M);
                    [M_unique,I1,I2]=unique(M_sort);
                    charM=char(M_unique);
                    if length(M_unique)>1
                        charM=charM(10:end-3);
                        charM=strsplit(charM,'], [');
                    else
                        charM={charM};
                    end
                    if any(strcmp(charM,'0')) %Get rid of '0' from lists!
                        i=find(strcmp(charM,'0'));
                        M_unique(i)=[];
                        charM(i)=[];
                        ind(find(I2==i))=[];
                        I2(find(I2==i))=[];
                        I2(find(I2>i))=I2(find(I2>i))-1;
                    end
                    obj.SymtypesChar=charM(:);
                    obj.Symsubs=charM(:);
                    obj.Symtypes=M_unique(:);

                    index=zeros(length(ind),3);
                    index(:,3)=I2;
                    [I,J]=ind2sub(dim,ind);
                    index(:,1)=I;
                    index(:,2)=J;
                    index=sortrows(index,[2 1]);
                    obj.Index=index;
                else
                    error('Input must be a symbolic or numeric matrix');
                end
            elseif nargin==2
                obj.Size=[varargin{1},varargin{2}];
            end
        end        
        
        function obj=subsasgn(obj,s,M) 
            switch s(1).type
                case '()'
                    dim=obj.Size;
                    ind=s.subs;
                    if length(ind)==2
                        if any(ind{1}>dim(1)) 
                            dim(1)=max(ind{1});
                        end
                        if any(ind{2}>dim(2))
                            dim(2)=max(ind{2});
                        end
                        indexes(:,1:2)=[repmat(reshape(ind{1},length(ind{1}),1),length(ind{2}),1),kron(reshape(ind{2},length(ind{2}),1),ones(length(ind{1}),1))];
                        if size(M)==[length(ind{1}) length(ind{2})]
                            corrsize=1;
                        else
                            corrsize=0;
                        end
                    elseif length(ind)==1
                        if any(ind{1}>dim(1)*dim(2))
                            error('Subscript assignment out of range.')
                        end
                        [I,J]=ind2sub(dim,ind{1});
                        indexes(:,1:2)=[I',J'];
                        if size(M)==[length(ind{1}) 1]
                            corrsize=1;
                        elseif size(M)==[1 length(ind{1})]
                            corrsize=1;
                        else
                            corrsize=0;
                        end
                    else
                        error('Not a valid index assignment.')
                    end
                    if corrsize==1
                        %Replace and add entries
                        index_old=obj.Index;
                        [~,i_change,i_B]=intersect(index_old(:,1:2),indexes,'rows');
                        %Is additional content symbolic
                        if isa(M,'numeric')
                            M=sym(M);
                        end
                        %Overwrite and add elements
                        %First delete elements that should be overwritten
                        %and rename indexes
                        index_old(i_change,:)=[];
                        obj.Index=index_old;
                        %Add new Elements and in case it's needed add new
                        %symtypes
                        %Transform ind to index_new
                        M=sparsesym(M);
                        indexes_new=indexes(sub2ind(M.size,M.Index(:,1),M.Index(:,2)),:);
                        indexes_new(:,3)=M.Index(:,3)+length(obj.Symtypes);
                        [obj,index_change_add]=AddSymtypes(obj,M.Symtypes);
                        [obj,index_change_unify]=UnifyDoubleSymtypes(obj);
                        indexes_new(:,3)=index_change_unify(indexes_new(:,3));
                        indexes=[obj.Index;indexes_new];
                        indexes=sortrows(indexes,[2 1 3]);
                        obj.Index=indexes;
                        obj.Size=dim;
                        obj=RemoveUnreferenced(obj);
                    else
                        error('Subscripted assignment dimension mismatch.')
                    end
                case '.'
                    error('Struct contents reference from a non-struct array object.')
                case '{}'
                    error('Cell contents reference from a non-cell array object.')
                otherwise
                    error('Not a valid indexing expression')
            end
        end
        
        function subobj=subsref(obj,s)
            switch s(1).type
                case '()'
                    if length(s) == 1
                        % Implement obj(indices)
                        ind=s.subs;
                        if length(ind)==2
                            dim=obj.Size;
                            if isa(ind{1},'char')
                                if strcmp(ind{1},':')
                                    ind{1}=1:prod(dim);
                                end
                            end
                            if isa(ind{2},'char')
                                if strcmp(ind{2},':')
                                    ind{2}=1:prod(dim);
                                end
                            end
                            subobj=sparsesym;
                            index_old=obj.Index;
                            index_new(:,1:2)=[repmat(reshape(ind{1},length(ind{1}),1),length(ind{2}),1),kron(reshape(ind{2},length(ind{2}),1),ones(length(ind{1}),1))];
                            subobj.Size=[length(ind{1}),length(ind{2})];
                            %Determine jnon trivial elements
                            [~,i_A,i_B] = intersect(index_old(:,1:2),index_new,'rows');
                            [I2,J2]=ind2sub(subobj.Size,i_B);
                            subobj.Index(:,1:2)=[I2,J2];
                            symtypechar=obj.SymtypesChar(index_old(i_A,3));
                            [symtypeunique,IA,IC]=unique(symtypechar);
                            symtype=obj.Symtypes(index_old(i_A,3));
                            symsubs=obj.Symsubs(index_old(i_A,3));
                            symtype=symtype(IA);
                            symsubs=symsubs(IA);
                            subobj.Index(:,3)=IC;
                            subobj.Symtypes=symtype;
                            subobj.SymtypesChar=symtypechar;
                            subobj.Symsubs=symsubs;
                        elseif length(ind)==1
                            dim=obj.Size;
                            if isa(ind{1},'char')
                                if strcmp(ind{1},':')
                                    ind{1}=1:prod(dim);
                                end
                            end
                            subobj=sparsesym;
                            index_old=obj.Index;
                            [I,J]=ind2sub(dim,ind{1});
                            index_new(:,1:2)=[I',J'];
                            subobj.Size=[length(ind{1}),1];
                            %Determine non trivial elements
                            [~,i_A,i_B] = intersect(index_old(:,1:2),index_new,'rows');
                            [I2,J2]=ind2sub(subobj.Size,i_B);
                            subobj.Index(:,1:2)=[I2,J2];
                            symtypechar=obj.SymtypesChar(index_old(i_A,3));
                            [symtypeunique,IA,IC]=unique(symtypechar);
                            symtype=obj.Symtypes(index_old(i_A,3));
                            symsubs=obj.Symsubs(index_old(i_A,3));
                            symtype=symtype(IA);
                            symsubs=symsubs(IA);
                            subobj.Index(:,3)=IC;
                            subobj.Symtypes=symtype;
                            subobj.SymtypesChar=symtypechar;
                            subobj.Symsubs=symsubs;
                        else
                            error('Tensors of order >2 are not (yet) supported for sparse symbolic matrices.')
                        end
                    else
                        error('Not a valid indexing expression.')
                    end
                case '.'
                    error('Struct contents reference from a non-struct array object.')
                case '{}'
                    error('Cell contents reference from a non-cell array object.')
                otherwise
                    error('Not a valid indexing expression.')
            end
        end
        
    end
    
    methods(Access = protected)
%% Create Output of Elements -------------------------------------------------------------

        function header = getHeader(obj)
            header =[];
            index=obj.Index();
            if length(index)>=1
                index=sortrows(index,[2 1]);
                symsubs=obj.Symsubs();
                dim=obj.Size();
                s_cell=cell(3,size(index,1));
                s_cell(1:2,:)=num2cell(index(:,1:2))';
                max_len=length(num2str(max(max(index(:,1:2)))));
                s_cell(1:2,:)=cellfun(@(x) sprintf('%*d',max_len,x),s_cell(1:2,:),'UniformOutput',false);
                s_cell(3,:)=symsubs(index(:,3));
                header=sprintf('\t(%s,%s) \t  %s\n',s_cell{:});
            else
                header=sprintf('\tEmpty SparseSym');
            end
        end
        
        function propertygroupelem = getPropertyGroups(obj)
            propgrp='';
            propertygroupelem=matlab.mixin.util.PropertyGroup(propgrp);
        end
   
        function footer = getFooter(obj)
            footer =sprintf('\b\b');
        end
        
    end
    
    methods 
        
%% GET Property Access Methods -------------------------------------------------------------
        
        function index=get.Index(obj)
            index=obj.Index;
        end
        
        function symtypes=get.Symtypes(obj)
            symtypes=obj.Symtypes;
        end
        
        function symsubs=get.Symsubs(obj)
            symsubs=obj.Symsubs;
        end
        
        function dim=get.Size(obj)
            dim=obj.Size;
        end
    end
    
    methods(Access = protected)
       
        %% Useful functions ----------------------------------------------------------------
        function subobj=RemoveSymtypes(obj,todelete)
            subobj=obj;
            if isa(todelete,'numeric')
                %Remove indexes
                subobj.Symtypes(todelete)=[];
                subobj.SymtypesChar(todelete)=[];
                subobj.Symsubs(todelete)=[];
            elseif isa(todelete,'char')
                todelete=strcmp(subobj.SymtypesChar,todelete);
                subobj.Symtypes(todelete)=[];
                subobj.SymtypesChar(todelete)=[];
                subobj.Symsubs(todelete)=[];
            elseif isa(todelete,'sym') && isa(todelete,'sparsesym')==0
                if length(todelete)>1
                    todelete=todelete(:);
                    todelete=char(todelete);
                    todelete=todelete(10:end-3);
                    todelete=strsplit(todelete,'], [');
                else
                    todelete=char(todelete);
                end
                [~,i_a,todelete]=intersect(todelete,subobj.SymtypesChar);
                subobj.SymtypesChar(todelete)=[];
                subobj.Symtypes(todelete)=[];
                subobj.Symsubs(todelete)=[];
            elseif isa(todelete,'cell')
                if isa(todelete{1},'sym') && isa(todelete,'sparsesym')==0
                    todelete=cellfun(@char,todelete,'UniformOutput',false);
                end
                if isa(todelete{1},'char')
                    [~,i_a,todelete]=intersect(todelete,subobj.SymtypesChar);
                    subobj.SymtypesChar(todelete)=[];
                    subobj.Symtypes(todelete)=[];
                    subobj.Symsubs(todelete)=[];
                else
                    error('Symtypes that shall be removed need to be specified via indexes, symbolic or character input or via a cell of character or symbolic inputs.')
                end
            else
                error('Symtypes that shall be removed need to be specified via indexes, symbolic or character input or via a cell of character or symbolic inputs.')
            end
        end

        function [obj,index_change]=UnifyDoubleSymtypes(obj)
            %Determine uniques
            symtypechar=obj.SymtypesChar;
            [symtypeunique,IA,IB]=unique(symtypechar);
            
            obj.Symtypes=obj.Symtypes(IA);
            obj.SymtypesChar=obj.SymtypesChar(IA);
            obj.Symsubs=obj.Symsubs(IA);
            
            %Now change the indexes of the remaining symtypes  
            indexes=obj.Index;
            indexes2=indexes;
            for i=1:length(IB)
                indexes2(find(indexes(:,3)==i),3)=IB(i);
            end
            indexes2=sortrows(indexes2,[2,1,3]);
            obj.Index=indexes2;
            index_change=IB;
        end
        
        function obj=RemoveUnreferenced(obj)
            %Determine referenced list
            indexes=obj.Index;   
            ref_list=unique(indexes(:,3));
            
            obj.SymtypesChar=obj.SymtypesChar(ref_list);
            obj.Symtypes=obj.Symtypes(ref_list);
            obj.Symsubs=obj.Symsubs(ref_list);
            
            indexes2=indexes;
            for i=1:length(ref_list)
                indexes2(find(indexes(:,3)==ref_list(i)),3)=i;
            end
            obj.Index=indexes2;
        end
        
        function [subobj index]=AddSymtypes(obj,symtypes)
            subobj=obj;
            if isa(symtypes,'sym') && isa(symtypes,'sparsesym')==0
                there_before=1:length(obj.Symtypes);

                symtypes=symtypes(:);
                if length(symtypes)>0
                    [M_unique,~,index]=unique(symtypes);
                    charM=char(M_unique);
                    if length(M_unique)>1
                        charM=charM(10:end-3);
                        charM=strsplit(charM,'], [');
                    else
                        charM={charM};
                    end
                    index=index+length(subobj.Symtypes(:));
                    
                    subobj.SymtypesChar=[subobj.SymtypesChar(:);charM(:)];
                    subobj.Symsubs=[subobj.Symsubs(:);charM(:)];
                    subobj.Symtypes=[subobj.Symtypes(:);M_unique(:)];
                else
                    index=[(1:length(subobj.Symtypes(:)))'+length(subobj.Symtypes(:))];
                end
            else
                error('The symtype to be added needs to be a symbolic array.')
            end
        end
        
        function obj=RemoveZeros(obj)
            charM=obj.SymtypesChar;
            if any(strcmp(charM,'0')) %Get rid of '0' from lists!
                j=find(strcmp(charM,'0'));
                indexes=obj.Index;
                for i=1:length(j)
                    indexes(indexes(:,3)==j(i),:)=[];
                end
                obj.Index=indexes;
                obj=RemoveUnreferenced(obj);
            end
        end
    end
    
    methods
        
        function M=full(obj)
                %Create full representation of sparse matrix
                M=sym(zeros(obj.Size));
                if length(obj.Index)>0
                    index=sub2ind(obj.Size,obj.Index(:,1),obj.Index(:,2));
                    M(index)=obj.Symtypes(obj.Index(:,3));
                end
        end
        
        function M=fullchar(obj)
                %Create full representation of sparse matrix
                M=cell(obj.Size);
                M(:)={'0'};
                index=sub2ind(obj.Size,obj.Index(:,1),obj.Index(:,2));
                M(index)=obj.SymtypesChar(obj.Index(:,3));
        end
        
        function M=fullsubchar(obj)
                %Create full representation of sparse matrix
                M=cell(obj.Size);
                M(:)={'0'};
                index=sub2ind(obj.Size,obj.Index(:,1),obj.Index(:,2));
                M(index)=obj.Symsubs(obj.Index(:,3));
        end
        
        function M=fullsub(obj)
                %Create full representation of sparse matrix
                M=sym(zeros(obj.Size));
                index=sub2ind(obj.Size,obj.Index(:,1),obj.Index(:,2));
                M(index)=sym(obj.Symsubs(obj.Index(:,3)));
        end
    
        function dim=size(obj,varargin)
            if length(varargin)==0
                dim=obj.Size;
            else
                dim=obj.Size(varargin{1});
            end
        end
        
        function num=numel(obj)
            dim=obj.Size;
            num=prod(dim(:));
        end
        
        function ind = end(obj,k,n)
            dim=obj.Size;
            if n==2
                ind=dim(k);
            elseif n==1
                ind=prod(dim(:));
            else
                error('Sparse symbolic matrices currently do not support tensors of orders beyond 2.')
            end
        end
        
        function subobj=horzcat(obj,varargin)
            Ai=varargin(:);
            Aisize=cell2mat(cellfun(@size,Ai,'UniformOutput',false));
            if all(size(obj,1)==Aisize(1:end,1))
                subobj=obj;
                for i=1:length(Ai)
                    dim=size(subobj);
                    s.type='()';
                    s.subs={[1:dim(1)] [dim(2)+1:dim(2)+Aisize(i,2)]};
                    subobj=subsasgn(subobj,s,Ai{i});
                end
            else
                error('Matrices to be horizontally concatenated need to be of same vertical size.')
            end
        end
        
        function subobj=vertcat(obj,varargin)
            Ai=varargin(:);
            Aisize=cell2mat(cellfun(@size,Ai,'UniformOutput',false));
            if all(size(obj,2)==Aisize(1:end,2))
                subobj=obj;
                for i=1:length(Ai)
                    dim=size(subobj);
                    s.type='()';
                    s.subs={[dim(1)+1:dim(1)+Aisize(i,1)] [1:dim(2)]};
                    subobj=subsasgn(subobj,s,Ai{i});
                end
            else
                error('Matrices to be horizontally concatenated need to be of same vertical size.')
            end
        end
        
        function obj=transpose(obj)
            index=obj.Index;
            obj.Index=sortrows(index(:,[2,1,3]),2)
            dim=obj.Size;
            obj.Size=dim([2,1]);
        end
        
        function obj=ctranspose(obj)
            index=obj.Index;
            obj.Index=sortrows(index(:,[2,1,3]),2)
            dim=obj.Size;
            obj.Size=dim([2,1]);
            %Complex conjugate
            symtypes=obj.Symtypes;
            symtypeschar=obj.SymtypesChar;
            symsubs=obj.Symsubs;
            symtypes=conj(symtypes);
            to_change=strcmp(symsubs,symtypeschar);
            
            symtypeschar=char(symtypes);
            if length(symtypes)>1
                symtypeschar=symtypeschar(10:end-3);
                symtypeschar=strsplit(symtypeschar,'], [');
            end
            
            symsubs(to_change)=symtypeschar(to_change);
            symsubs(to_change==0)=strcat('conj(',symsubs(to_change==0),')');
            e=strfind(symsubs,'conj(');
            who=cellfun(@length,e)==2;
            e=strfind(symsubs,'))');
            who2=cellfun(@length,e)>=1;
            if any(who & who2)
                symsubs{who & who2}=symsubs{who & who2}(11:end-2);
            end
            obj.Symtypes=symtypes;
            obj.SymtypesChar=symtypeschar;
            obj.Symsubs=symsubs;
        end
          
        function obj=uplus(obj)
            obj;
        end
        
        function obj=uminus(obj)
            symtypes=obj.Symtypes;
            symtypeschar=obj.SymtypesChar;
            symsubs=obj.Symsubs;
            symtypes=-symtypes;
            to_change=strcmp(symsubs,symtypeschar);
            
            symtypeschar=char(symtypes);
            if length(symtypes)>1
                symtypeschar=symtypeschar(10:end-3);
                symtypeschar=strsplit(symtypeschar,'], [');
            end
            
            symsubs(to_change)=symtypeschar(to_change);
            symsubs(to_change==0)=strcat('-(',symsubs(to_change==0),')');
            e=strfind(symsubs,'-(');
            who=cellfun(@length,e)==2;
            e=strfind(symsubs,'))');
            who2=cellfun(@length,e)>=1;
            if any(who & who2)
                symsubs{who & who2}=symsubs{who & who2}(5:end-2);
            end
            obj.Symtypes=symtypes;
            obj.SymtypesChar=symtypeschar;
            obj.Symsubs=symsubs;
        end
        
        function new_obj=plus(a,b)
            if isa(a,'sparsesym')==0
                a=sparsesym(a);
            end
            if isa(b,'sparsesym')==0
                b=sparsesym(b);
            end
            if a.Size==b.Size
                new_obj=sparsesym(a.Size(1),a.Size(2));
                %Overlapping elements
                indexA=a.Index;
                indexB=b.Index;
                [inter_ind,inter_a,inter_b]=intersect(indexA(:,1:2),indexB(:,1:2),'rows');
                
                index_add=[indexA(inter_a,3) indexB(inter_b,3)];
                [index_unique,ia,ib]=unique(index_add,'rows');
                symtypes=a.Symtypes(index_unique(:,1))+b.Symtypes(index_unique(:,2));
                if length(symtypes)>0
                    index=[inter_ind ib];
                    [new_obj index_overlap]=AddSymtypes(new_obj,symtypes);

                    index(:,3)=index_overlap(index(:,3));
                    new_obj.Index=index;
                end
                
                %Non-overlapping elements
                [diff_ind_a,diff_a]=setdiff(indexA(:,1:2),indexB(:,1:2),'rows');
                [diff_ind_b,diff_b]=setdiff(indexB(:,1:2),indexA(:,1:2),'rows');
                
                [index_unique_a,~,ia]=unique(indexA(diff_a,3));
                symtypes=a.Symtypes(index_unique_a(:,1));
                ia=ia+length(new_obj.Symtypes);
                if length(symtypes)>0
                    len_symtypes=length(new_obj.Symtypes);
                    index=[new_obj.Index;diff_ind_a ia];
                    len_ind=length(index(:,1));
                    [new_obj index_overlap]=AddSymtypes(new_obj,symtypes);
                    index_overlap=[1:len_symtypes,index_overlap'];
                    
                    index(len_ind+1:end,3)=index_overlap(index(len_ind+1:end,3));
                    new_obj.Index=index;
                end
                
                [index_unique_b,~,ib]=unique(indexB(diff_b,3));
                symtypes=b.Symtypes(index_unique_b(:,1));
                ib=ib+length(new_obj.Symtypes);
                if length(symtypes)>0
                    len_symtypes=length(new_obj.Symtypes);
                    index=[new_obj.Index;diff_ind_b ib];
                    len_ind=length(index(:,1));
                    [new_obj index_overlap]=AddSymtypes(new_obj,symtypes);
                    index_overlap=[1:len_symtypes,index_overlap'];
                    
                    index(len_ind+1:end,3)=index_overlap(index(len_ind+1:end,3));
                    new_obj.Index=index;
                end
                
                new_obj=RemoveZeros(new_obj);
                new_obj=UnifyDoubleSymtypes(new_obj);
            else
                error('Sparse symbolic matrices must be of same size.')
            end
        end
        
        function new_obj=minus(a,b)
            if isa(a,'sparsesym')==0
                a=sparsesym(a);
            end
            if isa(b,'sparsesym')==0
                b=sparsesym(b);
            end
            if a.Size==b.Size
                new_obj=sparsesym(a.Size(1),a.Size(2));
                %Overlapping elements
                indexA=a.Index;
                indexB=b.Index;
                [inter_ind,inter_a,inter_b]=intersect(indexA(:,1:2),indexB(:,1:2),'rows');
                
                index_add=[indexA(inter_a,3) indexB(inter_b,3)];
                [index_unique,ia,ib]=unique(index_add,'rows');
                symtypes=a.Symtypes(index_unique(:,1))-b.Symtypes(index_unique(:,2));
                if length(symtypes)>0
                    index=[inter_ind ib];
                    [new_obj index_overlap]=AddSymtypes(new_obj,symtypes);

                    index(:,3)=index_overlap(index(:,3));
                    new_obj.Index=index;
                end

                %Non-overlapping elements
                [diff_ind_a,diff_a]=setdiff(indexA(:,1:2),indexB(:,1:2),'rows');
                [diff_ind_b,diff_b]=setdiff(indexB(:,1:2),indexA(:,1:2),'rows');
                
                [index_unique_a,~,ia]=unique(indexA(diff_a,3));
                symtypes=a.Symtypes(index_unique_a(:,1));
                ia=ia+length(new_obj.Symtypes);
                if length(symtypes)>0
                    len_symtypes=length(new_obj.Symtypes);
                    len_ind=length(index(:,1));
                    index=[new_obj.Index;diff_ind_a ia];
                    [new_obj index_overlap]=AddSymtypes(new_obj,symtypes);
                    index_overlap=[1:len_symtypes,index_overlap'];
                    
                    index(len_ind+1:end,3)=index_overlap(index(len_ind+1:end,3));
                    new_obj.Index=index;
                end
                
                [index_unique_b,~,ib]=unique(indexB(diff_b,3));
                symtypes=-b.Symtypes(index_unique_b(:,1));
                ib=ib+length(new_obj.Symtypes);
                if length(symtypes)>0
                    len_symtypes=length(new_obj.Symtypes);
                    len_ind=length(index(:,1));
                    index=[new_obj.Index;diff_ind_b ib];
                    [new_obj index_overlap]=AddSymtypes(new_obj,symtypes);
                    index_overlap=[1:len_symtypes,index_overlap'];
                    
                    index(len_ind+1:end,3)=index_overlap(index(len_ind+1:end,3));
                    new_obj.Index=index;
                end
                new_obj=RemoveZeros(new_obj);
                new_obj=UnifyDoubleSymtypes(new_obj);
            else
                error('Sparse symbolic matrices must be of same size.')
            end
        end
        
        function new_obj=times(a,b)
            if isa(a,'sparsesym')==0
                a=sparsesym(a);
            end
            if isa(b,'sparsesym')==0
                b=sparsesym(b);
            end
            if a.Size==b.Size
                new_obj=sparsesym(a.Size(1),a.Size(2));
                %Overlapping elements
                indexA=a.Index;
                indexB=b.Index;
                [inter_ind,inter_a,inter_b]=intersect(indexA(:,1:2),indexB(:,1:2),'rows');
                
                index_add=[indexA(inter_a,3) indexB(inter_b,3)];
                [index_unique,ia,ib]=unique(index_add,'rows');
                symtypes=a.Symtypes(index_unique(:,1)).*b.Symtypes(index_unique(:,2));
                if length(symtypes)>0
                    index=[inter_ind ib];
                    [new_obj index_overlap]=AddSymtypes(new_obj,symtypes);

                    index(:,3)=index_overlap(index(:,3));
                    new_obj.Index=index;
                end

                new_obj=RemoveZeros(new_obj);
                new_obj=UnifyDoubleSymtypes(new_obj);
            else
                error('Sparse symbolic matrices must be of same size.')
            end
        end
        
        function new_obj=mtimes(a,b)
            if isa(a,'sparsesym')==0
                a=sparsesym(a);
            end
            if isa(b,'sparsesym')==0
                b=sparsesym(b);
            end
            if a.Size(2)==b.Size(1)
                new_obj=sparsesym(a.Size(1),b.Size(2));
                %Overlapping elements
                indexA=sortrows(a.Index,[2,1]);
                indexB=sortrows(b.Index,[1,2]);
                uniqueA=unique(indexA(:,2));
                uniqueB=unique(indexB(:,1));
                column_ind=unique(indexB(:,2));
                row_ind=unique(indexA(:,1));
                indexA=indexA(ismember(indexA(:,2),column_ind),:);
                indexB=indexB(ismember(indexB(:,1),row_ind),:);
                symsA=a.Symtypes;
                symsB=b.Symtypes;
                
                index=[];
                symtypes=[];
                
                for i=1:length(column_ind)
                    for j=1:length(row_ind)
                        As=indexA(find(indexA(:,1)==row_ind(j)),:)
                        Bs=indexB(find(indexB(:,2)==column_ind(i)),:)
                        [inter,ia,ib]=intersect(As(:,2),Bs(:,1))
                        new_sym=symsA(As(ia,3)).'*symsB(Bs(ib,3));
                        symtypes=[symtypes,new_sym]
                        index=[index;row_ind(j),column_ind(i),length(symtypes)];
                    end
                end
                if length(symtypes)>0
                    [new_obj index_overlap]=AddSymtypes(new_obj,symtypes);

                    index(:,3)=index_overlap(index(:,3));
                    new_obj.Index=index;
                    
                    new_obj=RemoveZeros(new_obj);
                    new_obj=UnifyDoubleSymtypes(new_obj);
                end
            else
                error('Sparse symbolic matrices must be of same size.')
            end
        end
        
        function new_obj=power(a,b)
            if isa(a,'sparsesym')==0
                a=sparsesym(a);
            end
            if isa(b,'sparsesym')==0
                b=sparsesym(b);
            end
            if a.Size==b.Size
                new_obj=sparsesym(a.Size(1),a.Size(2));
                %Overlapping elements
                indexA=a.Index;
                indexB=b.Index;
                [inter_ind,inter_a,inter_b]=intersect(indexA(:,1:2),indexB(:,1:2),'rows');
                
                index_add=[indexA(inter_a,3) indexB(inter_b,3)];
                [index_unique,ia,ib]=unique(index_add,'rows');
                symtypes=a.Symtypes(index_unique(:,1))^b.Symtypes(index_unique(:,2));
                if length(symtypes)>0
                    index=[inter_ind ib];
                    [new_obj index_overlap]=AddSymtypes(new_obj,symtypes);

                    index(:,3)=index_overlap(index(:,3));
                    new_obj.Index=index;
                end

                new_obj=RemoveZeros(new_obj);
                new_obj=UnifyDoubleSymtypes(new_obj);
            else
                error('Sparse symbolic matrices must be of same size.')
            end
         end
        
         function new_obj=ldivide(a,b)
            if isa(a,'sparsesym')==0
                a=sparsesym(a);
            end
            if isa(b,'sparsesym')==0
                b=sparsesym(b);
            end
            if a.Size==b.Size
                new_obj=sparsesym(a.Size(1),a.Size(2));
                %Overlapping elements
                indexA=a.Index;
                indexB=b.Index;
                [inter_ind,inter_a,inter_b]=intersect(indexA(:,1:2),indexB(:,1:2),'rows');
                
                index_add=[indexA(inter_a,3) indexB(inter_b,3)];
                [index_unique,ia,ib]=unique(index_add,'rows');
                symtypes=a.Symtypes(index_unique(:,1)).\b.Symtypes(index_unique(:,2));
                if length(symtypes)>0
                    index=[inter_ind ib];
                    [new_obj index_overlap]=AddSymtypes(new_obj,symtypes);

                    index(:,3)=index_overlap(index(:,3));
                    new_obj.Index=index;
                end

                new_obj=RemoveZeros(new_obj);
                new_obj=UnifyDoubleSymtypes(new_obj);
            else
                error('Sparse symbolic matrices must be of same size.')
            end
         end
        
         function new_obj=rdivide(a,b)
            if isa(a,'sparsesym')==0
                a=sparsesym(a);
            end
            if isa(b,'sparsesym')==0
                b=sparsesym(b);
            end
            if a.Size==b.Size
                new_obj=sparsesym(a.Size(1),a.Size(2));
                %Overlapping elements
                indexA=a.Index;
                indexB=b.Index;
                [inter_ind,inter_a,inter_b]=intersect(indexA(:,1:2),indexB(:,1:2),'rows');
                
                index_add=[indexA(inter_a,3) indexB(inter_b,3)];
                [index_unique,ia,ib]=unique(index_add,'rows');
                symtypes=a.Symtypes(index_unique(:,1))./b.Symtypes(index_unique(:,2));
                if length(symtypes)>0
                    index=[inter_ind ib];
                    [new_obj index_overlap]=AddSymtypes(new_obj,symtypes);

                    index(:,3)=index_overlap(index(:,3));
                    new_obj.Index=index;
                end

                new_obj=RemoveZeros(new_obj);
                new_obj=UnifyDoubleSymtypes(new_obj);
            else
                error('Sparse symbolic matrices must be of same size.')
            end
         end
        
        function obj=Substitute(obj,to_substitute,substitution)
            if isa(to_substitute,'numeric')
                %Remove indexes
                if isa(substitution,'cell')
                    obj.Symsubs(to_substitute)=substitution;
                elseif isa(substitution,'char')
                    obj.Symsubs(to_substitute)={substitution};
                end
            elseif isa(to_substitute,'char')
                to_del_ind=strcmp(obj.SymtypesChar,to_substitute);
                if isa(substitution,'cell')
                    obj.Symsubs(to_del_ind)=substitution;
                elseif isa(substitution,'char')
                    obj.Symsubs(to_del_ind)={substitution};
                end
            elseif isa(to_substitute,'sym') && isa(to_substitute,'sparsesym')==0
                if length(to_substitute)>1
                    to_substitute=to_substitute(:);
                    to_substitute=char(to_substitute);
                    to_substitute=to_substitute(10:end-3);
                    to_substitute=strsplit(to_substitute,'], [');
                else
                    to_substitute=char(to_substitute);
                end
                [~,i_a,i_b]=intersect(to_substitute,obj.SymtypesChar);
                if isa(substitution,'cell')
                    obj.Symsubs(i_b)=substitution;
                elseif isa(substitution,'char')
                    obj.Symsubs(i_b)={substitution};
                end
            elseif isa(to_substitute,'cell')
                if isa(to_substitute{1},'sym')
                    to_substitute=cellfun(@char,to_substitute,'UniformOutput',false);
                end
                if isa(to_substitute{1},'char')
                    [~,i_a,i_b]=intersect(to_substitute,obj.SymtypesChar);
                    obj.Symsubs(i_b)=substitution;
                else
                    error('Symtypes that shall be removed need to be specified via indexes, symbolic or character input or via a cell of character or symbolic inputs.')
                end
            else
                error('Symtypes that shall be removed need to be specified via indexes, symbolic or character input or via a cell of character or symbolic inputs.')
            end
        end
    end
    
end

