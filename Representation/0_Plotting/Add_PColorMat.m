function [ bra_str ket_str indexes ] = Add_PColorMat( ax,H,mode )
%PCOLORMAT creates pcolor plots of matrices (of any kind) adding support
%for fock and index labels and adding elements to the matrix to represent
%the full matrix
%mode = 'fock_bin' binary ket's and bra's
%       {'fock_bin_real_abs_imag',n,N} Number ket's and bra's -> add png or tex for
%       no_edge -> no edge in surf
%       phase   -> color determined by phase , brightness by absolute value
%           phase overrides real, abs (default), imag

if isa(mode,'char')
    mode={mode};
end
if length(mode)<3 && length(mode)>1
    error('Mode must also specify the maximum number of particles per mode (n) and the number of modes N.')
elseif length(mode)==1
    mode{2}=1;
    mode{3}=log2(size(H,1));
end


if length(mode)<5
    max_ticks=16;
else
    max_ticks=mode{5};
end
tick_num=min([length(H),max_ticks]);

if any(strfind(mode{1},'fock'))>0
    if tick_num~=length(H) 
        div=length(H)/tick_num;  
        indexes=1:div:length(H); 

        if any(strfind(mode{1},'tex'))==0
            bra_str=bra_fock_str(indexes,mode);
            ket_str=ket_fock_str(indexes,mode);
        else
            bra=bra_fock_tex(1:length(H),mode);
            ket=ket_fock_tex(1:length(H),mode);        
        end
    else
        indexes=1:length(H);

        if any(strfind(mode{1},'tex'))==0
            bra_str=bra_fock_str(1:length(H),mode);
            ket_str=ket_fock_str(1:length(H),mode);
        else
            bra=bra_fock_tex(1:length(H),mode);
            ket=ket_fock_tex(1:length(H),mode);
        end
    end
    if any(strfind(mode{1},'tex'))>0
        for i=1:length(bra)
            bra_str{i}=['$',bra{i}(2:end),'$'];
            ket_str{i}=['$',ket{i}(2:end),'$'];
        end
    end
elseif any(strfind(mode{1},'none'))
    if tick_num~=length(H)
        div=length(H)/tick_num;  
        indexes=1:div:length(H);
    else
        indexes=1:length(H);
    end
    bra_str=cell(length(indexes),1);
    ket_str=cell(length(indexes),1);
    for i=1:length(indexes)
        ket_str{i}='';
        bra_str{i}='';
    end
else %if length(strfind(mode{1},'index'))==1
    if tick_num~=length(H)
        div=length(H)/tick_num;  
        indexes=1:div:length(H);
    else
        indexes=1:length(H);
    end
    
    max_length=1;
    bra_str=cell(length(indexes),1);
    ket_str=cell(length(indexes),1);
    long_str=num2str(max(indexes));
    type=['%0' num2str(length(long_str)) 'd'];
    for i=1:length(indexes)
        ket_str{i}=num2str(indexes(i),type);
    end
    ket_str=strcat({'['},ket_str,{']'});
    bra_str=ket_str;
end

    if any(strfind(mode{1},'phase'))
        U=H;
        cmap=circshift(hsv(256),0);%floor(256/2));  %phasecolormap
        bg=[0.5 0.5 1]*0.20;
        mapper=linspace(0,2*pi*255/256,256);
        an=angle(U);
        an(an<0)=an(an<0)+2*pi;
        m=max(abs(U(:)));
        for i=1:size(an,1)*size(an,2)
            [alpha,beta]=min(abs(an(i)-mapper));
            an(i)=beta;
        end
        in=abs(H)/m;
        C1(:)=cmap(an(:),1).*in(:)+bg(1)*(1-in(:));
        C2(:)=cmap(an(:),2).*in(:)+bg(2)*(1-in(:));
        C3(:)=cmap(an(:),3).*in(:)+bg(3)*(1-in(:));
        C1=reshape(C1,size(an,1),size(an,2));
        C2=reshape(C2,size(an,1),size(an,2));
        C3=reshape(C3,size(an,1),size(an,2));
        C(:,:,1)=C1;
        C(:,:,2)=C2;
        C(:,:,3)=C3;
        g=imagesc(ax,C);
    elseif any(strfind(mode{1},'real'))
        U=[H zeros(length(H),1); zeros(1,length(H)+1)];
        g=pcolor(ax,real(U));
    elseif any(strfind(mode{1},'imag'))
        U=[H zeros(length(H),1); zeros(1,length(H)+1)];
        g=pcolor(ax,imag(U));
    else
        U=[H zeros(length(H),1); zeros(1,length(H)+1)];
        g=pcolor(ax,abs(U));
    end
    mode2='';
    if max(size(H))>50
        if any(strfind(mode{1},'phase'))==0
            set(g, 'EdgeColor', 'none')
        end
    elseif length(strfind(mode{1},'no_edge'))>0
        mode2='noedge';
    end
    %Add axis ticks

        if tick_num==length(U) %Take all ticks
            ax.XTick=indexes+0.5;
            ax.YTick=indexes+0.5;
            ax.XTickLabel=ket_str(:);
            ax.YTickLabel=bra_str(:);
        else 
            ax.XTick=indexes+0.5;
            ax.YTick=indexes+0.5;
            ax.XTickLabel=ket_str(:);
            ax.YTickLabel=bra_str(:);
        end
        
    set(ax,'XTickLabelRotation',45)
    axis(ax,[1 length(U)+1 1 length(U)+1  0 1])
    set(ax,'XAxisLocation','top','YAxisLocation','left','ydir','reverse');
    axis(ax,'equal');
    view(ax,0,90);
    %drawnow();
end

