function [ filename ] = Unitary_Anim( H,t,mode,sizle,file_name,sizle2 )
%UNITARY_ANIM animates the unitary evolution of a Hamiltonian over time t
%mode = 'fock_bin' binary ket's and bra's
%       {'fock_bin',n,N} Number ket's and bra's -> add gif or mp4 for
%       movies
%       no_edge -> no edge in surf
if nargin==2
    mode='bin';
    sizle=600;
    file_name='';
    sizle2=1260;
elseif nargin==3
    sizle=600;
    file_name='';
    sizle2=1260;
elseif nargin==4
    file_name='';
    sizle2=1260;
elseif nargin==5
    sizle2=1260;
end

if isa(mode,'char')
    mode={mode};
end
if length(mode)<3 && length(mode)>1
    error('Mode must also specify the maximum number of particles per mode (n) and the number of modes N.')
elseif length(mode)==1
    mode{2}=1;
    mode{3}=log2(size(H,1));
end
%if length(strfind(mode{1},'fock'))==0
%    mode{1}=['fock_' mode{1}];
%end

%Create figure
f=figure('Position',[100 100 round(16/9*sizle) sizle]);
h=axes('Parent',f);

%Create more t steps if necessary
if length(t)==2
    t=linspace(min(t),max(t),100);
elseif length(t)==3 && t(3)==round(t(3))
    t=linspace(min(t([1:2])),max(t([1:2])),t(3));
end
  
%% Save Files
if any(strfind(mode{1},'gif')) | any(strfind(mode{1},'mp4'))
    curr_folder=pwd;
    save_folder=fullfile(curr_folder,'Animations');
    if exist(save_folder)==0
        mkdir('Animations');
    end

    if length(file_name)>0
        save_folder2=fullfile(save_folder,file_name);
        if exist(save_folder2)==0
            mkdir(save_folder2);
        end

        if any(strfind(mode{1},'gif'))==1
            filename=[fullfile(save_folder2,file_name),'.gif'];
        elseif any(strfind(mode{1},'mp4'))==1
            filename=[fullfile(save_folder2,file_name)];
        end
    else
        file_name='new_name';
        save_folder2=fullfile(save_folder,file_name);
        if exist(save_folder2)==0
            mkdir(save_folder2);
        end

        if any(strfind(mode{1},'gif'))==1
            filename=[fullfile(save_folder2,file_name),'.gif'];
        elseif any(strfind(mode{1},'mp4'))==1
            filename=[fullfile(save_folder2,file_name)];
        end
    end
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

        if any(strfind(mode{1},'mp4'))==0
            bra_str=bra_fock_str(indexes,mode);
            ket_str=ket_fock_str(indexes,mode);
        else
            bra=bra_fock_tex(1:length(H),mode);
            ket=ket_fock_tex(1:length(H),mode);        
        end
    else
        indexes=1:length(H);

        if any(strfind(mode{1},'mp4'))==0
            bra_str=bra_fock_str(1:length(H),mode);
            ket_str=ket_fock_str(1:length(H),mode);
        else
            bra=bra_fock_tex(1:length(H),mode);
            ket=ket_fock_tex(1:length(H),mode);
        end
    end
    if any(strfind(mode{1},'mp4'))>0
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

for i=1:length(t)
    U=[H2U(H,t(i)) zeros(length(H),1); zeros(1,length(H)+1)];
    if any(strfind(mode{1},'real'))
        g=pcolor(h,real(U));
    elseif any(strfind(mode{1},'imag'))
        g=pcolor(h,imag(U));
    else
        g=pcolor(h,abs(U));
    end
    mode2='';
    if max(size(H))>50
        set(g, 'EdgeColor', 'none')
    elseif length(strfind(mode{1},'no_edge'))>0
        mode2='noedge';
    end
    title(h,['t=' num2str(t(i))]);
    
    %Add axis ticks
        if tick_num==length(U) %Take all ticks
            h.XTick=indexes+0.5;
            h.YTick=indexes+0.5;
            h.XTickLabel=ket_str(:);
            h.YTickLabel=bra_str(:);
        else 
            h.XTick=indexes+0.5;
            h.YTick=indexes+0.5;
            h.XTickLabel=ket_str(:);
            h.YTickLabel=bra_str(:);
        end
        
    set(h,'XTickLabelRotation',45)
    axis(h,[1 length(U)+1 1 length(U)+1  0 1])
    set(h,'XAxisLocation','top','YAxisLocation','left','ydir','reverse');
    axis equal;
    view(0,90);
    drawnow;
    
    if any(strfind(mode{1},'gif'))==1
        fprintf([num2str(i),'/',num2str(length(t)),'\n']);
        if length(file_name)>0
            frame=getframe(h); 
            im=frame2im(frame); 
            [imind,cm]=rgb2ind(im,256); 
            % Write to the GIF File 
            if i==1 
                imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',0.05); 
            else 
                imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.05);
            end 
        else
            pause(0.2);
        end
    elseif any(strfind(mode{1},'mp4'))==1
        fprintf([num2str(i),'/',num2str(length(t)),'\n']);
        file=PColor2Tikz(fullfile(save_folder2,file_name),mode2);
        cd(save_folder2)
        A=lualatex2pdf2preview( file,0,600 );
        A=imresize(A,[sizle2,NaN]);
        cd(curr_folder);
        if i==1
            imwrite(A,fullfile(save_folder2,'frame_title.png'));
            v=VideoWriter(fullfile(save_folder2,[file_name '.avi']),'Uncompressed AVI');
            open(v)
        end
        writeVideo(v,A);
    end
end
if any(strfind(mode{1},'mp4'))==1
    close(v);
end
end

