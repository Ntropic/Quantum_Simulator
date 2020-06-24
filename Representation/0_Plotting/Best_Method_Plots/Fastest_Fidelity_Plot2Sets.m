function [ h ] = Fastest_Fidelity_Plot2Sets( h,mode,fid,data )
%FASTEST_FIDELITY_PLOT2SETS Plots the fastest Fidelity comparing different
%approaches for 2 datasets
%data={s_q_pred,s_q_pred2,F_pred,F_pred2,F,F2,s_q,s_q2,costs,names1,names2,N_m,N_m2,n_m,n_m2};
if length(h)==2
    ax=h(2);
    h=h(1);
end
%Extract data
s_q_pred=data{1};
s_q_pred2=data{2};
F_pred=data{3};
F_pred2=data{4};
F=data{5};
F2=data{6};
s_q=data{7};
s_q2=data{8};
costs=data{9};
names1=data{10};
names2=data{11};
N_m=data{12};
N_m2=data{13};
n_m=data{14};
n_m2=data{15};
cost_maps=data{16};
cost_maps2=data{17};
extrap_maps=data{18};
extrap_maps2=data{19};

%Create names
names={};
for i=1:length(names1)
    for j=1:length(names2)
        names={names{:},[names1{i} ' (' names2{j} ')']};
    end
end

if length(h)==0
    h=figure('Name','Fastest Fidelity Plot','Position',[100 100 1000 800]);
    if length(strfind(mode,'redraw'))>0
        ax=subplot('Position',[0.075 0.05 0.8 0.9]);
        fidelity_slider=uicontrol('Parent',h,'Style','slider','Tag','fid_slide','Position',[925 150 20 310],'Min',0,'Max',15,'Value',-log10(1-fid),'SliderStep',[1/100 1/100]);
        addlistener(fidelity_slider,'ContinuousValueChange',@(src,evt)Redraw(src,evt,ax,h,data,mode));
        error_text=uicontrol('Parent',h,'Style','text','Tag','error_text','Position',[885 90 100 50],'String',['1-F=' num2str(1-fid,'%1.1e')],'BackgroundColor',get(h,'Color'),'HorizontalAlignment','center');
    else
        ax=subplot(1,1,1);    
    end
end

skips=0; %How many colors are skipped between blocks
len=length(cost_maps);
sizes=size(cost_maps{1});
len2=length(cost_maps2);
sizes2=size(cost_maps2{1});

len_names_1=length(names1);
len_names_2=length(names2);
len_col=len+skips*(len_names_1-1);
cmap=jet(len_col); %Color coding for the different approaches

    %First plot some surfaces in the background to get legend entries
    for i=1:len
        cost=cost_maps{i};
        surf(ax,n_m(1:2),N_m(1:2),-1*[1 1;1 1],(i+skips*floor((i-1)/len_names_2))*ones(2))
        colormap(cmap);
        hold on;
    end
    
    %Find minimum costs
    for i=1:length(cost_maps)
        %Create one big matrix for comparison
        co_map(:,i)=cost_maps{i}(:);
    end
    [a minind]=min(co_map,[],2);
    col1=reshape(cmap((minind+skips*floor((minind-1)/len_names_2)),1),sizes);
    col2=reshape(cmap((minind+skips*floor((minind-1)/len_names_2)),2),sizes);
    col3=reshape(cmap((minind+skips*floor((minind-1)/len_names_2)),3),sizes);
    col(:,:,1)=col1;
    col(:,:,2)=col2;
    col(:,:,3)=col3;
    image(ax,n_m,N_m,col);

    %Find 2nd set minimum costs
    
    for i=1:length(cost_maps)
        %Create one big matrix for comparison
        co_map2(:,i)=cost_maps2{i}(:);
    end
    [a minind2]=min(co_map2,[],2);
    col12=reshape(cmap(minind2,1),sizes2);
    col22=reshape(cmap(minind2,2),sizes2);
    col32=reshape(cmap(minind2,3),sizes2);
    col_2(:,:,1)=col12;
    col_2(:,:,2)=col22;
    col_2(:,:,3)=col32;
    image(ax,n_m2,N_m2,col_2);
    
    if length(strfind(mode,'extrap'))>0
        %Hatchet Plot of extrapolated values
        hatchet_map=zeros(sizes);
        hatchet_map2=zeros(sizes2);
        for k=1:length(minind)
            hatchet_map(k)=extrap_maps{minind(k)}(k);
        end
        for k=1:length(minind2)
            hatchet_map2(k)=extrap_maps2{minind2(k)}(k);
        end
        u_n_m=unique([n_m n_m2]);
        u_N_m=unique([N_m N_m2]);
        hatch=zeros(max(u_N_m)-min(u_N_m),max(u_n_m)-min(u_n_m));
        hatch((min(N_m)-min(u_N_m)+1):(max(N_m)-min(u_N_m)+1),(min(n_m)-min(u_n_m)+1):(max(n_m)-min(u_n_m)+1))=hatchet_map;
        hatch((min(N_m2)-min(u_N_m)+1):(max(N_m2)-min(u_N_m)+1),(min(n_m2)-min(u_n_m)+1):(max(n_m2)-min(u_n_m)+1))=hatchet_map2;
        if any(hatch(:)>0)
            %contour(ax,u_n_m,u_N_m,hatch);
            hatch_mat(ax,'',hatch,5,u_n_m,u_N_m )
        end
    end
    legend(names,'Location', 'NorthEast');
    hold off;

xlabel('# Photons')
ylabel('# WGs')
set(ax,'ytick',unique([N_m N_m2]));
set(ax,'xtick',unique([n_m n_m2]));
axis equal;
view(0,90);

axis tight;
hold off;
drawnow();
end

function Redraw(a,b,ax,h,data,mode)
    %REDRAW creates the plots according to the settings in the interface
    %Find uicontrol elements
    fidelity_slider=findobj(h,'Tag','fid_slide');
    error_text=findobj(h,'Tag','error_text');
   
    
    %Extract values from uicontrol elements
    slide_num=get(fidelity_slider,'Value');   %Boundaries of the Plot                     %{'Mesh','Surface','PColor','Contour'}
    
    %Calculate new fidelity condition
    err=10^(-slide_num);
    fid=1-err;
    %Change ui text elements
    set(error_text,'String',['1-F=' num2str(err,'%1.1e')]);
    
    %Extract data
    s_q_pred=data{1};
    s_q_pred2=data{2};
    F_pred=data{3};
    F_pred2=data{4};
    F=data{5};
    F2=data{6};
    s_q=data{7};
    s_q2=data{8};
    costs=data{9};
    %names1=data{10};
    %names2=data{11};
    %N_m=data{12};
    %N_m2=data{13};
    %n_m=data{14};
    %n_m2=data{15};
    %cost_maps=data{16};
    %cost_maps2=data{17};
    %extrap_maps=data{18};
    %extrap_maps2=data{19};
    
    %Calculate new cost maps
    [cost_maps extrap_maps]=Fastest_Fidelity(fid,F,s_q,costs,F_pred,s_q_pred);
    [cost_maps2 extrap_maps2]=Fastest_Fidelity(fid,F2,s_q2,costs,F_pred2,s_q_pred2);
    %Renew Plot
    data{16}=cost_maps;
    data{17}=cost_maps2;
    data{18}=extrap_maps;
    data{19}=extrap_maps2;
    h=Fastest_Fidelity_Plot2Sets([h ax],mode,fid,data);
end