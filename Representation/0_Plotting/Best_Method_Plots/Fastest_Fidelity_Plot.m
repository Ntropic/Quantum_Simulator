function [ h ] = Fastest_Fidelity_Plot( h,mode,cost_maps,N_m,n_m,fid,data )
%FASTEST_FIDELITY_PLOT Plots the fastest Fidelity comparing different
%approaches
if length(h)==2
    ax=h(2);
    h=h(1);
end
names1=data{6};
names2=data{7};
names={};
for i=1:length(names1)
    for j=1:length(names2)
        names={names{:},[names1{i} ' (' names2{j} ')']};
    end
end

if length(h)==0
    h=figure('Name','Fastest Fidelity Plot','Position',[100 100 1000 600]);
    if length(strfind(mode,'redraw'))>0
        ax=subplot('Position',[0.075 0.05 0.8 0.9]);
        fidelity_slider=uicontrol('Parent',h,'Style','slider','Tag','fid_slide','Position',[925 150 20 310],'Min',0,'Max',10,'Value',-log10(1-fid),'SliderStep',[1/100 1/100]);
        addlistener(fidelity_slider,'ContinuousValueChange',@(src,evt)Redraw(src,evt,ax,h,data,mode));
        error_text=uicontrol('Parent',h,'Style','text','Tag','error_text','Position',[885 90 100 50],'String',['1-F=' num2str(1-fid,'%1.1e')],'BackgroundColor',get(h,'Color'),'HorizontalAlignment','center');
    else
        ax=subplot(1,1,1);    
    end
end

len=length(cost_maps);
sizes=size(cost_maps{1});

cmap=jet(len); %Color coding for the different approaches

if length(strfind(mode,'pixel'))>0
    %First plot some surfaces in the background to get legend entries
    for i=1:len
        cost=cost_maps{i};
        surf(ax,n_m(1:2),N_m(1:2),-1*[1 1;1 1],i*ones(2))
        colormap(cmap);
        hold on;
    end
    
    %Find minimum costs
    sizes=size(cost_maps{1});
    len=prod(sizes);
    for i=1:length(cost_maps)
        %Create one big matrix for comparison
        co_map(:,i)=cost_maps{i}(:);
    end
    [a minind]=min(co_map,[],2);
    col1=reshape(cmap(minind,1),sizes);
    col2=reshape(cmap(minind,2),sizes);
    col3=reshape(cmap(minind,3),sizes);
    col(:,:,1)=col1;
    col(:,:,2)=col2;
    col(:,:,3)=col3;
    image(ax,n_m,N_m,col);
    legend(names,'Location', 'NorthEast');
    hold off;
else %Surface plot
    min_cost=min(min(cost_maps{1}));
    for i=2:len
        if min(min(cost_maps{i}))>min_cost
            min_cost=min(min(cost_maps{i}));
        end
    end
    for i=1:len
        cost=cost_maps{i};
        surf(ax,n_m,N_m,1./cost*min_cost,i*ones(size(cost)))
        colormap(cmap);
        hold on;
    end
    legend(names,'Location', 'NorthEast');
end
xlabel('# Photons')
ylabel('# WGs')
set(ax,'ytick',N_m);
set(ax,'xtick',n_m);
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
    slide_num=get(fidelity_slider,'Value');                       
    %{'Mesh','Surface','PColor','Contour'}
    
    err=10^(-slide_num);
    fid=1-err;
    %Change ui text elements
    set(error_text,'String',['1-F=' num2str(err,'%1.1e')]);
    
    %Renew Plot
    s_q_pred=data{1};
    F_pred=data{2};
    F=data{3};
    s_q=data{4};
    costs=data{5};
    %names1=data{6};
    %names2=data{7};
    N_m=data{8};
    n_m=data{9};
    cost_maps=Fastest_Fidelity(fid,F,s_q,costs,F_pred,s_q_pred);
    h=Fastest_Fidelity_Plot([h ax],mode,cost_maps,N_m,n_m,fid,data);
end