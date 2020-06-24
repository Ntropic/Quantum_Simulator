function h=Compare_Fidelity_Distribution_Circ_XY_All_Interaction()
    %Compare_Fidelity_Distribution_Circ_XY_All_Interaction.m
    clear all;
    close all;
    clc;

    p=pwd;
    if any(strfind(p,'\'));
        elem=strsplit(p,'\');
    else
        elem=strsplit(p,'/');
    end
    shortened=fullfile(elem{1:end-3});
    addpath(genpath(shortened));
    %Load gates
    load('../../../Gates_Table/elem_gates.mat','-mat')
    load('../../../Gates_Table/comp_gates.mat','-mat')
    load('../XY_quantum_gates.mat','-mat')

    %Trotter step number
    steps=20;

    %Iteration time
    t=pi/4;

    n_m=3:10;


    bins=200;
    how_many=10^4;

    load(['comp_fid_dist_N_' num2str(max(n_m)) '.mat']);
    steps=data{1};
    x=data{2};
    xo=data{3};
    Order_f1_av=data{4};
    Order_f1_min=data{5};
    Order_f1_noon=data{6};
    order_hist=data{7};
    f1_bins=data{8};
    f1_av=data{9};
    f1_min=data{10};
    f1_noon=data{11};
    f1_bins=data{12};
    n_m=data{13};

    %% Plot the figures
    i=1;
    %% 1-Fidelity Plot
    h=figure('Position',[100 100 1200 800]);
    set(h,'color',0.98*[1 1 1]);

    ax1=subplot(1,2,1);
    lin_x=linspace(min(x{i}),max(x{i}),length(x{i}));
    dx=lin_x(2)-lin_x(1);
    lin_x=lin_x+dx/2-(x{i}(end-1)-x{i}(end));
    imagesc(ax1,1:steps,lin_x,f1_bins{i}(end:-1:1,:),'HandleVisibility','off');
    je=gray(512);
    je=je(255:512,:);
    je=je(end:-1:1,:);
    colormap(je)
    caxis([0 max(max(f1_bins{i}(1:end-1,:)))])
    set(ax1,'YDir','normal')
    set(ax1, 'YScale', 'log')
    hold on;
    plot3(ax1,1:steps,f1_av{i},ones(size(f1_av{i})),'k')
    plot3(ax1,1:steps,f1_min{i},ones(size(f1_min{i})),'--k')
    plot3(ax1,1:steps,f1_noon{i},ones(size(f1_noon{i})),'--r')
    xlabel(ax1,'Trotter steps')
    ylabel(ax1,'1-Fidelity')
    legend(ax1,{'Average Fidelity','Minimum Fidelity','NOON Fidelity'},'Location','northeast')
    axis(ax1,[1 steps min(x{i}) max(x{i}) 0 1]);
    view(ax1,0,90)
    hold off;
    drawnow;

    %% Order of Convergence of 1-Fidelity Plot
    ax2=subplot(1,2,2);
    lin_xo=linspace(min(xo{i}),max(xo{i}),length(xo{i}));
    dxo=lin_xo(2)-lin_xo(1);
    lin_xo=lin_xo-dxo/2;
    imagesc(ax2,1:(steps-1),lin_xo(end:-1:1),order_hist{i}(end:-1:1,:),'HandleVisibility','off');

    colormap(je)
    caxis([0 max(max(order_hist{i}))])
    set(ax2,'YDir','normal')
    hold on;
    plot3(ax2,1:(steps-1),Order_f1_av{i},ones(size(Order_f1_av{i})),'k')
    plot3(ax2,1:(steps-1),Order_f1_min{i},ones(size(Order_f1_min{i})),'--k')
    plot3(ax2,1:(steps-1),Order_f1_noon{i},ones(size(Order_f1_noon{i})),'--r')
    xlabel(ax2,'Trotter steps')
    ylabel(ax2,'Order of Convergence')
    legend(ax2,{'Average Fidelity','Minimum Fidelity','NOON Fidelity'},'Location','northeast')
    axis(ax2,[1 steps-1 min(xo{i}) max(xo{i}) 0 1]);
    view(ax2,0,90)
    hold off;
    drawnow;
    
    %% UI Elements

    photon_slider=uicontrol('Parent',h,'Style','slider','Tag','ph_slide','Position',[1125 150 20 310],'Min',1,'Max',length(n_m),'Value',1,'SliderStep',[1/length(n_m) 1/length(n_m)]);
    addlistener(photon_slider,'ContinuousValueChange',@(src,evt)Redraw(src,evt,[ax1 ax2],h,data));
    photon_text=uicontrol('Parent',h,'Style','text','Tag','photon_text','Position',[1085 90 100 50],'String',['N=' num2str(n_m(1)-1)],'BackgroundColor',get(h,'Color'),'HorizontalAlignment','center');
    
    
end
function Redraw(a,b,ax,h,data)
    ax1=ax(1);
    ax2=ax(2);
    
    %data={steps,x,xo,Order_f1_av,Order_f1_min,Order_f1_noon,order_hist,f1_bins,f1_av,f1_min,f1_noon,f1_bins};
    steps=data{1};
    x=data{2};
    xo=data{3};
    Order_f1_av=data{4};
    Order_f1_min=data{5};
    Order_f1_noon=data{6};
    order_hist=data{7};
    f1_bins=data{8};
    f1_av=data{9};
    f1_min=data{10};
    f1_noon=data{11};
    f1_bins=data{12};
    n_m=data{13};

    %Find uicontrol elements
    photon_slider=findobj(h,'Tag','ph_slide');
    photon_text=findobj(h,'Tag','photon_text');

    slide_num=round(get(photon_slider,'Value'));  
    set(photon_text,'String',['N=' num2str(n_m(slide_num)-1)]);
    
    N=n_m(slide_num)-1;
    i=slide_num;
    %% 1-Fidelity Plot
    lin_x=linspace(min(x{i}),max(x{i}),length(x{i}));
    dx=lin_x(2)-lin_x(1);
    lin_x=lin_x+dx/2-(x{i}(end-1)-x{i}(end));
    hold(ax1,'off');

    imagesc(ax1,1:steps,lin_x,f1_bins{i}(end:-1:1,:),'HandleVisibility','off');
    caxis([0 max(max(f1_bins{i}(10:end-10,2:end-2)))])
    je=gray(512);
    je=je(255:512,:);
    je=je(end:-1:1,:);
    colormap(je)

    set(ax1,'YDir','normal')
    set(ax1, 'YScale', 'log')
    hold(ax1,'on');
    
    plot3(ax1,1:steps,f1_av{i},ones(size(f1_av{i})),'k')
    plot3(ax1,1:steps,f1_min{i},ones(size(f1_min{i})),'--k')
    plot3(ax1,1:steps,f1_noon{i},ones(size(f1_noon{i})),'--r')
    xlabel(ax1,'Trotter steps')
    ylabel(ax1,'1-Fidelity')
    legend(ax1,{'Average Fidelity','Minimum Fidelity','NOON Fidelity'},'Location','northeast')
    axis(ax1,[1 steps min(x{i}) max(x{i}) 0 1]);
    view(ax1,0,90)
    hold(ax1,'off');

    %% Order of Convergence of 1-Fidelity Plot
    lin_xo=linspace(min(xo{i}),max(xo{i}),length(xo{i}));
    dxo=lin_xo(2)-lin_xo(1);
    lin_xo=lin_xo-dxo/2;
    hold(ax2,'off');
    caxis([0 max(max(order_hist{i}))])
    imagesc(ax2,1:(steps-1),lin_xo(end:-1:1),order_hist{i}(end:-1:1,:),'HandleVisibility','off');
    colormap(je)

    set(ax2,'YDir','normal')
    hold(ax2,'on');
    plot3(ax2,1:(steps-1),Order_f1_av{i},ones(size(Order_f1_av{i})),'k')
    plot3(ax2,1:(steps-1),Order_f1_min{i},ones(size(Order_f1_min{i})),'--k')
    plot3(ax2,1:(steps-1),Order_f1_noon{i},ones(size(Order_f1_noon{i})),'--r')
    xlabel(ax2,'Trotter steps')
    ylabel(ax2,'Order of Convergence')
    legend(ax2,{'Average Fidelity','Minimum Fidelity','NOON Fidelity'},'Location','northeast')
    axis(ax2,[1 steps-1 min(xo{i}) max(xo{i}) 0 1]);
    view(ax2,0,90)
    hold(ax2,'off');

    %matlab2tikz(['Fidelity_Density/order_fidelity_density_N_' num2str(N) '.tex'],'standalone',true)
end