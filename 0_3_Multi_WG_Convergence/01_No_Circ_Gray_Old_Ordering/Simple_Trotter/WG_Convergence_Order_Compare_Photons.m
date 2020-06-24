%Convergence_Order_Compare_Trotter.m
function h=Convergence_Order_Compare_Trotter( )
clear all;
close all;
clc;

mode=4;

p=pwd;
if any(strfind(p,'\'));
    elem=strsplit(p,'\');
else
    elem=strsplit(p,'/');
end
shortened=fullfile(elem{1:end-3});
addpath(genpath(shortened));

%Photon number
n_m=1:7;
n_m2=1:3;

%Trotter Iteration number
steps=10;

%Waveguide number
N_m=2:4;
N_m2=4:6; %Loading only depends on maximum

%Iteration time
t=pi/4;


name=['trot_compare_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_N_' num2str(max(N_m)) '_s_' num2str(steps) '_t_' num2str(t) '_pauli_gray'];
name2=['trot_compare_n_' num2str(min(n_m2)) '_to_' num2str(max(n_m2)) '_N_' num2str(max(N_m2)) '_s_' num2str(steps) '_t_' num2str(t) '_pauli_gray'];

a=load([name '.mat']);
a2=load([name2 '.mat']);

F_av_ind=a.F_av_ind;
F_av_ind2=a2.F_av_ind;

%% Plot Convergence 
% For Waveguides
if mode==1  %Evolution of Error (1-F)
    trotts=10; %Initial Plot
    
    h=figure('Name','Error Rates','Position',[100,100,1000,600]);
    set(h,'color',0.98*[1 1 1]);
    ax=subplot('Position',[0.075 0.05 0.8 0.9]);
    
    c=[];
    c2=[];
    %(#waveguides , #photons , #Trotter steps)
    c=1-F_av_ind;
    c2=1-F_av_ind2;
    c2=c2(3:end,:,:);
    
    trotter_slider=uicontrol('Parent',h,'Style','slider','Tag','trott_slide','Position',[925 150 20 310],'Min',1,'Max',10,'Value',steps,'SliderStep',1./((steps-1)*[1 1]));
    trotter_text=uicontrol('Parent',h,'Style','text','Tag','trotter_text','Position',[885 90 100 50],'String',['#T=' num2str(trotts)],'BackgroundColor',get(h,'Color'),'HorizontalAlignment','center');
    addlistener(trotter_slider,'ContinuousValueChange',@(src,evt)Redraw(src,evt,ax,h,n_m,N_m,n_m2,N_m2,c,c2));

    
    s1=surf(ax,n_m,N_m,c(:,:,trotts));
    hold on;
    s2=surf(ax,n_m2,N_m2,c2(:,:,trotts));
    
    ax.YTick=[min([N_m N_m2]):max([N_m N_m2])];
    ax.XTick=[min([n_m n_m2]):max([n_m n_m2])];

    view(-45,60)
    ylabel('# WGs');
    xlabel('# Photons');
    zlabel('1-F');
    axis tight;
    drawnow();
    hold off;
    
elseif mode==2
    
    trotts=10; %Initial Plot
    
    h=figure('Name','Error Rates','Position',[100,100,1000,600]);
    set(h,'color',0.98*[1 1 1]);
    ax=subplot(1,1,1);
    
    c=[];
    c2=[];
    %(#waveguides , #photons , #Trotter steps)
    c=1-F_av_ind;
    c2=1-F_av_ind2;
    c2=c2(3:end,:,:);
    
    %trotter_slider=uicontrol('Parent',h,'Style','slider','Tag','trott_slide','Position',[925 150 20 310],'Min',1,'Max',10,'Value',trotts,'SliderStep',1./[9 9]);
    %trotter_text=uicontrol('Parent',h,'Style','text','Tag','trotter_text','Position',[885 90 100 50],'String',['#T=' num2str(trotts)],'BackgroundColor',get(h,'Color'),'HorizontalAlignment','center');
    %addlistener(trotter_slider,'ContinuousValueChange',@(src,evt)Redraw(src,evt,ax,h,n_m,N_m,n_m2,N_m2,c,c2));
    
    s1=surf(ax,n_m,N_m,c(:,:,trotts));
    hold on;
    s2=surf(ax,n_m2,N_m2,c2(:,:,trotts));
    
    ax.YTick=[min([N_m N_m2]):max([N_m N_m2])];
    ax.XTick=[min([n_m n_m2]):max([n_m n_m2])];

    axis([min([n_m n_m2]) max([n_m n_m2]) min([N_m N_m2]) max([N_m N_m2]) 0 max([c(:);c2(:)])])
    view(-45,60)
    ylabel('# WGs');
    xlabel('# Photons');
    zlabel('1-F');
    axis tight;
    drawnow();
    hold off;
    
elseif mode==3 %Convergence_Order in Number of Waveguides
    
    h_c_wg=figure('Position',[100,100,1400,600]);
    c=[];
    c2=[];
    for l=1:steps
        for j=2:size(F_av_ind,2)
            %(#waveguides , #photons , #Trotter steps)
            c(:,j-1,l)=Order_Of_Convergence(1-F_av_ind(1:end,j,l),N_m+0.5);
        end
        for j=2:size(F_av_ind2,2)
            c2(:,j-1,l)=Order_Of_Convergence(1-F_av_ind2(1:end,j,l),(2:6)+0.5);
        end
    end

    for i=1:size(c,1)
        for j=1:size(c,2)
            surf(j+1+(0:steps-1)*1/9,[i i+1]+2,reshape(c(i:i+1,j,:),2,steps));
            hold on;
        end
    end
    for i=2:size(c2,1)
        for j=1:size(c2,2)
            surf(j+1+(0:steps-1)*1/9,[i i+1]+2,reshape(c2(i:i+1,j,:),2,steps));
            hold on;
        end
    end

    xlabel('Number of photons')
    zlabel({'Convergence Order','\fontsize{9}(# of Waveguides)'})
    ylabel('Waveguides')

    ax=get(h_c_wg,'children');
    ax.YTick=[3:6];
    ax.XTick=[1:7];

    z=c(2,4,1)-c(2,4,10);

    quiver3(5,4.05,0,1,0,-z,'k')

    te=text(5.5,4.15,0,'Trotter steps','HorizontalAlignment','center','FontSize',10);
    set(te,'Rotation',19)

    view(-45,60)
    ylabel('# WGs');
    xlabel('# Photons');
    zlabel('Order of Convergence (# WGs)');
    axis tight;
    set(h_c_wg,'color',0.98*[1 1 1])
    drawnow();

    A=getframe(h_c_wg);
    imwrite(A.cdata,'wg_convergence_compare.png');
    matlab2tikz('Order_of_Convergence_Number_of_Photons.tex','standalone',true)
elseif mode==4 %Convergence_Order in Number of Photons
    
    h_c_wg=figure('Position',[100,100,1400,600]);
    
    F_av_ind2=F_av_ind2(3:end,:,:);
    c=[];
    for l=1:steps
        for j=2:size(F_av_ind,1)
            c(j-1,:,l)=Order_Of_Convergence(1-F_av_ind(j,1:end,l),(2:7)+0.5);
        end
    end
    c2=[];
    for l=1:steps
        for j=2:size(F_av_ind2,1)
            c2(j-1,:,l)=Order_Of_Convergence(1-F_av_ind2(j,1:end,l),(2:4)+0.5);
        end
    end

    c=permute(c,[3,1,2]);
    c2=permute(c2,[3,1,2]);
    for i=1:size(c,2)
        for j=size(c,3)-1:-1:1
            surf([j+1 j+2],i+1+(0:steps-1)*1/9,reshape(c(:,i,j:j+1),steps,2));
            hold on;
        end
    end
    for i=1:size(c2,2)
        for j=1:size(c2,3)-1
            surf([j+1 j+2],i+3+(0:steps-1)*1/9,reshape(c2(:,i,j:j+1),steps,2));
            hold on;
        end
    end


    ax=get(h_c_wg,'children');
    ax.YTick=[3:6];
    ax.XTick=[1:7];

    z=c(1,2,4)-c(10,2,4);

    %quiver3(5,4.05,0,1,0,-z,'k')

    %te=text(5.5,4.15,0,'Trotter steps','HorizontalAlignment','center','FontSize',10);
    %set(te,'Rotation',19)

    view(-45,60)
    ylabel('# WGs');
    xlabel('# Photons');
    zlabel('Order of Convergence (# Photons)');
    axis tight;
    set(h_c_wg,'color',0.98*[1 1 1])
    drawnow();

    A=getframe(h_c_wg);
    imwrite(A.cdata,'wg_convergence_compare.png');
    matlab2tikz('Order_of_Convergence_Number_of_Photons.tex','standalone',true)
end
end

function Redraw(a,b,ax,h,n_m,N_m,n_m2,N_m2,c,c2)
    trotter_slider=findobj(h,'Tag','trott_slide');
    trotter_text=findobj(h,'Tag','trotter_text');

    %Extract values from uicontrol elements
    slide_num=round(get(trotter_slider,'Value'));   %Boundaries of the Plot 
    set(trotter_text,'String',['#T=' num2str(slide_num)]);
    
    s1=surf(ax,n_m,N_m,c(:,:,slide_num));
    hold on;
    s2=surf(ax,n_m2,N_m2,c2(:,:,slide_num));
    
    ax.YTick=[min([N_m N_m2]):max([N_m N_m2])];
    ax.XTick=[min([n_m n_m2]):max([n_m n_m2])];

    view(-45,60)
    ylabel('# WGs');
    xlabel('# Photons');
    zlabel('1-F');
    axis tight;
    axis([min([n_m n_m2]) max([n_m n_m2]) min([N_m N_m2]) max([N_m N_m2]) 0 max([c(:);c2(:)])])
    set(h,'color',0.98*[1 1 1])
    drawnow();
    hold off;
end