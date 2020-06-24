function [n_F_u F_u]=pretty_histogram( F,F_e_u,F_g,n_bins,legender,title_str )

f=figure();
h=histogram(F,n_bins,'BinLimits',[0 1],'EdgeColor','none','Normalization','probability');
axis tight;
ax=gca;
x_lim=ax.XLim;
y_lim=ax.YLim;
binE=h.BinEdges;
binV=h.Values;
clear h;

x=[binE(1)];
y=[0];
for i=1:length(binV)
    x=[x binE(i) binE(i+1)];
    y=[y binV(i) binV(i)];
end
x=[x binE(end)];
y=[y 0];

fill(x,y,[ 1.0000 0.5490 0]);
axis tight;
hold on;    
if length(F_e_u)>0
    f_e_u=plot(F_e_u*[1 1],[y_lim(1) y_lim(2)*1.125],'r--');
end
if length(F_g)>0
    f_g=plot(F_g*[1 1],[y_lim(1) y_lim(2)*1.125],'b--');
end
y_lim2=[y_lim(1) y_lim(2)*1.05];
axis([0 1 y_lim(1) y_lim(2)*1.125]);
xlabel('$F_1$');
ylabel('$p(F_1)$');
set(gca,'YTick',[])
set(gca,'YTickLabel',[])
if nargin==6
    title(title_str);
end

F_u=uniquetol(F,10^-14);
n_F_u=length(F_u);
if nargin>4
    if 1==legender(1)
        if length(F_e_u)>0 && length(F_g)>0
            legend([f_e_u f_g],{'Even-Uneven Ordering','Gray Ordering'});
        end
    end
    if length(legender)>1
        fy=find(y);
        fx=x(fy);
        dx=max(diff(x))/2;
        minimax=[min(fx)-dx max(fx)+dx];
        filler=fill([minimax minimax([2 1 1])],[y_lim(1) y_lim(1) y_lim(2)*1.05 y_lim(2)*1.05 y_lim(1)],0.5*[1 1 1]);
        set(filler,'facealpha',.25);
        set(filler,'EdgeColor','none');
        plot([minimax minimax([2 1 1])],[y_lim(1) y_lim(1) y_lim(2)*1.05 y_lim(2)*1.05 y_lim(1)],':k');
       
        %Make new axis
        % left or right
        [~,lr]=min([minimax(1),1-minimax(2)]);
        ax1=gca;
        pos=get(ax1,'Position'); %xmin ymin dx dy
        dpx=pos(3);
        dpy=pos(4);
        if lr==2    %(left)
            pos2=[pos(1)+0.025*dpx , pos(2)+dpy*0.3 , dpx*(0.4) , dpy*(0.675)];
        else         %(right)
            pos2=[pos(1)+0.575*dpx , pos(2)+dpy*0.3 , dpx*(0.4) , dpy*(0.675)];
        end
        ax2=axes('Position',pos2,'Color','none');
        h2=histogram(ax2,F,round(n_bins/2),'BinLimits',minimax,'EdgeColor','none','Normalization','probability');
        x_lim=ax2.XLim;
        y_lim=ax2.YLim;
        binE=h2.BinEdges;
        binV=h2.Values;
        clear h2;

        x2=[binE(1)];
        y2=[0];
        for i=1:length(binV)
            x2=[x2 binE(i) binE(i+1)];
            y2=[y2 binV(i) binV(i)];
        end
        x2=[x2 binE(end)];
        y2=[y2 0];

        fill(ax2,x2,y2,[ 1.0000 0.5490 0]);
        axis tight;
        hold on;   
        set(ax2,'Color','none');
        set(ax2,'box','off');
        set(ax2,'YTick',[]);
        set(ax2,'YTickLabel',[]);
        set(ax2,'XTick',[]);
        set(ax2,'XTickLabel',[]);
        filler2=fill(ax2,minimax([1 2 2 1 1]),[y_lim(1) y_lim(1) y_lim(2)*1.05 y_lim(2)*1.05 y_lim(1)],0.5*[1 1 1]);
        set(filler2,'facealpha',.25);
        set(filler2,'EdgeColor','none');
        plot(minimax([1 2 2 1 1]),[y_lim(1) y_lim(1) y_lim(2)*1.05 y_lim(2)*1.05 y_lim(1)],':k');
        if length(F_e_u)>0
            f_e_u=plot(ax2,F_e_u*[1 1],[y_lim(1) y_lim(2)*1.05],'r--');
        end
        if length(F_g)>0
            f_g=plot(ax2,F_g*[1 1],[y_lim(1) y_lim(2)*1.05],'b--');
        end
        axis([minimax(1) minimax(2) y_lim(1) y_lim(2)*1.05]);

        %Connect plots
        p1=ax.Position;
        p1x=ax.XLim;
        p1x(2)=p1x(2)-p1x(1);
        p1y=ax.YLim;
        p1y(2)=p1y(2)-p1y(1); %piy=[pymin dpy];
        p2=ax2.Position;

        p2_xl=(p2(1)-p1(1))/p1(3);
        p2_yl=(p2(2)-p1(2))/p1(4);
        p2_xr=(p2(1)+p2(3)-p1(1))/p1(3);
        p2_yr=(p2(2)+p2(4)-p1(2))/p1(4);
        xl=p1x(1)+p2_xl*p1x(2);
        yl=p1y(1)+p2_yl*p1y(2);
        xr=p1x(1)+p2_xr*p1x(2);
        yr=p1y(1)+p2_yr*p1y(2);

        if lr==2
            plot(ax1,[minimax(1) xl],[0 yl],':k');
            plot(ax1,[minimax(2) xr],[y_lim2(2) yr],':k');
        else
            plot(ax1,[minimax(1) xl],[y_lim2(2) yr],':k');
            plot(ax1,[minimax(2) xr],[0 yl],':k');
        end
    end
end
hold off;
end

