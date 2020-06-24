function [ f ] = WG_Plot2PColor_Compare( ket_string,ket_string2,phi_i_red,phi_a_red,ts,ta,interp,horz_lines,how_often)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
if nargin==6
    horz_lines=0;
    how_often=10;
end
if nargin==7
    how_often=10;
end

f=figure('Position', [10 10 1300 400]);
ax1=subplot(1,2,1);
[surf_mat xi ti dx]=WG_Plot(phi_i_red,ts,1,10,1);
h=pcolor(ti,xi,surf_mat);
set(h, 'EdgeColor', 'none');
hold on;
x_lim=xlim(gca);
for i=1:length(ket_string)-1
    a=patch([x_lim x_lim(end:-1:1)],[[1 1]*(i+dx) [1 1]*(i+1-dx)], 'k');
    a.FaceAlpha = 0.4;
    a.EdgeColor='None';
    plot(x_lim,[1 1]*(i-dx),'Color',[1 1 1]*0.5);
    plot(x_lim,[1 1]*(i+dx),'Color',[1 1 1]*0.5);
end
a=patch([x_lim x_lim(end:-1:1)],[[1 1]*(0) [1 1]*(1-dx)], 'k');
a.FaceAlpha = 0.25;
a.EdgeColor='None';
a=patch([x_lim x_lim(end:-1:1)],[[1 1]*(i+1+dx) [1 1]*(i+1.5)], 'k');
a.FaceAlpha = 0.25;
a.EdgeColor='None';

ax1.YTick=[1:length(ket_string)];
ax1.YTickLabel=ket_string;

%% Second Plot
ax2=subplot(1,2,2);
[surf_mat xi ti dx]=WG_Plot(phi_a_red,ta,interp,how_often,1);
h=pcolor(ti,xi,surf_mat);
set(h, 'EdgeColor', 'none');
hold on;
x_lim=xlim(gca);
for i=1:length(ket_string)-1
    a=patch([x_lim x_lim(end:-1:1)],[[1 1]*(i+dx) [1 1]*(i+1-dx)], 'k');
    a.FaceAlpha = 0.4;
    a.EdgeColor='None';
    plot(x_lim,[1 1]*(i-dx),'Color',[1 1 1]*0.5);
    plot(x_lim,[1 1]*(i+dx),'Color',[1 1 1]*0.5);
end
a=patch([x_lim x_lim(end:-1:1)],[[1 1]*(0) [1 1]*(1-dx)], 'k');
a.FaceAlpha = 0.25;
a.EdgeColor='None';
a=patch([x_lim x_lim(end:-1:1)],[[1 1]*(i+1+dx) [1 1]*(i+1.5)], 'k');
a.FaceAlpha = 0.25;
a.EdgeColor='None';

if horz_lines==1
    y_lim=ylim(gca);
    for i=2:length(ta)-1
        plot([1 1]*ta(i),ylim,'Color',[1 1 1]*0.5);
    end
end
if length(ket_string)==0
    ax2.YTick=[];
else
    ax2.YTick=[1:length(ket_string2)];
    ax2.YTickLabel=ket_string2;
end
%cb=colorbar();
%ylabel(cb,'p')
set(gcf,'color',0.98*[1 1 1])
end

