function [ f ] = WG_Plot2PColor( ket_string,phi_i_red,ts )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
f=figure();
[surf_mat xi ti dx]=WG_Plot(phi_i_red,ts,1,10,1);
h=pcolor(ti,xi,surf_mat);
set(h, 'EdgeColor', 'none');
hold on;
x_lim=xlim(gca);
for i=1:length(ket_string)-1
    a=patch([x_lim x_lim(end:-1:1)],[[1 1]*(i+dx) [1 1]*(i+1-dx)], 'k');
    a.FaceAlpha = 0.25;
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
end

