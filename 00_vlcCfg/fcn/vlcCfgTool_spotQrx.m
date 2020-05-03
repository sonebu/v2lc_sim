function [ theta_l, theta_r ] = vlcCfgTool_spotQrx( handles )
%VLCCFGTOOL_SPOTQRX Summary of this function goes here
%   Detailed explanation goes here

%%% Plot the spot on the quad rectangles, 

dDL = handles.dDL;
dSL = handles.dSL;
dFL = handles.dFL;
dPX = handles.dPX;

dr = (dDL/2)*(dFL-dSL)/dFL;

lat_l = handles.lat -1;
fwd_l = handles.fwd;

lat_r = handles.lat + 1;
fwd_r = handles.fwd;

theta_l = atan2(lat_l,fwd_l);
theta_r = atan2(lat_r,fwd_r);

qpd_x_l = tan(theta_l)*dSL;
qpd_x_r = tan(theta_r)*dSL;
cla(handles.spotQrxAxis)
%%% visuals
hc = [-dPX -dPX dPX     dPX];
rectangle(handles.spotQrxAxis,'Position',hc,'LineWidth',1);
ha = [-dPX  0   dPX     dPX];
rectangle(handles.spotQrxAxis,'Position',ha,'LineWidth',1);
hd = [0    -dPX dPX     dPX];
rectangle(handles.spotQrxAxis,'Position',hd,'LineWidth',1);
hb = [0     0   dPX     dPX];
rectangle(handles.spotQrxAxis,'Position',hb,'LineWidth',1);
h_circ_l = viscircles(handles.spotQrxAxis,[0+qpd_x_l 0],dr,'Color','b','LineWidth',1);
h_circ_r = viscircles(handles.spotQrxAxis,[0+qpd_x_r 0],dr,'Color','r','LineWidth',1);
axis(handles.spotQrxAxis,[-dPX*1.5 dPX*1.5 -dPX*1.5 dPX*1.5])
clear ha hb hc hd h_circ
set(handles.spotQrxAxis,'YTickLabel',[]);
set(handles.spotQrxAxis,'XTickLabel',[]);

end

