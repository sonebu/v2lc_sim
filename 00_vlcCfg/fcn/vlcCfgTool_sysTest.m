function [ output_args ] = vlcCfgTool_sysTest( handles )
%VLCCFGTOOL_SYSTEST Summary of this function goes here
%   Detailed explanation goes here
h = handles;

%%% Plot the ego/target test configuration. 

cla(h.sysTestAxis);
plot(h.sysTestAxis,-1,0,'r*','LineWidth',5) %qrx1
hold(h.sysTestAxis,'on');
plot(h.sysTestAxis, 1,0,'b*','LineWidth',5) %qrx2
plot(h.sysTestAxis,[0 0],[-200 500],'k');
plot(h.sysTestAxis,[-1 h.lat],[0 h.fwd],'k--');
plot(h.sysTestAxis,[1 h.lat],[0 h.fwd],'k--');
plot(h.sysTestAxis,h.lat,h.fwd,'g*','LineWidth',5); %tx
quiver(h.sysTestAxis,h.lat,h.fwd,cosd(h.hdg-90)*4,sind(h.hdg-90)*4,'g','MaxHeadSize',3)
text(h.sysTestAxis,h.lat,h.fwd,'  \leftarrow TX');
text(h.sysTestAxis,1,0,'  \leftarrow QRX 1&2');
xlim(h.sysTestAxis,[-12 12])
ylim(h.sysTestAxis,[-1 35])
grid(h.sysTestAxis,'on');
hold(h.sysTestAxis,'off');

%%% This part became pretty obsolete, we've updated the ray optics
%%% simulations according to the new rectangle approx, dropped Erdem's
%%% simit.
%%% @self: check drawings at Sariyer
% dir_vect = [-cosd(h.hdg-90);-sind(h.hdg-90)];
% dir_vect = dir_vect/(norm(dir_vect));
% norm_dir_vect = [dir_vect(2);-dir_vect(1)];
% 
% P2_t0 = [h.lat h.fwd];
% P2_prime_t0 = P2_t0 - dir_vect'*(P2_t0*dir_vect);
% P1_t0 = P2_prime_t0 - norm_dir_vect'*(norm(P2_t0-P2_prime_t0)/tan(deg2rad(90)-deg2rad(120)));
% P3_t0 = P2_prime_t0 + norm_dir_vect'*(norm(P2_t0-P2_prime_t0)/tan(deg2rad(90)-deg2rad(120)));
% patch([[P1_t0(1) P2_t0(1) P3_t0(1)] [P1_t0(1) P3_t0(1)]], [[P1_t0(2) P2_t0(2) P3_t0(2)] [P1_t0(2) P3_t0(2)]],'y--','FaceAlpha',0.1);

end

