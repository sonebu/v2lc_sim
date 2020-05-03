clear all; close all; clc

%% Inputs

% Simulation Parameters
t_stop  = 1;    % seconds
t_dt    = 1e-4; % seconds
vehCfg_t = [t_dt:t_dt:t_stop]; % seconds

% Define waypoints, rxv absolute, txv relative to rxv, waypoints are for
% front center bumper, vehicle dimensions are decided later, in 02_v2lcDataGen
vehCfg_tgtV_x_rel = [0.1458  3.5];
vehCfg_tgtV_y_rel = [ 11 11];
vehCfg_egoV_x = [0 0];
vehCfg_egoV_y = [0 120];

% Interpolate method for waypoints
% be careful with this, pchip can generate physically meaningless 
% vehicle movements, check with test.avi at v2lcDataGen.
interp_method = 'linear'; 

%% Simulation. 

% This part is basically SE(2) operations to ensure that the vehicles go
% through the waypoints. Trivial. The extra rotations and translations are 
% used for matching SUMO's weird heading representation.
% Vehicle dimensions are decided later, in 02_v2lcDataGen, this simulates
% the center front bumper position and orientation 

tmp = linspace(1,length(vehCfg_tgtV_x_rel),length(vehCfg_t));
vehCfg_tgtV_x_rel = interp1(1:length(vehCfg_tgtV_x_rel),vehCfg_tgtV_x_rel,tmp,interp_method);
vehCfg_tgtV_y_rel = interp1(1:length(vehCfg_tgtV_y_rel),vehCfg_tgtV_y_rel,tmp,interp_method);
tmp = linspace(1,length(vehCfg_egoV_x),length(vehCfg_t));
vehCfg_egoV_x = interp1(1:length(vehCfg_egoV_x),vehCfg_egoV_x,tmp,interp_method);
vehCfg_egoV_y = interp1(1:length(vehCfg_egoV_y),vehCfg_egoV_y,tmp,interp_method);

vehCfg_egoV_hdg = zeros(1,length(vehCfg_egoV_x));
for i=2:length(vehCfg_egoV_hdg)-1
    vehCfg_egoV_hdg(i) = rad2deg( (atan2(vehCfg_egoV_x(i)-vehCfg_egoV_x(i-1),vehCfg_egoV_y(i)-vehCfg_egoV_y(i-1)) + atan2(vehCfg_egoV_x(i+1)-vehCfg_egoV_x(i),vehCfg_egoV_y(i+1)-vehCfg_egoV_y(i)))/2 );
end
vehCfg_egoV_hdg(1) = vehCfg_egoV_hdg(2); vehCfg_egoV_hdg(end) = vehCfg_egoV_hdg(end-1);

rot_mtx = zeros(2,2,length(vehCfg_t));
for i=1:length(vehCfg_t)
    rot_mtx(1,1,i) =  cosd(-vehCfg_egoV_hdg(i));
	rot_mtx(1,2,i) = -sind(-vehCfg_egoV_hdg(i));
    rot_mtx(2,1,i) =  sind(-vehCfg_egoV_hdg(i));
	rot_mtx(2,2,i) =  cosd(-vehCfg_egoV_hdg(i));
    rel = rot_mtx(:,:,i)*[vehCfg_tgtV_x_rel(i);vehCfg_tgtV_y_rel(i)];
    vehCfg_tgtV_x_rel(i) = rel(1);
    vehCfg_tgtV_y_rel(i) = rel(2);
end
clear i rot_mtx rel

vehCfg_tgtV_x = vehCfg_egoV_x + vehCfg_tgtV_x_rel;
vehCfg_tgtV_y = vehCfg_egoV_y + vehCfg_tgtV_y_rel;

clear tmp

vehCfg_tgtV_hdg = zeros(1,length(vehCfg_tgtV_x));
for i=2:length(vehCfg_egoV_hdg)-1
    vehCfg_tgtV_hdg(i) = rad2deg( (atan2(vehCfg_tgtV_x(i)-vehCfg_tgtV_x(i-1),vehCfg_tgtV_y(i)-vehCfg_tgtV_y(i-1)) + atan2(vehCfg_tgtV_x(i+1)-vehCfg_tgtV_x(i),vehCfg_tgtV_y(i+1)-vehCfg_tgtV_y(i)))/2 );
end
vehCfg_tgtV_hdg(1) = vehCfg_tgtV_hdg(2); vehCfg_tgtV_hdg(end) = vehCfg_tgtV_hdg(end-1);

vehCfg_egoV_spd = zeros(1,length(vehCfg_egoV_x));
vehCfg_tgtV_spd = zeros(1,length(vehCfg_tgtV_x));
for i=2:length(vehCfg_egoV_spd)-1
    vehCfg_egoV_spd(i) = (sqrt((vehCfg_egoV_x(i)-vehCfg_egoV_x(i-1))*(vehCfg_egoV_x(i)-vehCfg_egoV_x(i-1)) + (vehCfg_egoV_y(i)-vehCfg_egoV_y(i-1))*(vehCfg_egoV_y(i)-vehCfg_egoV_y(i-1))) +...
                          sqrt((vehCfg_egoV_x(i+1)-vehCfg_egoV_x(i))*(vehCfg_egoV_x(i+1)-vehCfg_egoV_x(i)) + (vehCfg_egoV_y(i+1)-vehCfg_egoV_y(i))*(vehCfg_egoV_y(i+1)-vehCfg_egoV_y(i))) )/(2*t_dt);
    vehCfg_tgtV_spd(i) = (sqrt((vehCfg_tgtV_x(i)-vehCfg_tgtV_x(i-1))*(vehCfg_tgtV_x(i)-vehCfg_tgtV_x(i-1)) + (vehCfg_tgtV_y(i)-vehCfg_tgtV_y(i-1))*(vehCfg_tgtV_y(i)-vehCfg_tgtV_y(i-1))) +...
                          sqrt((vehCfg_tgtV_x(i+1)-vehCfg_tgtV_x(i))*(vehCfg_tgtV_x(i+1)-vehCfg_tgtV_x(i)) + (vehCfg_tgtV_y(i+1)-vehCfg_tgtV_y(i))*(vehCfg_tgtV_y(i+1)-vehCfg_tgtV_y(i))) )/(2*t_dt);
end
vehCfg_egoV_spd(1) = vehCfg_egoV_spd(2); vehCfg_egoV_spd(end) = vehCfg_egoV_spd(end-1);
vehCfg_tgtV_spd(1) = vehCfg_tgtV_spd(2); vehCfg_tgtV_spd(end) = vehCfg_tgtV_spd(end-1);

vehCfg.tgt.x = vehCfg_tgtV_x(:);
vehCfg.tgt.y = vehCfg_tgtV_y(:);
vehCfg.tgt.spd = vehCfg_tgtV_spd(:);
vehCfg.tgt.hdg = vehCfg_tgtV_hdg(:);
vehCfg.ego.x = vehCfg_egoV_x(:);
vehCfg.ego.y = vehCfg_egoV_y(:);
vehCfg.ego.spd = vehCfg_egoV_spd(:);
vehCfg.ego.hdg = vehCfg_egoV_hdg(:);
vehCfg.t = vehCfg_t(:);

answer = inputdlg('Enter filename for the vehicular configuration file','Vehicular Configuration Filename',[1 50],{'vehCfg_<explanation>.mat'});
save(strcat('data/',answer{1}),'vehCfg')

