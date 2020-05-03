clear all
close all
clc

% choose sumo simulation csv after running run_<simName>.bat
[fn, pt] = uigetfile(fullfile(pwd,'sumo_ws','*.csv'),'Choose SUMO Simulation CSV');

% choose vehicle IDs inside that simulation as ego and target, this
% probably requires manually watching the simulation and choosing vehicles
prompt = {'Enter ego vehicle ID:','Enter target vehicle ID:'};
dlgtitle = 'Ego and Target vehicle IDs from SUMO run';
dims = [1 30];
definput = {'veh13','veh12'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
clearvars -except answer fn pt

% read CSV, start testing for mishaps
T = readtable(strcat(pt,fn));

% check for columns vehicle_speed, vehicle_angle, vehicle_x, vehicle_y in csv
if((~any(strcmp(T.Properties.VariableNames,'vehicle_x')))     || ...
   (~any(strcmp(T.Properties.VariableNames,'vehicle_y')))     || ...
   (~any(strcmp(T.Properties.VariableNames,'vehicle_speed'))) || ...
   (~any(strcmp(T.Properties.VariableNames,'vehicle_angle'))) || ...
   (~any(strcmp(T.Properties.VariableNames,'vehicle_id')))    || ...
   (~any(strcmp(T.Properties.VariableNames,'timestep_time'))))
	error('The CSV you selected does not contain either one of these columns: vehicle_x, vehicle_y, vehicle_speed, vehicle_angle, vehicle_id, timestep_time')
end

% check for vehicle IDs (row) in sim in csv
id_idx = find(strcmp([T.Properties.VariableNames], 'vehicle_id')); % single line engine
if((~any(strcmp(table2cell(T(:,id_idx)),answer{1})))     || ...
   (~any(strcmp(table2cell(T(:,id_idx)),answer{2}))) )
	error('The CSV you selected does not contain data for the ego or the target vehicle you selected')
end

% Start creating ego and target arrays
t_idx = find(strcmp([T.Properties.VariableNames], 'timestep_time'));
x_idx = find(strcmp([T.Properties.VariableNames], 'vehicle_x')); 
y_idx = find(strcmp([T.Properties.VariableNames], 'vehicle_y'));
spd_idx = find(strcmp([T.Properties.VariableNames], 'vehicle_speed'));
hdg_idx = find(strcmp([T.Properties.VariableNames], 'vehicle_angle'));

ego_rows = find(strcmp(table2cell(T(:,id_idx)),answer{1}));
tgt_rows = find(strcmp(table2cell(T(:,id_idx)),answer{2}));
ego_T = sortrows(T(ego_rows,:),t_idx);
tgt_T = sortrows(T(tgt_rows,:),t_idx);
clearvars -except ego_T tgt_T t_idx id_idx x_idx y_idx spd_idx hdg_idx

ego_time = table2array(ego_T(:,t_idx));
tgt_time = table2array(tgt_T(:,t_idx));
t_dt = tgt_time(2)-tgt_time(1);
common_t_min = max(min(ego_time),min(tgt_time));
common_t_max = min(max(ego_time),max(tgt_time));

% check for mutual simulation time
if(common_t_min>common_t_max)
    error('Your ego and target vehicles are not in the simulation at the same time window');
end

% enter bounds for time window
dlgtitle = 'Simulation Time Window:';
prompt = {sprintf('Common ego-target vehicle start time: %.2f s, stop time: %.2f s\nEnter start time in absolute seconds\nEnter 0 for choosing the complete available window',common_t_min, common_t_max),sprintf('Enter start time in absolute seconds\nEnter 0 for choosing the complete available window')};
dims = [1 80];
definput = {'0','0'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
clear definput dims prompt dlgtitle
if((str2num(answer{1}) ~= 0) && (str2num(answer{1}) >= common_t_min))
    common_t_min = round(str2num(answer{1}),round(log10(1/t_dt)));
end
ego_t_minidx = find(ego_time==common_t_min);
tgt_t_minidx = find(tgt_time==common_t_min);

if((str2num(answer{2}) ~= 0) && (str2num(answer{2}) <= common_t_max))
    common_t_max = round(str2num(answer{2}),round(log10(1/t_dt)));
end
ego_t_maxidx = find(ego_time==common_t_max);
tgt_t_maxidx = find(tgt_time==common_t_max);

ego_T = ego_T(ego_t_minidx:ego_t_maxidx,:);
tgt_T = tgt_T(tgt_t_minidx:tgt_t_maxidx,:);
clearvars -except ego_T tgt_T t_idx x_idx y_idx spd_idx hdg_idx common_t_min common_t_max

vehCfg_egoV_x   = table2array(ego_T(:,x_idx));
vehCfg_egoV_y   = table2array(ego_T(:,y_idx));
vehCfg_egoV_spd = table2array(ego_T(:,spd_idx));
vehCfg_egoV_hdg = table2array(ego_T(:,hdg_idx));

vehCfg_tgtV_x   = table2array(tgt_T(:,x_idx));
vehCfg_tgtV_y   = table2array(tgt_T(:,y_idx));
vehCfg_tgtV_spd = table2array(tgt_T(:,spd_idx));
vehCfg_tgtV_hdg = table2array(tgt_T(:,hdg_idx));

vehCfg_t   = table2array(ego_T(:,t_idx));

%%% This rotation is not actually necessary but keeping 
%%% it as a template here
% vehCfg_egoV_hdg = zeros(1,length(vehCfg_egoV_x));
% vehCfg_tgtV_hdg = zeros(1,length(vehCfg_tgtV_x));
% for i=2:length(vehCfg_egoV_hdg)-1
%     vehCfg_egoV_hdg(i) = rad2deg( (atan2(vehCfg_egoV_x(i)-vehCfg_egoV_x(i-1),vehCfg_egoV_y(i)-vehCfg_egoV_y(i-1)) + atan2(vehCfg_egoV_x(i+1)-vehCfg_egoV_x(i),vehCfg_egoV_y(i+1)-vehCfg_egoV_y(i)))/2 );
% 	vehCfg_tgtV_hdg(i) = rad2deg( (atan2(vehCfg_tgtV_x(i)-vehCfg_tgtV_x(i-1),vehCfg_tgtV_y(i)-vehCfg_tgtV_y(i-1)) + atan2(vehCfg_tgtV_x(i+1)-vehCfg_tgtV_x(i),vehCfg_tgtV_y(i+1)-vehCfg_tgtV_y(i)))/2 );
% end
% vehCfg_egoV_hdg(1) = vehCfg_egoV_hdg(2); vehCfg_egoV_hdg(end) = vehCfg_egoV_hdg(end-1);
% vehCfg_tgtV_hdg(1) = vehCfg_tgtV_hdg(2); vehCfg_tgtV_hdg(end) = vehCfg_tgtV_hdg(end-1);

vehCfg.tgt.x = vehCfg_tgtV_x(:);
vehCfg.tgt.y = vehCfg_tgtV_y(:);
vehCfg.tgt.spd = vehCfg_tgtV_spd(:);
vehCfg.tgt.hdg = vehCfg_tgtV_hdg(:);
vehCfg.ego.x = vehCfg_egoV_x(:);
vehCfg.ego.y = vehCfg_egoV_y(:);
vehCfg.ego.spd = vehCfg_egoV_spd(:);
vehCfg.ego.hdg = vehCfg_egoV_hdg(:);
vehCfg.t = vehCfg_t(:)-vehCfg_t(1);

% save the configuration
answer = inputdlg('Enter filename for the vehicular configuration file','Vehicular Configuration Filename',[1 50],{'vehCfg_sumo_<explanation>.mat'});
save(strcat('data/',answer{1}),'vehCfg')
