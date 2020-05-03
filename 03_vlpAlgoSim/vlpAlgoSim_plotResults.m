clear all
close all
clc

answer = questdlg('Which simulation scenario?', 'Choose scenario', 'Localization SoA comp.','Platoon SUMO','Collision Avoidance - Lane Change', 'Localization SoA comp.');
% Handle response
switch answer
    case 'Localization SoA comp.'

load('results/vlpAlgoRun_bechadg_-45dBm_2000Hz.mat');
e_x = vlpAlgoRun.actl_x-vlpAlgoRun.calc_x;
e_y = vlpAlgoRun.actl_y-vlpAlgoRun.calc_y;

%%%% XY Grid Plot of Actual versus Estimated
figure
ca = plot(vlpAlgoRun.calc_x,vlpAlgoRun.calc_y, 'r','LineWidth',0.5);
grid on, hold on,xlabel('Ego x [m]'),ylabel('Ego y [m]'); 
ac = plot(vlpAlgoRun.actl_x,vlpAlgoRun.actl_y, 'b','LineWidth',1);
legend([ca ac],'Location','Northwest','Estimated', 'True');
xlim([-4.5 4.5]);
ylim([-1 8]);

%%% Histogram of errors in lateral and longitudinal directions
figure
subplot(1,2,1),histfit(e_x*100),title('Histogram of x errors'), xlabel('Error [cm]'), ylabel('Frequency'), legend('Location','Northwest',sprintf('std=%.2f',std(e_x*100)));
subplot(1,2,2),histfit(e_y*100),title('Histogram of y errors'), xlabel('Error [cm]'), ylabel('Frequency'), legend('Location','Northwest',sprintf('std=%.2f',std(e_y*100)));

%%% Individual x and y estimates versus actuals
figure
plot(vlpAlgoRun.vehCfg.t(1:vlpAlgoRun.DecimRate:end),vlpAlgoRun.actl_x, 'b-' , 'LineWidth',1),grid on, hold on
plot(vlpAlgoRun.vehCfg.t(1:vlpAlgoRun.DecimRate:end),vlpAlgoRun.calc_x, 'r--', 'LineWidth',0.5),grid on, hold on
plot(vlpAlgoRun.vehCfg.t(1:vlpAlgoRun.DecimRate:end),vlpAlgoRun.actl_y, 'b-' , 'LineWidth',1),grid on, hold on
plot(vlpAlgoRun.vehCfg.t(1:vlpAlgoRun.DecimRate:end),vlpAlgoRun.calc_y, 'r--', 'LineWidth',0.5),grid on, hold off
xlabel('Time [s]'), ylabel('Coordinate [m]');
legend('Location','West','X True','X Est.','Y True','Y Est.'); 
title('X/Y Estimation Performance');

    case 'Platoon SUMO'

load('results/vlpAlgoRun_kocbalik_-50dBm_500Hz.mat');
vlpAlgoRun_hiNs_hiRt = vlpAlgoRun;
hiNs = -50;
hiRt = round(1/(vlpAlgoRun.vehCfg.t_dt*vlpAlgoRun.DecimRate));
e_x_hiNs_hiRt = vlpAlgoRun_hiNs_hiRt.actl_x-vlpAlgoRun_hiNs_hiRt.calc_x;
e_y_hiNs_hiRt = vlpAlgoRun_hiNs_hiRt.actl_y-vlpAlgoRun_hiNs_hiRt.calc_y;
e_hiNs_hiRt = sqrt(e_x_hiNs_hiRt.^2 + e_y_hiNs_hiRt.^2);

load('results/vlpAlgoRun_kocbalik_-70dBm_500Hz.mat');
vlpAlgoRun_loNs_hiRt = vlpAlgoRun;
loNs = -70;
e_x_loNs_hiRt = vlpAlgoRun_loNs_hiRt.actl_x-vlpAlgoRun_loNs_hiRt.calc_x;
e_y_loNs_hiRt = vlpAlgoRun_loNs_hiRt.actl_y-vlpAlgoRun_loNs_hiRt.calc_y;
e_loNs_hiRt = sqrt(e_x_loNs_hiRt.^2 + e_y_loNs_hiRt.^2);

load('results/vlpAlgoRun_kocbalik_-50dBm_50Hz.mat');
vlpAlgoRun_hiNs_loRt = vlpAlgoRun;
loRt = round(1/(vlpAlgoRun.vehCfg.t_dt*vlpAlgoRun.DecimRate));
e_x_hiNs_loRt = vlpAlgoRun_hiNs_loRt.actl_x-vlpAlgoRun_hiNs_loRt.calc_x;
e_y_hiNs_loRt = vlpAlgoRun_hiNs_loRt.actl_y-vlpAlgoRun_hiNs_loRt.calc_y;
e_hiNs_loRt = sqrt(e_x_hiNs_loRt.^2 + e_y_hiNs_loRt.^2);

figure
ap = plot(vlpAlgoRun_hiNs_hiRt.actl_x, vlpAlgoRun_hiNs_hiRt.actl_y, 'Color', [0.4 0.9 0.9],'LineWidth',4);
grid on, hold on,xlabel('Ego x [m]'),ylabel('Ego y [m]'); 
p_hiNs_hiRt = plot(vlpAlgoRun_hiNs_hiRt.calc_x, vlpAlgoRun_hiNs_hiRt.calc_y, 'b','LineWidth',1.5);
p_loNs_hiRt = plot(vlpAlgoRun_loNs_hiRt.calc_x, vlpAlgoRun_loNs_hiRt.calc_y, 'r', 'LineWidth',1.5);
p_hiNs_loRt  = plot(vlpAlgoRun_hiNs_loRt.calc_x,  vlpAlgoRun_hiNs_loRt.calc_y,  'g','LineWidth',1.5);
legend([p_loNs_hiRt p_hiNs_hiRt p_hiNs_loRt ap],'Location','Northeast',sprintf('VLP: %dHz, %ddBm',hiRt, loNs)...
                                                                      ,sprintf('VLP: %dHz, %ddBm',hiRt, hiNs)...
                                                                      ,sprintf('VLP: %dHz, %ddBm',loRt, hiNs)...
                                                                      ,'True');
xlim([-8 8]);
ylim([0 22]);
hold off; 

figure
p_hiNs_hiRt = semilogy(vlpAlgoRun_hiNs_hiRt.vehCfg.t(1:vlpAlgoRun_hiNs_hiRt.DecimRate:end), e_hiNs_hiRt, 'b', 'LineWidth', 1);
grid on, hold on, xlabel('Time [s]'), ylabel('Error [m]'); 
p_loNs_hiRt = semilogy(vlpAlgoRun_loNs_hiRt.vehCfg.t(1:vlpAlgoRun_loNs_hiRt.DecimRate:end), e_loNs_hiRt, 'r', 'LineWidth', 1);
p_hiNs_loRt = semilogy(vlpAlgoRun_hiNs_loRt.vehCfg.t(1:vlpAlgoRun_hiNs_loRt.DecimRate:end),   e_hiNs_loRt,  'g', 'LineWidth', 1);
p10 = semilogy([0 vlpAlgoRun_hiNs_loRt.vehCfg.t(end)],[0.1 0.1], 'y--', 'LineWidth', 4);
legend([p_loNs_hiRt p_hiNs_hiRt p_hiNs_loRt p10],'Location','Northwest',sprintf('VLP: %dHz, %ddBm',hiRt, loNs)...
                                                                       ,sprintf('VLP: %dHz, %ddBm',hiRt, hiNs)...
                                                                       ,sprintf('VLP: %dHz, %ddBm',loRt, hiNs)...
                                                                       ,'10-cm Line');

    case 'Collision Avoidance - Lane Change'

load('results/vlpAlgoRun_laneChg_-40dBm_1000Hz.mat');
vlpAlgoRun_hiHz = vlpAlgoRun;
hiHz = 1/(vlpAlgoRun.vehCfg.t_dt*vlpAlgoRun.DecimRate);
e_x_hiHz = vlpAlgoRun_hiHz.actl_x-vlpAlgoRun_hiHz.calc_x;
e_y_hiHz = vlpAlgoRun_hiHz.actl_y-vlpAlgoRun_hiHz.calc_y;
e_hiHz = sqrt(e_x_hiHz.^2 + e_y_hiHz.^2);

load('results/vlpAlgoRun_laneChg_-40dBm_250Hz.mat');
vlpAlgoRun_loHz = vlpAlgoRun;
loHz = 1/(vlpAlgoRun.vehCfg.t_dt*vlpAlgoRun.DecimRate);
e_x_loHz = vlpAlgoRun_loHz.actl_x-vlpAlgoRun_loHz.calc_x;
e_y_loHz = vlpAlgoRun_loHz.actl_y-vlpAlgoRun_loHz.calc_y;
e_loHz = sqrt(e_x_loHz.^2 + e_y_loHz.^2);

%%%% XY Grid Plot of Actual versus Esimated
figure
hold on;
for i=1:length(vlpAlgoRun_hiHz.actl_x)
    xx = [vlpAlgoRun_hiHz.actl_x(i)-3*sqrt(vlpAlgoRun_hiHz.bound_var_x(i)) ...
          vlpAlgoRun_hiHz.actl_x(i)+3*sqrt(vlpAlgoRun_hiHz.bound_var_x(i)) ...
          vlpAlgoRun_hiHz.actl_x(i)+3*sqrt(vlpAlgoRun_hiHz.bound_var_x(i)) ...
          vlpAlgoRun_hiHz.actl_x(i)-3*sqrt(vlpAlgoRun_hiHz.bound_var_x(i))];
    yy = [vlpAlgoRun_hiHz.actl_y(i)-3*sqrt(vlpAlgoRun_hiHz.bound_var_y(i)) ...
          vlpAlgoRun_hiHz.actl_y(i)-3*sqrt(vlpAlgoRun_hiHz.bound_var_y(i)) ...
          vlpAlgoRun_hiHz.actl_y(i)+3*sqrt(vlpAlgoRun_hiHz.bound_var_y(i)) ...
          vlpAlgoRun_hiHz.actl_y(i)+3*sqrt(vlpAlgoRun_hiHz.bound_var_y(i))];
    p_hi = patch(xx,yy,[0.9290 0.6940 0.1250],'EdgeColor','none'); % ,'FaceAlpha',.1
end
for i=1:length(vlpAlgoRun_loHz.actl_x)
    xx = [vlpAlgoRun_loHz.actl_x(i)-3*sqrt(vlpAlgoRun_loHz.bound_var_x(i)) ...
          vlpAlgoRun_loHz.actl_x(i)+3*sqrt(vlpAlgoRun_loHz.bound_var_x(i)) ...
          vlpAlgoRun_loHz.actl_x(i)+3*sqrt(vlpAlgoRun_loHz.bound_var_x(i)) ...
          vlpAlgoRun_loHz.actl_x(i)-3*sqrt(vlpAlgoRun_loHz.bound_var_x(i))];
    yy = [vlpAlgoRun_loHz.actl_y(i)-3*sqrt(vlpAlgoRun_loHz.bound_var_y(i)) ...
          vlpAlgoRun_loHz.actl_y(i)-3*sqrt(vlpAlgoRun_loHz.bound_var_y(i)) ...
          vlpAlgoRun_loHz.actl_y(i)+3*sqrt(vlpAlgoRun_loHz.bound_var_y(i)) ...
          vlpAlgoRun_loHz.actl_y(i)+3*sqrt(vlpAlgoRun_loHz.bound_var_y(i))];
    p_lo = patch(xx,yy,[0.3010 0.7450 0.9330],'EdgeColor','none'); % ,'FaceAlpha',.1
end

ca2 = plot(vlpAlgoRun_hiHz.calc_x,vlpAlgoRun_hiHz.calc_y, 'g','LineWidth',2);
ca1 = plot(vlpAlgoRun_loHz.calc_x,vlpAlgoRun_loHz.calc_y, 'r','LineWidth',2);
grid on, hold on,xlabel('Ego x [m]'),ylabel('Ego y [m]'); 
ac = plot(vlpAlgoRun_loHz.actl_x,vlpAlgoRun_loHz.actl_y, 'b','LineWidth',2);
legend([ca1 ca2 ac p_lo p_hi],'Location','Southwest',sprintf('Estimated, %dHz',loHz)...
                                                    ,sprintf('Estimated, %dHz',hiHz)...
                                                    , 'True'...
                                                    , sprintf('Bound, %dHz',loHz)...
                                                    , sprintf('Bound, %dHz',hiHz));
xlim([-5 5]);
xticks = [-5 -4 -3 -2 -1 0 1 2 3 4 5];
set(gca,'XTick',xticks)
ylim([-1 9]);
hold off; 

%%% Filter results a bit, they're too erratic, no additional insight gained
%%% by looking at full rate errors, we already have the variance by bounds...
fl = 10;
e_loHz_f = movmean(e_loHz,fl);
e_hiHz_f = movmean(e_hiHz,fl);

var_lo = sqrt(3*sqrt(vlpAlgoRun_loHz.bound_var_x).^2 + 3*sqrt(vlpAlgoRun_loHz.bound_var_y).^2);
var_hi = sqrt(3*sqrt(vlpAlgoRun_hiHz.bound_var_x).^2 + 3*sqrt(vlpAlgoRun_hiHz.bound_var_y).^2);
figure
plo = semilogy(vlpAlgoRun_loHz.vehCfg.t(1:vlpAlgoRun_loHz.DecimRate:end), e_loHz_f, 'r', 'LineWidth', 2);
grid on, hold on, xlabel('Time [s]'), ylabel('Error [m]'); 
phi = semilogy(vlpAlgoRun_hiHz.vehCfg.t(1:vlpAlgoRun_hiHz.DecimRate:end), e_hiHz_f, 'g', 'LineWidth', 2);
bhi = semilogy(vlpAlgoRun_hiHz.vehCfg.t(1:vlpAlgoRun_hiHz.DecimRate:end), var_hi, 'Color', [0.9290 0.6940 0.1250], 'LineWidth', 3);
blo = semilogy(vlpAlgoRun_loHz.vehCfg.t(1:vlpAlgoRun_loHz.DecimRate:end), var_lo, 'Color', [0.3010 0.7450 0.9330], 'LineWidth', 3);
p_10 = semilogy([0 vlpAlgoRun_hiHz.vehCfg.t(end)],[0.1 0.1], 'y--', 'LineWidth', 4);
legend([plo phi blo bhi p_10],'Location','Southwest',sprintf('Error, %dHz',loHz)...
                                                        ,sprintf('Error, %dHz',hiHz)...
                                                        ,sprintf('Bound, %dHz',loHz)...
                                                        ,sprintf('Bound, %dHz',hiHz)...
                                                        ,'10-cm Line');
                                                    
%%% We've used the magnifyOnFigure package on top of these figures. Source:
%%% https://www.mathworks.com/matlabcentral/fileexchange/26007-on-figure-magnifier
                                                    
end

