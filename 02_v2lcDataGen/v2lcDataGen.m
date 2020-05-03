clear all;  close all;  clc;
addpath('fcn');

%% Load VLC Configuration 
[fn, pt] = uigetfile('../00_vlcCfg/data/*.mat','Load VLC Configuration');
load(strcat(pt,fn));

vlcCfg_qrx.dsl  = vlcCfg_rx_dSL;                  % meters
vlcCfg_qrx.dpx  = vlcCfg_rx_dPX;                  % meters
vlcCfg_qrx.ddl  = vlcCfg_rx_dDL;                  % meters
vlcCfg_qrx.dfl  = vlcCfg_rx_dFL;                  % meters
vlcCfg_qrx.dres = vlcCfg_rx_dRes;                 % meters
vlcCfg_qrx.rSpot = (vlcCfg_qrx.ddl/2)*(vlcCfg_qrx.dfl-vlcCfg_qrx.dsl)/vlcCfg_qrx.dfl; % meters
vlcCfg_qrx.thetaSwpRes = vlcCfg_rx_thetaSwpRes;   % degrees % not needed in data generation
vlcCfg_qrx.thetaLut = vlcCfg_rx_thetaLut;  % not needed in data generation
vlcCfg_qrx.lum2pkpk = vlcCfg_rx_lum2pkpk;         % V/lum
vlcCfg_tx.plrPatFit = vlcCfg_tx_plrPatFit;
vlcCfg_tx.plrPatNum = vlcCfg_tx_plrPatNum; % not needed in data generation
vlcCfg_tx.maxLum = vlcCfg_tx_maxLum;              % lum

%% Load Vehicular Configuration
[fn, pt] = uigetfile('../01_vehCfg/data/*.mat','Load Vehicular Configuration');
load(strcat(pt,fn));

clearvars -except vlcCfg_tx vlcCfg_qrx vehCfg

vehCfg.t_dt = vehCfg.t(4)-vehCfg.t(3);

% Rotate to fit our coordinate frame in this simulation. 
rot_mtx = zeros(2,2,length(vehCfg.t));
for i=1:length(vehCfg.t)
    rot_mtx(1,1,i) =  cosd(vehCfg.ego.hdg(i));
	rot_mtx(1,2,i) = -sind(vehCfg.ego.hdg(i));
    rot_mtx(2,1,i) =  sind(vehCfg.ego.hdg(i));
	rot_mtx(2,2,i) =  cosd(vehCfg.ego.hdg(i));
    diff = [vehCfg.tgt.x(i)-vehCfg.ego.x(i);vehCfg.tgt.y(i)-vehCfg.ego.y(i)];
    rel = rot_mtx(:,:,i)*diff;    
    vehCfg.tgt.lat_rel(i) = rel(1);
    vehCfg.tgt.fwd_rel(i) = rel(2);
end
clear i rot_mtx rel

vehCfg.tgt.hdg_rel = vehCfg.tgt.hdg - vehCfg.ego.hdg;

% Hardcoded, change if sumo configuration changes
vehCfg.ego.vehWidth = 1.6; % actual vehicle width is 1.8, lights are assumed to be 10cm inside from each end
vehCfg.tgt.vehWidth = 1.6; % actual vehicle width is 1.8, lights are assumed to be 10cm inside from each end
vehCfg.ego.vehLength = 5.0;
vehCfg.tgt.vehLength = 5.0;

%% Trajectory video, check test.avi in the same folder!!

% Comment this block if you don't want a video, this is for checking before
% running so that you don't have to wait until the end to see if you did
% something blatantly wrong..

%%%%%%% VIDEO %%%%%%%
v=VideoWriter('test.avi');
v.FrameRate = 30;
open(v);
figure
for i=round(1:((1/vehCfg.t_dt)/v.FrameRate):length(vehCfg.tgt.x))
    beta_line_x = [vehCfg.tgt.x(i)-2*sind(vehCfg.tgt.hdg(i)) vehCfg.tgt.x(i)];
    beta_line_y = [vehCfg.tgt.y(i)-2*cosd(vehCfg.tgt.hdg(i)) vehCfg.tgt.y(i)];
    alpha_line_x = [vehCfg.tgt.x(i)-2*sind(vehCfg.tgt.hdg_rel(i)) vehCfg.tgt.x(i)];
    alpha_line_y = [vehCfg.tgt.y(i)-2*cosd(vehCfg.tgt.hdg_rel(i)) vehCfg.tgt.y(i)];
    plot(vehCfg.tgt.x(1:i),vehCfg.tgt.y(1:i),'b','LineWidth',1.5),hold on
 	plot(vehCfg.ego.x(1:i),vehCfg.ego.y(1:i),'k','LineWidth',1.5),hold on
    plot(beta_line_x,beta_line_y,'g--','LineWidth',1.5), hold off
    title('Black: RX , Blue: TX, Green: Headlight Direction')
    grid on
    xl = [min(min(vehCfg.tgt.x),min(vehCfg.ego.x))-1 max(max(vehCfg.tgt.x),max(vehCfg.ego.x))+1];
    yl = [min(min(vehCfg.tgt.y),min(vehCfg.ego.y))-1 max(max(vehCfg.tgt.y),max(vehCfg.ego.y))+1];
    szx = xl(2)-xl(1);
	szy = yl(2)-yl(1);
    if(szx)>=(szy)
        xlim(xl);
        ylim([yl(1)-(szx-szy)/2 yl(2)+(szx-szy)/2]);
    else
        xlim([xl(1)-(szy-szx)/2 xl(2)+(szy-szx)/2]);
        ylim(yl);
    end
    drawnow
    frame = getframe(gcf);
    writeVideo(v,frame);
end
close(v);
%%%%%%% VIDEO %%%%%%%

answer = questdlg('The video looks OK?', 'Continue?', 'Yes','No (stops script)','No (stops script)');
% Handle response
switch answer
    case 'Yes'
    case 'No (stops script)'
        clear all;
        return;
end

%% Simulate
qrx1_gain_txR = zeros(length(vehCfg.t),1);
qrx2_gain_txR = zeros(length(vehCfg.t),1);
qrx1_gain_txL = zeros(length(vehCfg.t),1);
qrx2_gain_txL = zeros(length(vehCfg.t),1);
qrx1_gain_txR_pt = zeros(length(vehCfg.t),1);
qrx2_gain_txR_pt = zeros(length(vehCfg.t),1);
qrx1_gain_txL_pt = zeros(length(vehCfg.t),1);
qrx2_gain_txL_pt = zeros(length(vehCfg.t),1);
pwrRatio_qrx1_txR = zeros(length(vehCfg.t),1);
pwrRatio_qrx2_txR = zeros(length(vehCfg.t),1);
pwrRatio_qrx1_txL = zeros(length(vehCfg.t),1);
pwrRatio_qrx2_txL = zeros(length(vehCfg.t),1);
f = waitbar(0, 'Generating data for this trajectory');
A_sh_qrx1_txR = zeros(length(vehCfg.t),1); B_sh_qrx1_txR = zeros(length(vehCfg.t),1); C_sh_qrx1_txR = zeros(length(vehCfg.t),1); D_sh_qrx1_txR = zeros(length(vehCfg.t),1);
A_sh_qrx2_txR = zeros(length(vehCfg.t),1); B_sh_qrx2_txR = zeros(length(vehCfg.t),1); C_sh_qrx2_txR = zeros(length(vehCfg.t),1); D_sh_qrx2_txR = zeros(length(vehCfg.t),1);
A_sh_qrx1_txL = zeros(length(vehCfg.t),1); B_sh_qrx1_txL = zeros(length(vehCfg.t),1); C_sh_qrx1_txL = zeros(length(vehCfg.t),1); D_sh_qrx1_txL = zeros(length(vehCfg.t),1);
A_sh_qrx2_txL = zeros(length(vehCfg.t),1); B_sh_qrx2_txL = zeros(length(vehCfg.t),1); C_sh_qrx2_txL = zeros(length(vehCfg.t),1); D_sh_qrx2_txL = zeros(length(vehCfg.t),1);

tic
for i=1:length(vehCfg.t)
    waitbar(i/length(vehCfg.t),f, 'Generating data for this trajectory');

    % front bumper center to rear bumper center
    x_1 = vehCfg.tgt.lat_rel(i) - vehCfg.tgt.vehLength*sind(vehCfg.tgt.hdg_rel(i));
    y_1 = vehCfg.tgt.fwd_rel(i) - vehCfg.tgt.vehLength*cosd(vehCfg.tgt.hdg_rel(i));

    % rear bumper center to left tail_light
    x_2_L = x_1 - (vehCfg.tgt.vehWidth/2)*cosd(vehCfg.tgt.hdg_rel(i));
    y_2_L = y_1 + (vehCfg.tgt.vehWidth/2)*sind(vehCfg.tgt.hdg_rel(i));
    
    % rear bumper center to right tail_light
    x_2_R = x_1 + (vehCfg.tgt.vehWidth/2)*cosd(vehCfg.tgt.hdg_rel(i));
    y_2_R = y_1 - (vehCfg.tgt.vehWidth/2)*sind(vehCfg.tgt.hdg_rel(i));

    vehCfg.tgt.leftTailPosX(i) = x_2_L;
    vehCfg.tgt.leftTailPosY(i) = y_2_L;
    vehCfg.tgt.rightTailPosX(i) = x_2_R;
    vehCfg.tgt.rightTailPosY(i) = y_2_R;
    
    % Compute qrx gains (optical!)
    qrx1_gain_txR(i) = v2lcDataGen_rxGain(vlcCfg_tx.plrPatFit, x_2_R + (vehCfg.ego.vehWidth/2), y_2_R, 0, vlcCfg_qrx.ddl);
    qrx2_gain_txR(i) = v2lcDataGen_rxGain(vlcCfg_tx.plrPatFit, x_2_R - (vehCfg.ego.vehWidth/2), y_2_R, 0, vlcCfg_qrx.ddl);
    qrx1_gain_txL(i) = v2lcDataGen_rxGain(vlcCfg_tx.plrPatFit, x_2_L + (vehCfg.ego.vehWidth/2), y_2_L, 0, vlcCfg_qrx.ddl);
    qrx2_gain_txL(i) = v2lcDataGen_rxGain(vlcCfg_tx.plrPatFit, x_2_L - (vehCfg.ego.vehWidth/2), y_2_L, 0, vlcCfg_qrx.ddl);

    % Compute shares for each qrx cell, on each qrx, for each TX, so 4 runs in total
 	[ A_sh_qrx1_txR(i), B_sh_qrx1_txR(i), C_sh_qrx1_txR(i), D_sh_qrx1_txR(i) ] = v2lcDataGen_qdMeas(0, y_2_R, x_2_R + (vehCfg.ego.vehWidth/2), vlcCfg_qrx.dsl, vlcCfg_qrx.dres, vlcCfg_qrx.dpx, vlcCfg_qrx.rSpot);
 	[ A_sh_qrx2_txR(i), B_sh_qrx2_txR(i), C_sh_qrx2_txR(i), D_sh_qrx2_txR(i) ] = v2lcDataGen_qdMeas(0, y_2_R, x_2_R - (vehCfg.ego.vehWidth/2), vlcCfg_qrx.dsl, vlcCfg_qrx.dres, vlcCfg_qrx.dpx, vlcCfg_qrx.rSpot);
 	[ A_sh_qrx1_txL(i), B_sh_qrx1_txL(i), C_sh_qrx1_txL(i), D_sh_qrx1_txL(i) ] = v2lcDataGen_qdMeas(0, y_2_L, x_2_L + (vehCfg.ego.vehWidth/2), vlcCfg_qrx.dsl, vlcCfg_qrx.dres, vlcCfg_qrx.dpx, vlcCfg_qrx.rSpot);
 	[ A_sh_qrx2_txL(i), B_sh_qrx2_txL(i), C_sh_qrx2_txL(i), D_sh_qrx2_txL(i) ] = v2lcDataGen_qdMeas(0, y_2_L, x_2_L - (vehCfg.ego.vehWidth/2), vlcCfg_qrx.dsl, vlcCfg_qrx.dres, vlcCfg_qrx.dpx, vlcCfg_qrx.rSpot);
end
a = toc;
disp(['took ' num2str(a/60) ' minutes'])
clear a; close(f);
qdA_lumens_qrx1_txR = qrx1_gain_txR.*A_sh_qrx1_txR.*vlcCfg_tx.maxLum; qdB_lumens_qrx1_txR = qrx1_gain_txR.*B_sh_qrx1_txR.*vlcCfg_tx.maxLum;
qdC_lumens_qrx1_txR = qrx1_gain_txR.*C_sh_qrx1_txR.*vlcCfg_tx.maxLum; qdD_lumens_qrx1_txR = qrx1_gain_txR.*D_sh_qrx1_txR.*vlcCfg_tx.maxLum;
qdA_lumens_qrx2_txR = qrx2_gain_txR.*A_sh_qrx2_txR.*vlcCfg_tx.maxLum; qdB_lumens_qrx2_txR = qrx2_gain_txR.*B_sh_qrx2_txR.*vlcCfg_tx.maxLum;
qdC_lumens_qrx2_txR = qrx2_gain_txR.*C_sh_qrx2_txR.*vlcCfg_tx.maxLum; qdD_lumens_qrx2_txR = qrx2_gain_txR.*D_sh_qrx2_txR.*vlcCfg_tx.maxLum;
qdA_lumens_qrx1_txL = qrx1_gain_txL.*A_sh_qrx1_txL.*vlcCfg_tx.maxLum; qdB_lumens_qrx1_txL = qrx1_gain_txL.*B_sh_qrx1_txL.*vlcCfg_tx.maxLum;
qdC_lumens_qrx1_txL = qrx1_gain_txL.*C_sh_qrx1_txL.*vlcCfg_tx.maxLum; qdD_lumens_qrx1_txL = qrx1_gain_txL.*D_sh_qrx1_txL.*vlcCfg_tx.maxLum;
qdA_lumens_qrx2_txL = qrx2_gain_txL.*A_sh_qrx2_txL.*vlcCfg_tx.maxLum; qdB_lumens_qrx2_txL = qrx2_gain_txL.*B_sh_qrx2_txL.*vlcCfg_tx.maxLum;
qdC_lumens_qrx2_txL = qrx2_gain_txL.*C_sh_qrx2_txL.*vlcCfg_tx.maxLum; qdD_lumens_qrx2_txL = qrx2_gain_txL.*D_sh_qrx2_txL.*vlcCfg_tx.maxLum;

% Convert lumens to pkpk, we're assuming a 50ohm output impedance here,
% check the paper
qdA_pkpk_qrx1_txR = qdA_lumens_qrx1_txR*vlcCfg_qrx.lum2pkpk; qdB_pkpk_qrx1_txR = qdB_lumens_qrx1_txR*vlcCfg_qrx.lum2pkpk;
qdC_pkpk_qrx1_txR = qdC_lumens_qrx1_txR*vlcCfg_qrx.lum2pkpk; qdD_pkpk_qrx1_txR = qdD_lumens_qrx1_txR*vlcCfg_qrx.lum2pkpk;
qdA_pkpk_qrx2_txR = qdA_lumens_qrx2_txR*vlcCfg_qrx.lum2pkpk; qdB_pkpk_qrx2_txR = qdB_lumens_qrx2_txR*vlcCfg_qrx.lum2pkpk;
qdC_pkpk_qrx2_txR = qdC_lumens_qrx2_txR*vlcCfg_qrx.lum2pkpk; qdD_pkpk_qrx2_txR = qdD_lumens_qrx2_txR*vlcCfg_qrx.lum2pkpk;
qdA_pkpk_qrx1_txL = qdA_lumens_qrx1_txL*vlcCfg_qrx.lum2pkpk; qdB_pkpk_qrx1_txL = qdB_lumens_qrx1_txL*vlcCfg_qrx.lum2pkpk;
qdC_pkpk_qrx1_txL = qdC_lumens_qrx1_txL*vlcCfg_qrx.lum2pkpk; qdD_pkpk_qrx1_txL = qdD_lumens_qrx1_txL*vlcCfg_qrx.lum2pkpk;
qdA_pkpk_qrx2_txL = qdA_lumens_qrx2_txL*vlcCfg_qrx.lum2pkpk; qdB_pkpk_qrx2_txL = qdB_lumens_qrx2_txL*vlcCfg_qrx.lum2pkpk;
qdC_pkpk_qrx2_txL = qdC_lumens_qrx2_txL*vlcCfg_qrx.lum2pkpk; qdD_pkpk_qrx2_txL = qdD_lumens_qrx2_txL*vlcCfg_qrx.lum2pkpk;

figure
plot(vehCfg.t,qdA_pkpk_qrx1_txR), hold on
plot(vehCfg.t,qdB_pkpk_qrx1_txR,'r'), hold on
plot(vehCfg.t,qdC_pkpk_qrx1_txR,'g'), hold on
plot(vehCfg.t,qdD_pkpk_qrx1_txR,'c'), hold off

figure
plot(vehCfg.t,qdA_pkpk_qrx2_txR), hold on
plot(vehCfg.t,qdB_pkpk_qrx2_txR,'r'), hold on
plot(vehCfg.t,qdC_pkpk_qrx2_txR,'g'), hold on
plot(vehCfg.t,qdD_pkpk_qrx2_txR,'c'), hold off

% Pack for algorithm simulations
v2lcRun.qrxRight.qdA_pkpk.txR = qdA_pkpk_qrx1_txR; v2lcRun.qrxRight.qdB_pkpk.txR = qdB_pkpk_qrx1_txR;
v2lcRun.qrxRight.qdC_pkpk.txR = qdC_pkpk_qrx1_txR; v2lcRun.qrxRight.qdD_pkpk.txR = qdD_pkpk_qrx1_txR;
v2lcRun.qrxLeft.qdA_pkpk.txR = qdA_pkpk_qrx2_txR; v2lcRun.qrxLeft.qdB_pkpk.txR = qdB_pkpk_qrx2_txR;
v2lcRun.qrxLeft.qdC_pkpk.txR = qdC_pkpk_qrx2_txR; v2lcRun.qrxLeft.qdD_pkpk.txR = qdD_pkpk_qrx2_txR;
v2lcRun.qrxRight.qdA_pkpk.txL = qdA_pkpk_qrx1_txL; v2lcRun.qrxRight.qdB_pkpk.txL = qdB_pkpk_qrx1_txL;
v2lcRun.qrxRight.qdC_pkpk.txL = qdC_pkpk_qrx1_txL; v2lcRun.qrxRight.qdD_pkpk.txL = qdD_pkpk_qrx1_txL;
v2lcRun.qrxLeft.qdA_pkpk.txL = qdA_pkpk_qrx2_txL; v2lcRun.qrxLeft.qdB_pkpk.txL = qdB_pkpk_qrx2_txL;
v2lcRun.qrxLeft.qdC_pkpk.txL = qdC_pkpk_qrx2_txL; v2lcRun.qrxLeft.qdD_pkpk.txL = qdD_pkpk_qrx2_txL;

v2lcRun.qrxRight.gainPt.txR = qrx1_gain_txR_pt;
v2lcRun.qrxLeft.gainPt.txR = qrx2_gain_txR_pt; 
v2lcRun.qrxRight.gainPt.txL = qrx1_gain_txL_pt;
v2lcRun.qrxLeft.gainPt.txL = qrx2_gain_txL_pt;

answer = inputdlg('Enter filename for the simulation output file','Simulation Output Filename',[1 50],{'v2lcRun_<explanation>.mat'});
save(strcat('data/',answer{1}),'vehCfg','vlcCfg_qrx','vlcCfg_tx','v2lcRun');

