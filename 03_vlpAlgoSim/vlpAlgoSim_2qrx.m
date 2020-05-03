clear all; clc;
% close all;
addpath('fcn');

%% Load simulation output file (the algorithm only uses QRX measurements and the QRX configuration)
[fn, pt] = uigetfile('../02_v2lcDataGen/data/*.mat','Load Simulation Output');
load(strcat(pt,fn));
clear vlcCfg_tx

%% Algorithm Parameters
vlpAlgoSimCfg.adcSmplRate = 1e5; % samples/sec
vlpAlgoSimCfg.adcBitDepth = 16; % bits
vlpAlgoSimCfg.thetaLutFit = vlcCfg_qrx.thetaLut;
vlpAlgoSimCfg.vlpRate = 50; % Hz

% For faster simulation and less storage requirement "vlp_sim.m" only
% computes the lumens reaching QRX, the signal is constructed here when
% it's needed, with the addition of noise. The SNR specified here is used
% while constructing the signal. This is OK since the chn model is AWGN.
vlpAlgoSimCfg.ns_dBm = -50; % -40: 10 mVpkpk , -60: 1 mVpkpk 

%% Simulation post-processing

% Simulation: save actual teta values for reference
thetaAct_qrx1_txL = atan2(vehCfg.tgt.leftTailPosY,vehCfg.tgt.leftTailPosX + (vehCfg.ego.vehWidth)/2);
thetaAct_qrx2_txL = atan2(vehCfg.tgt.leftTailPosY,vehCfg.tgt.leftTailPosX - (vehCfg.ego.vehWidth)/2);
thetaAct_qrx1_txR = atan2(vehCfg.tgt.rightTailPosY,vehCfg.tgt.rightTailPosX + (vehCfg.ego.vehWidth)/2);
thetaAct_qrx2_txR = atan2(vehCfg.tgt.rightTailPosY,vehCfg.tgt.rightTailPosX - (vehCfg.ego.vehWidth)/2);

% Simulation: interpolate the pkpk measurements for each qrx to signal
interp_method = 'pchip';
vlcCfg_qrx.t_dt = 1/vlpAlgoSimCfg.adcSmplRate;
vlcCfg_qrx.t = [0:vlcCfg_qrx.t_dt:(vehCfg.t(end))-vlcCfg_qrx.t_dt]; % seconds
tmp = linspace(1,length(v2lcRun.qrxRight.qdA_pkpk.txR),length(vlcCfg_qrx.t));
qdA_pkpk_qrx1_txR_sig = interp1(1:length(v2lcRun.qrxRight.qdA_pkpk.txR),v2lcRun.qrxRight.qdA_pkpk.txR,tmp,interp_method);
qdB_pkpk_qrx1_txR_sig = interp1(1:length(v2lcRun.qrxRight.qdB_pkpk.txR),v2lcRun.qrxRight.qdB_pkpk.txR,tmp,interp_method);
qdC_pkpk_qrx1_txR_sig = interp1(1:length(v2lcRun.qrxRight.qdC_pkpk.txR),v2lcRun.qrxRight.qdC_pkpk.txR,tmp,interp_method);
qdD_pkpk_qrx1_txR_sig = interp1(1:length(v2lcRun.qrxRight.qdD_pkpk.txR),v2lcRun.qrxRight.qdD_pkpk.txR,tmp,interp_method);
qdA_pkpk_qrx2_txR_sig = interp1(1:length(v2lcRun.qrxLeft.qdA_pkpk.txR), v2lcRun.qrxLeft.qdA_pkpk.txR,tmp,interp_method);
qdB_pkpk_qrx2_txR_sig = interp1(1:length(v2lcRun.qrxLeft.qdB_pkpk.txR), v2lcRun.qrxLeft.qdB_pkpk.txR,tmp,interp_method);
qdC_pkpk_qrx2_txR_sig = interp1(1:length(v2lcRun.qrxLeft.qdC_pkpk.txR), v2lcRun.qrxLeft.qdC_pkpk.txR,tmp,interp_method);
qdD_pkpk_qrx2_txR_sig = interp1(1:length(v2lcRun.qrxLeft.qdD_pkpk.txR), v2lcRun.qrxLeft.qdD_pkpk.txR,tmp,interp_method);
qdA_pkpk_qrx1_txL_sig = interp1(1:length(v2lcRun.qrxRight.qdA_pkpk.txL),v2lcRun.qrxRight.qdA_pkpk.txL,tmp,interp_method);
qdB_pkpk_qrx1_txL_sig = interp1(1:length(v2lcRun.qrxRight.qdB_pkpk.txL),v2lcRun.qrxRight.qdB_pkpk.txL,tmp,interp_method);
qdC_pkpk_qrx1_txL_sig = interp1(1:length(v2lcRun.qrxRight.qdC_pkpk.txL),v2lcRun.qrxRight.qdC_pkpk.txL,tmp,interp_method);
qdD_pkpk_qrx1_txL_sig = interp1(1:length(v2lcRun.qrxRight.qdD_pkpk.txL),v2lcRun.qrxRight.qdD_pkpk.txL,tmp,interp_method);
qdA_pkpk_qrx2_txL_sig = interp1(1:length(v2lcRun.qrxLeft.qdA_pkpk.txL), v2lcRun.qrxLeft.qdA_pkpk.txL,tmp,interp_method);
qdB_pkpk_qrx2_txL_sig = interp1(1:length(v2lcRun.qrxLeft.qdB_pkpk.txL), v2lcRun.qrxLeft.qdB_pkpk.txL,tmp,interp_method);
qdC_pkpk_qrx2_txL_sig = interp1(1:length(v2lcRun.qrxLeft.qdC_pkpk.txL), v2lcRun.qrxLeft.qdC_pkpk.txL,tmp,interp_method);
qdD_pkpk_qrx2_txL_sig = interp1(1:length(v2lcRun.qrxLeft.qdD_pkpk.txL), v2lcRun.qrxLeft.qdD_pkpk.txL,tmp,interp_method);
qrx1.qdA_eps_txR = qdA_pkpk_qrx1_txR_sig; qrx1.qdB_eps_txR = qdB_pkpk_qrx1_txR_sig; 
qrx1.qdC_eps_txR = qdC_pkpk_qrx1_txR_sig; qrx1.qdD_eps_txR = qdD_pkpk_qrx1_txR_sig;
qrx2.qdA_eps_txR = qdA_pkpk_qrx2_txR_sig; qrx2.qdB_eps_txR = qdB_pkpk_qrx2_txR_sig; 
qrx2.qdC_eps_txR = qdC_pkpk_qrx2_txR_sig; qrx2.qdD_eps_txR = qdD_pkpk_qrx2_txR_sig;
qrx1.qdA_eps_txL = qdA_pkpk_qrx1_txL_sig; qrx1.qdB_eps_txL = qdB_pkpk_qrx1_txL_sig; 
qrx1.qdC_eps_txL = qdC_pkpk_qrx1_txL_sig; qrx1.qdD_eps_txL = qdD_pkpk_qrx1_txL_sig;
qrx2.qdA_eps_txL = qdA_pkpk_qrx2_txL_sig; qrx2.qdB_eps_txL = qdB_pkpk_qrx2_txL_sig; 
qrx2.qdC_eps_txL = qdC_pkpk_qrx2_txL_sig; qrx2.qdD_eps_txL = qdD_pkpk_qrx2_txL_sig;

% Simulation: Pass the interpolated pkpk measurements through an adc
tmp = round((2^vlpAlgoSimCfg.adcBitDepth)*(qdA_pkpk_qrx1_txR_sig));
qdA_pkpk_qrx1_txR_sig = 5*(tmp)/(2^vlpAlgoSimCfg.adcBitDepth);
tmp = round((2^vlpAlgoSimCfg.adcBitDepth)*(qdB_pkpk_qrx1_txR_sig));
qdB_pkpk_qrx1_txR_sig = 5*(tmp)/(2^vlpAlgoSimCfg.adcBitDepth);
tmp = round((2^vlpAlgoSimCfg.adcBitDepth)*(qdC_pkpk_qrx1_txR_sig));
qdC_pkpk_qrx1_txR_sig = 5*(tmp)/(2^vlpAlgoSimCfg.adcBitDepth);
tmp = round((2^vlpAlgoSimCfg.adcBitDepth)*(qdD_pkpk_qrx1_txR_sig));
qdD_pkpk_qrx1_txR_sig = 5*(tmp)/(2^vlpAlgoSimCfg.adcBitDepth);
tmp = round((2^vlpAlgoSimCfg.adcBitDepth)*(qdA_pkpk_qrx2_txR_sig));
qdA_pkpk_qrx2_txR_sig = 5*(tmp)/(2^vlpAlgoSimCfg.adcBitDepth);
tmp = round((2^vlpAlgoSimCfg.adcBitDepth)*(qdB_pkpk_qrx2_txR_sig));
qdB_pkpk_qrx2_txR_sig = 5*(tmp)/(2^vlpAlgoSimCfg.adcBitDepth);
tmp = round((2^vlpAlgoSimCfg.adcBitDepth)*(qdC_pkpk_qrx2_txR_sig));
qdC_pkpk_qrx2_txR_sig = 5*(tmp)/(2^vlpAlgoSimCfg.adcBitDepth);
tmp = round((2^vlpAlgoSimCfg.adcBitDepth)*(qdD_pkpk_qrx2_txR_sig));
qdD_pkpk_qrx2_txR_sig = 5*(tmp)/(2^vlpAlgoSimCfg.adcBitDepth);

tmp = round((2^vlpAlgoSimCfg.adcBitDepth)*(qdA_pkpk_qrx1_txL_sig));
qdA_pkpk_qrx1_txL_sig = 5*(tmp)/(2^vlpAlgoSimCfg.adcBitDepth);
tmp = round((2^vlpAlgoSimCfg.adcBitDepth)*(qdB_pkpk_qrx1_txL_sig));
qdB_pkpk_qrx1_txL_sig = 5*(tmp)/(2^vlpAlgoSimCfg.adcBitDepth);
tmp = round((2^vlpAlgoSimCfg.adcBitDepth)*(qdC_pkpk_qrx1_txL_sig));
qdC_pkpk_qrx1_txL_sig = 5*(tmp)/(2^vlpAlgoSimCfg.adcBitDepth);
tmp = round((2^vlpAlgoSimCfg.adcBitDepth)*(qdD_pkpk_qrx1_txL_sig));
qdD_pkpk_qrx1_txL_sig = 5*(tmp)/(2^vlpAlgoSimCfg.adcBitDepth);
tmp = round((2^vlpAlgoSimCfg.adcBitDepth)*(qdA_pkpk_qrx2_txL_sig));
qdA_pkpk_qrx2_txL_sig = 5*(tmp)/(2^vlpAlgoSimCfg.adcBitDepth);
tmp = round((2^vlpAlgoSimCfg.adcBitDepth)*(qdB_pkpk_qrx2_txL_sig));
qdB_pkpk_qrx2_txL_sig = 5*(tmp)/(2^vlpAlgoSimCfg.adcBitDepth);
tmp = round((2^vlpAlgoSimCfg.adcBitDepth)*(qdC_pkpk_qrx2_txL_sig));
qdC_pkpk_qrx2_txL_sig = 5*(tmp)/(2^vlpAlgoSimCfg.adcBitDepth);
tmp = round((2^vlpAlgoSimCfg.adcBitDepth)*(qdD_pkpk_qrx2_txL_sig));
qdD_pkpk_qrx2_txL_sig = 5*(tmp)/(2^vlpAlgoSimCfg.adcBitDepth);
clear tmp

% Simulation: Generate TX waveform and use qrx qd gains to
%             obtain the final measurements
rng(1);
txR_sym = 2500; % bps
txR_wavFreqs = [5000 7500]; % Hz
txR_wav = zeros(length(vlcCfg_qrx.t),1);
txR_ttlNumSyms = txR_sym*vlcCfg_qrx.t(end);
txR_symNumSmpl = 1/(txR_sym*(vlcCfg_qrx.t_dt));
for i=1:txR_ttlNumSyms
    txR_wav((1+(i-1)*txR_symNumSmpl):(i*txR_symNumSmpl)) = ...
        sin(2*pi*txR_wavFreqs(1+(rand>0.5))*vlcCfg_qrx.t((1+(i-1)*txR_symNumSmpl):(i*txR_symNumSmpl)));
end

txL_sym = 2500; % bps
txL_wavFreqs = [10000 12500]; % Hz
txL_wav = zeros(length(vlcCfg_qrx.t),1);
txL_ttlNumSyms = txL_sym*vlcCfg_qrx.t(end);
txL_symNumSmpl = 1/(txL_sym*(vlcCfg_qrx.t_dt));
for i=1:txL_ttlNumSyms
    txL_wav((1+(i-1)*txL_symNumSmpl):(i*txL_symNumSmpl)) = ...
        sin(2*pi*txL_wavFreqs(1+(rand>0.5))*vlcCfg_qrx.t((1+(i-1)*txL_symNumSmpl):(i*txL_symNumSmpl)));
end

qrx1_qdA_txR = qdA_pkpk_qrx1_txR_sig'.*(txR_wav/2);
qrx1_qdB_txR = qdB_pkpk_qrx1_txR_sig'.*(txR_wav/2); 
qrx1_qdC_txR = qdC_pkpk_qrx1_txR_sig'.*(txR_wav/2); 
qrx1_qdD_txR = qdD_pkpk_qrx1_txR_sig'.*(txR_wav/2); 

qrx2_qdA_txR = qdA_pkpk_qrx2_txR_sig'.*(txR_wav/2); 
qrx2_qdB_txR = qdB_pkpk_qrx2_txR_sig'.*(txR_wav/2); 
qrx2_qdC_txR = qdC_pkpk_qrx2_txR_sig'.*(txR_wav/2); 
qrx2_qdD_txR = qdD_pkpk_qrx2_txR_sig'.*(txR_wav/2); 

qrx1_qdA_txL = qdA_pkpk_qrx1_txL_sig'.*(txL_wav/2); 
qrx1_qdB_txL = qdB_pkpk_qrx1_txL_sig'.*(txL_wav/2); 
qrx1_qdC_txL = qdC_pkpk_qrx1_txL_sig'.*(txL_wav/2); 
qrx1_qdD_txL = qdD_pkpk_qrx1_txL_sig'.*(txL_wav/2); 

qrx2_qdA_txL = qdA_pkpk_qrx2_txL_sig'.*(txL_wav/2); 
qrx2_qdB_txL = qdB_pkpk_qrx2_txL_sig'.*(txL_wav/2); 
qrx2_qdC_txL = qdC_pkpk_qrx2_txL_sig'.*(txL_wav/2); 
qrx2_qdD_txL = qdD_pkpk_qrx2_txL_sig'.*(txL_wav/2); 

clear qd*_pkpk_qrx*_tx*_sig

% AC intensities add up, no dc bias (yet)
qrx1_qdA = (qrx1_qdA_txR + qrx1_qdA_txL)/2; qrx1_qdB = (qrx1_qdB_txR + qrx1_qdB_txL)/2;
qrx1_qdC = (qrx1_qdC_txR + qrx1_qdC_txL)/2; qrx1_qdD = (qrx1_qdD_txR + qrx1_qdD_txL)/2;
qrx2_qdA = (qrx2_qdA_txR + qrx2_qdA_txL)/2; qrx2_qdB = (qrx2_qdB_txR + qrx2_qdB_txL)/2;
qrx2_qdC = (qrx2_qdC_txR + qrx2_qdC_txL)/2; qrx2_qdD = (qrx2_qdD_txR + qrx2_qdD_txL)/2;

clear qrx*_qd*_tx*

% Noise power computation
ns_dBW = vlpAlgoSimCfg.ns_dBm - 30;
ns_W = 10^(ns_dBW/10);
ns_vRMS = sqrt(ns_W * 50); % 50 ohm system
ns_n = randn(length(txR_wav),1);
ns_mult = ns_vRMS/std(ns_n);

% For AVERAGE SNR computation, just for debugging
qrx_pwr = (var(qrx1_qdA)+var(qrx2_qdA)+var(qrx1_qdB)+var(qrx2_qdB)+var(qrx1_qdC)+var(qrx2_qdC)+var(qrx1_qdD)+var(qrx2_qdD))/8;
snr = 10*log10(qrx_pwr/var(ns_mult*(randn(length(txR_wav),1))));

% Add DC bias and noise
qrx1_qdA = 2.5 + qrx1_qdA + ns_mult*(randn(length(txR_wav),1)); qrx1_qdA(qrx1_qdA>5) = 5; qrx1_qdA(qrx1_qdA<0) = 0;
qrx1_qdB = 2.5 + qrx1_qdB + ns_mult*(randn(length(txR_wav),1)); qrx1_qdB(qrx1_qdB>5) = 5; qrx1_qdB(qrx1_qdB<0) = 0;
qrx1_qdC = 2.5 + qrx1_qdC + ns_mult*(randn(length(txR_wav),1)); qrx1_qdC(qrx1_qdC>5) = 5; qrx1_qdC(qrx1_qdC<0) = 0;
qrx1_qdD = 2.5 + qrx1_qdD + ns_mult*(randn(length(txR_wav),1)); qrx1_qdD(qrx1_qdD>5) = 5; qrx1_qdD(qrx1_qdD<0) = 0;

qrx2_qdA = 2.5 + qrx2_qdA + ns_mult*(randn(length(txR_wav),1)); qrx2_qdA(qrx2_qdA>5) = 5; qrx2_qdA(qrx2_qdA<0) = 0;
qrx2_qdB = 2.5 + qrx2_qdB + ns_mult*(randn(length(txR_wav),1)); qrx2_qdB(qrx2_qdB>5) = 5; qrx2_qdB(qrx2_qdB<0) = 0;
qrx2_qdC = 2.5 + qrx2_qdC + ns_mult*(randn(length(txR_wav),1)); qrx2_qdC(qrx2_qdC>5) = 5; qrx2_qdC(qrx2_qdC<0) = 0;
qrx2_qdD = 2.5 + qrx2_qdD + ns_mult*(randn(length(txR_wav),1)); qrx2_qdD(qrx2_qdD>5) = 5; qrx2_qdD(qrx2_qdD<0) = 0;

% Pack qrx data
qrx1.qdA = qrx1_qdA; qrx1.qdB = qrx1_qdB; qrx1.qdC = qrx1_qdC; qrx1.qdD = qrx1_qdD; 
qrx2.qdA = qrx2_qdA; qrx2.qdB = qrx2_qdB; qrx2.qdC = qrx2_qdC; qrx2.qdD = qrx2_qdD; 
qrx1.thetaAct_txL = thetaAct_qrx1_txL; qrx1.thetaAct_txR = thetaAct_qrx1_txR;
qrx2.thetaAct_txL = thetaAct_qrx2_txL; qrx2.thetaAct_txR = thetaAct_qrx2_txR;

clear fn i interp_method pt qd*_pkpk_qrx* qrx*_qd*_* tx_maxLum tx_plrPatFit tx_plrPatNum tx_sym tx_symNumSmpl tx_ttlNumSyms tx_wavFreqs

%% Algorithm

fii = waitbar(0, 'Calculating');
vlpDecimRate = round(1/(vlpAlgoSimCfg.vlpRate*vehCfg.t_dt)); % algorithm gives output every *this* many cycles of movement
numTetaSmplPerPos = int32(round(vlpAlgoSimCfg.adcSmplRate*vehCfg.t_dt)*vlpDecimRate);
vlpLen = int32(round(length(vehCfg.tgt.lat_rel)/vlpDecimRate));
thetaCalc_qrx1_txR = zeros(vlpLen,1); thetaCalc_qrx2_txR = zeros(vlpLen,1);
thetaCalc_qrx1_txL = zeros(vlpLen,1); thetaCalc_qrx2_txL = zeros(vlpLen,1);
calc_lat_txR_rel = zeros(vlpLen,1); calc_fwd_txR_rel = zeros(vlpLen,1);
calc_lat_txL_rel = zeros(vlpLen,1); calc_fwd_txL_rel = zeros(vlpLen,1);
g_QRX =  vlpAlgoSimCfg.thetaLutFit;
for i=1:vlpLen
    waitbar(double(i)/double(vlpLen),fii, 'Calculating');

    % Algo: Find theta, QRX 1
    qdA = qrx1_qdA(1+(i-1)*numTetaSmplPerPos:i*numTetaSmplPerPos)-2.5;
    qdB = qrx1_qdB(1+(i-1)*numTetaSmplPerPos:i*numTetaSmplPerPos)-2.5;
    qdC = qrx1_qdC(1+(i-1)*numTetaSmplPerPos:i*numTetaSmplPerPos)-2.5;
    qdD = qrx1_qdD(1+(i-1)*numTetaSmplPerPos:i*numTetaSmplPerPos)-2.5;
    [~, thetaCalc_qrx1_txR(i), ~, a] = vlpAlgoSim_getTeta( g_QRX, qdA, qdB, qdC, qdD, txR_wav(1+(i-1)*numTetaSmplPerPos:i*numTetaSmplPerPos), 2e-2);
    [~, thetaCalc_qrx1_txL(i), ~, b] = vlpAlgoSim_getTeta( g_QRX, qdA, qdB, qdC, qdD, txL_wav(1+(i-1)*numTetaSmplPerPos:i*numTetaSmplPerPos), 2e-2);
    vldFlag_qrx1 = a*b;
    
    % Algo: Find theta, QRX 2
    qdA = qrx2_qdA(1+(i-1)*numTetaSmplPerPos:i*numTetaSmplPerPos)-2.5;
    qdB = qrx2_qdB(1+(i-1)*numTetaSmplPerPos:i*numTetaSmplPerPos)-2.5;
    qdC = qrx2_qdC(1+(i-1)*numTetaSmplPerPos:i*numTetaSmplPerPos)-2.5;
    qdD = qrx2_qdD(1+(i-1)*numTetaSmplPerPos:i*numTetaSmplPerPos)-2.5;
    [~, thetaCalc_qrx2_txR(i), ~, a] = vlpAlgoSim_getTeta( g_QRX, qdA, qdB, qdC, qdD, txR_wav(1+(i-1)*numTetaSmplPerPos:i*numTetaSmplPerPos), 2e-2);
    [~, thetaCalc_qrx2_txL(i), ~, b] = vlpAlgoSim_getTeta( g_QRX, qdA, qdB, qdC, qdD, txL_wav(1+(i-1)*numTetaSmplPerPos:i*numTetaSmplPerPos), 2e-2);
    vldFlag_qrx2 = a*b;
    
    % When the target is out of the FoV, obviously it's an invalid result.
    % we detect this in the function and assign a NaN to those results here
    % so that they don't appear on the plot. It's just extra clutter when
    % they're plotted.
	if( (vldFlag_qrx1 == 0) || (vldFlag_qrx2 == 0))
		calc_lat_txR_rel(i) = NaN;
		calc_fwd_txR_rel(i) = NaN;
		calc_lat_txL_rel(i) = NaN;
		calc_fwd_txL_rel(i) = NaN;
	else		
			
		% Algo: find position, use sine law!
		a_prio_d = vehCfg.ego.vehWidth;
		calc_lat_txR_rel(i) = a_prio_d*(1/2 + sind(thetaCalc_qrx2_txR(i))*cosd(thetaCalc_qrx1_txR(i))/sind(thetaCalc_qrx1_txR(i)-thetaCalc_qrx2_txR(i)));
		calc_fwd_txR_rel(i) = a_prio_d*cosd(thetaCalc_qrx2_txR(i))*cosd(thetaCalc_qrx1_txR(i))/sind(thetaCalc_qrx1_txR(i)-thetaCalc_qrx2_txR(i));
		calc_lat_txL_rel(i) = a_prio_d*(1/2 + sind(thetaCalc_qrx2_txL(i))*cosd(thetaCalc_qrx1_txL(i))/sind(thetaCalc_qrx1_txL(i)-thetaCalc_qrx2_txL(i)));
		calc_fwd_txL_rel(i) = a_prio_d*cosd(thetaCalc_qrx2_txL(i))*cosd(thetaCalc_qrx1_txL(i))/sind(thetaCalc_qrx1_txL(i)-thetaCalc_qrx2_txL(i));
    end

end
close(fii)

%%% Post-process (show midpoint for clarity)
actl_bump_center_x = (vehCfg.tgt.rightTailPosX(1:int32(vlpDecimRate):end) + vehCfg.tgt.leftTailPosX(1:int32(vlpDecimRate):end))/2;
actl_bump_center_y = (vehCfg.tgt.rightTailPosY(1:int32(vlpDecimRate):end) + vehCfg.tgt.leftTailPosY(1:int32(vlpDecimRate):end))/2;
calc_bump_center_x = (calc_lat_txR_rel + calc_lat_txL_rel)/2;
calc_bump_center_y = (calc_fwd_txR_rel + calc_fwd_txL_rel)/2;
calc_bump_center_x(end) = calc_bump_center_x(end-1); % workaround for NaN issue
calc_bump_center_y(end) = calc_bump_center_y(end-1); % workaround for NaN issue
e_x = actl_bump_center_x-calc_bump_center_x';
e_y = actl_bump_center_y-calc_bump_center_y';
[ var_x, var_y ] = vlpAlgoSim_2qrx_calcVarBound( ns_mult, vehCfg, vlcCfg_qrx, txL_wav, txR_wav, qrx1, qrx2, vlpDecimRate);

% Pack for results script
vlpAlgoRun.DecimRate = vlpDecimRate;
vlpAlgoRun.rx_ns_pkpk = ns_mult;
vlpAlgoRun.vehCfg = vehCfg;
vlpAlgoRun.vlcCfg_qrx = vlcCfg_qrx;
vlpAlgoRun.v2lcRun= v2lcRun;
vlpAlgoRun.actl_x = actl_bump_center_x;
vlpAlgoRun.actl_y = actl_bump_center_y;
vlpAlgoRun.calc_x = calc_bump_center_x';
vlpAlgoRun.calc_y = calc_bump_center_y';
vlpAlgoRun.bound_var_x = var_x;
vlpAlgoRun.bound_var_y = var_y;

%%%% XY Grid Plot of Actual versus Esimated
figure
ca = plot(calc_bump_center_x,calc_bump_center_y, 'r','LineWidth',0.5);
grid on, hold on,xlabel('Ego x [m]'),ylabel('Ego y [m]'); 
ac = plot(actl_bump_center_x,actl_bump_center_y, 'b','LineWidth',1);
legend([ca ac],'Location','Northwest','Estimated', 'True');
xlim_curr = get(gca,'xlim'); 
mx_x = max(abs(xlim_curr(1)),abs(xlim_curr(2)));
xlim_curr(1) = -(mx_x+0.1*mx_x);
xlim_curr(2) = +(mx_x+0.1*mx_x);
set(gca,'xlim',xlim_curr)
ylim_curr = get(gca,'ylim'); 
mx_y = max(abs(ylim_curr(1)),abs(ylim_curr(2)));
ylim_curr(1) = -1;
ylim_curr(2) = +(mx_y+0.1*mx_y);
set(gca,'ylim',ylim_curr)
hold off; 

answer = questdlg('Save simulation output?', 'Save?', 'Yes','No','No');
% Handle response
switch answer
    case 'Yes'
        answer = inputdlg('Enter filename for the output file','Output Filename',[1 50],{'vlpAlgoRun_<explanation>.mat'});
        save(strcat('data/',answer{1}),'vlpAlgoRun');
    case 'No'
        close all;
end


