function [ var_x, var_y ] = vlpAlgoSim_2qrx_calcVarBound( rx_ns_PkPk, vehCfg, vlcCfg_qrx, txL_wav, txR_wav, qrx1, qrx2, vlpDecimRate)
%VLPALGOSIM_CALCVARBOUND Summary of this function goes here
%   Detailed explanation goes here

%%% See the original submission of the manuscript for these computations,
%%% nothing new. 

var_q = var(rx_ns_PkPk*randn(1,length(qrx1.qdA)));
tbuf = int32(vehCfg.t_dt/vlcCfg_qrx.t_dt)*vlpDecimRate;

for ii = 1:length(txL_wav)/tbuf

    %%%%%%%%% For TX L %%%%%%%%%
    var_E = 0;
    for i=1:tbuf
        delE_delQ = txL_wav(i+(ii-1)*tbuf)/double(tbuf);
        var_E = var_E + ((delE_delQ)^2)*var_q;
    end
    var_E_arr(ii) = var_E;

    E_A_qrx1 = mean(qrx1.qdA_eps_txL( (1+(ii-1)*tbuf):ii*tbuf));
    E_B_qrx1 = mean(qrx1.qdB_eps_txL( (1+(ii-1)*tbuf):ii*tbuf));
    E_C_qrx1 = mean(qrx1.qdC_eps_txL( (1+(ii-1)*tbuf):ii*tbuf));
    E_D_qrx1 = mean(qrx1.qdD_eps_txL( (1+(ii-1)*tbuf):ii*tbuf));
    Phi_qrx1 = ((E_B_qrx1 + E_D_qrx1)-(E_A_qrx1 + E_C_qrx1))/(E_A_qrx1 + E_B_qrx1 + E_C_qrx1 + E_D_qrx1);
    delFLut_delPhi_qrx1 = deg2rad((feval(vlcCfg_qrx.thetaLut,Phi_qrx1+0.00001) - feval(vlcCfg_qrx.thetaLut,Phi_qrx1-0.00001))/(0.00002));
    % for A/C
    delPhi_delE_A_qrx1 = -2*(E_B_qrx1 + E_D_qrx1)/((E_A_qrx1+E_B_qrx1+E_C_qrx1+E_D_qrx1)^2);
    delPhi_delE_C_qrx1 = -2*(E_B_qrx1 + E_D_qrx1)/((E_A_qrx1+E_B_qrx1+E_C_qrx1+E_D_qrx1)^2);
    % for B/D
    delPhi_delE_B_qrx1 = 2*(E_A_qrx1 + E_C_qrx1)/((E_A_qrx1+E_B_qrx1+E_C_qrx1+E_D_qrx1)^2);
    delPhi_delE_D_qrx1 = 2*(E_A_qrx1 + E_C_qrx1)/((E_A_qrx1+E_B_qrx1+E_C_qrx1+E_D_qrx1)^2);
    
    var_Theta_qrx1 = ((delFLut_delPhi_qrx1*delPhi_delE_A_qrx1)^2)*var_E + ...
                     ((delFLut_delPhi_qrx1*delPhi_delE_B_qrx1)^2)*var_E + ...
                     ((delFLut_delPhi_qrx1*delPhi_delE_C_qrx1)^2)*var_E + ...
                     ((delFLut_delPhi_qrx1*delPhi_delE_D_qrx1)^2)*var_E;
    
    var_Theta_qrx1_arr(ii) = var_Theta_qrx1;
    
    E_A_qrx2 = mean(qrx2.qdA_eps_txL( (1+(ii-1)*tbuf):ii*tbuf));
    E_B_qrx2 = mean(qrx2.qdB_eps_txL( (1+(ii-1)*tbuf):ii*tbuf));
    E_C_qrx2 = mean(qrx2.qdC_eps_txL( (1+(ii-1)*tbuf):ii*tbuf));
    E_D_qrx2 = mean(qrx2.qdD_eps_txL( (1+(ii-1)*tbuf):ii*tbuf));
    Phi_qrx2 = ((E_B_qrx2 + E_D_qrx2)-(E_A_qrx2 + E_C_qrx2))/(E_A_qrx2 + E_B_qrx2 + E_C_qrx2 + E_D_qrx2);
    delFLut_delPhi_qrx2 = deg2rad((feval(vlcCfg_qrx.thetaLut,Phi_qrx2+0.00001) - feval(vlcCfg_qrx.thetaLut,Phi_qrx2-0.00001))/(0.00002));
    % for A/C
    delPhi_delE_A_qrx2 = -2*(E_B_qrx2 + E_D_qrx2)/((E_A_qrx2+E_B_qrx2+E_C_qrx2+E_D_qrx2)^2);
    delPhi_delE_C_qrx2 = -2*(E_B_qrx2 + E_D_qrx2)/((E_A_qrx2+E_B_qrx2+E_C_qrx2+E_D_qrx2)^2);
    % for B/D
    delPhi_delE_B_qrx2 = 2*(E_A_qrx2 + E_C_qrx2)/((E_A_qrx2+E_B_qrx2+E_C_qrx2+E_D_qrx2)^2);
    delPhi_delE_D_qrx2 = 2*(E_A_qrx2 + E_C_qrx2)/((E_A_qrx2+E_B_qrx2+E_C_qrx2+E_D_qrx2)^2);
    
    var_Theta_qrx2 = ((delFLut_delPhi_qrx2*delPhi_delE_A_qrx2)^2)*var_E + ...
                     ((delFLut_delPhi_qrx2*delPhi_delE_B_qrx2)^2)*var_E + ...
                     ((delFLut_delPhi_qrx2*delPhi_delE_C_qrx2)^2)*var_E + ...
                     ((delFLut_delPhi_qrx2*delPhi_delE_D_qrx2)^2)*var_E;

    var_Theta_qrx2_arr(ii) = var_Theta_qrx2;
    
    tR = qrx1.thetaAct_txL;
    tL = qrx2.thetaAct_txL;
    delX_delThetaQrx1 = vehCfg.ego.vehWidth*( sin(tL(ii))*cos(tL(ii))/((sin(tL(ii))*cos(tR(ii))-cos(tL(ii))*sin(tR(ii)))^2));
    delX_delThetaQrx2 = vehCfg.ego.vehWidth*(-sin(tR(ii))*cos(tR(ii))/((sin(tL(ii))*cos(tR(ii))-cos(tL(ii))*sin(tR(ii)))^2));
    delY_delThetaQrx1 = vehCfg.ego.vehWidth*(( cos(tL(ii))^2)/((sin(tL(ii))*cos(tR(ii))-cos(tL(ii))*sin(tR(ii)))^2));
    delY_delThetaQrx2 = vehCfg.ego.vehWidth*((-cos(tR(ii))^2)/((sin(tL(ii))*cos(tR(ii))-cos(tL(ii))*sin(tR(ii)))^2));
    
    var_x_txL(ii) = (delX_delThetaQrx1^2)*var_Theta_qrx1 + (delX_delThetaQrx2^2)*var_Theta_qrx2;
    var_y_txL(ii) = (delY_delThetaQrx1^2)*var_Theta_qrx1 + (delY_delThetaQrx2^2)*var_Theta_qrx2;
    
	%%%%%%%%% For TX R %%%%%%%%%
    var_E = 0;
    for i=1:tbuf
        delE_delQ = txR_wav(i+(ii-1)*tbuf)/double(tbuf);
        var_E = var_E + ((delE_delQ)^2)*var_q;
    end
    var_E_arr(ii) = var_E;
    
    E_A_qrx1 = mean(qrx1.qdA_eps_txR( (1+(ii-1)*tbuf):ii*tbuf));
    E_B_qrx1 = mean(qrx1.qdB_eps_txR( (1+(ii-1)*tbuf):ii*tbuf));
    E_C_qrx1 = mean(qrx1.qdC_eps_txR( (1+(ii-1)*tbuf):ii*tbuf));
    E_D_qrx1 = mean(qrx1.qdD_eps_txR( (1+(ii-1)*tbuf):ii*tbuf));
    Phi_qrx1 = ((E_B_qrx1 + E_D_qrx1)-(E_A_qrx1 + E_C_qrx1))/(E_A_qrx1 + E_B_qrx1 + E_C_qrx1 + E_D_qrx1);
    delFLut_delPhi_qrx1 = deg2rad((feval(vlcCfg_qrx.thetaLut,Phi_qrx1+0.00001) - feval(vlcCfg_qrx.thetaLut,Phi_qrx1-0.00001))/(0.00002));     % just an extra divide by zero guard
    % for A/C
    delPhi_delE_A_qrx1 = -2*(E_B_qrx1 + E_D_qrx1)/((E_A_qrx1+E_B_qrx1+E_C_qrx1+E_D_qrx1)^2);
    delPhi_delE_C_qrx1 = -2*(E_B_qrx1 + E_D_qrx1)/((E_A_qrx1+E_B_qrx1+E_C_qrx1+E_D_qrx1)^2);
    % for B/D
    delPhi_delE_B_qrx1 = 2*(E_A_qrx1 + E_C_qrx1)/((E_A_qrx1+E_B_qrx1+E_C_qrx1+E_D_qrx1)^2);
    delPhi_delE_D_qrx1 = 2*(E_A_qrx1 + E_C_qrx1)/((E_A_qrx1+E_B_qrx1+E_C_qrx1+E_D_qrx1)^2);
    
    var_Theta_qrx1 = ((delFLut_delPhi_qrx1*delPhi_delE_A_qrx1)^2)*var_E + ...
                     ((delFLut_delPhi_qrx1*delPhi_delE_B_qrx1)^2)*var_E + ...
                     ((delFLut_delPhi_qrx1*delPhi_delE_C_qrx1)^2)*var_E + ...
                     ((delFLut_delPhi_qrx1*delPhi_delE_D_qrx1)^2)*var_E;
    
    var_Theta_qrx1_arr(ii) = var_Theta_qrx1;
    
    E_A_qrx2 = mean(qrx2.qdA_eps_txR( (1+(ii-1)*tbuf):ii*tbuf));
    E_B_qrx2 = mean(qrx2.qdB_eps_txR( (1+(ii-1)*tbuf):ii*tbuf));
    E_C_qrx2 = mean(qrx2.qdC_eps_txR( (1+(ii-1)*tbuf):ii*tbuf));
    E_D_qrx2 = mean(qrx2.qdD_eps_txR( (1+(ii-1)*tbuf):ii*tbuf));
    Phi_qrx2 = ((E_B_qrx2 + E_D_qrx2)-(E_A_qrx2 + E_C_qrx2))/(E_A_qrx2 + E_B_qrx2 + E_C_qrx2 + E_D_qrx2);
    delFLut_delPhi_qrx2 = deg2rad((feval(vlcCfg_qrx.thetaLut,Phi_qrx2+0.00001) - feval(vlcCfg_qrx.thetaLut,Phi_qrx2-0.00001))/(0.00002));     % just an extra divide by zero guard
    % for A/C
    delPhi_delE_A_qrx2 = -2*(E_B_qrx2 + E_D_qrx2)/((E_A_qrx2+E_B_qrx2+E_C_qrx2+E_D_qrx2)^2);
    delPhi_delE_C_qrx2 = -2*(E_B_qrx2 + E_D_qrx2)/((E_A_qrx2+E_B_qrx2+E_C_qrx2+E_D_qrx2)^2);
    % for B/D
    delPhi_delE_B_qrx2 = 2*(E_A_qrx2 + E_C_qrx2)/((E_A_qrx2+E_B_qrx2+E_C_qrx2+E_D_qrx2)^2);
    delPhi_delE_D_qrx2 = 2*(E_A_qrx2 + E_C_qrx2)/((E_A_qrx2+E_B_qrx2+E_C_qrx2+E_D_qrx2)^2);
    
    var_Theta_qrx2 = ((delFLut_delPhi_qrx2*delPhi_delE_A_qrx2)^2)*var_E + ...
                     ((delFLut_delPhi_qrx2*delPhi_delE_B_qrx2)^2)*var_E + ...
                     ((delFLut_delPhi_qrx2*delPhi_delE_C_qrx2)^2)*var_E + ...
                     ((delFLut_delPhi_qrx2*delPhi_delE_D_qrx2)^2)*var_E;

    var_Theta_qrx2_arr(ii) = var_Theta_qrx2;
    
    tR = qrx1.thetaAct_txR;
    tL = qrx2.thetaAct_txR;
    delX_delThetaQrx1 = vehCfg.ego.vehWidth*( sin(tL(ii))*cos(tL(ii))/((sin(tL(ii))*cos(tR(ii))-cos(tL(ii))*sin(tR(ii)))^2));
    delX_delThetaQrx2 = vehCfg.ego.vehWidth*(-sin(tR(ii))*cos(tR(ii))/((sin(tL(ii))*cos(tR(ii))-cos(tL(ii))*sin(tR(ii)))^2));
    delY_delThetaQrx1 = vehCfg.ego.vehWidth*(( cos(tL(ii))^2)/((sin(tL(ii))*cos(tR(ii))-cos(tL(ii))*sin(tR(ii)))^2));
    delY_delThetaQrx2 = vehCfg.ego.vehWidth*((-cos(tR(ii))^2)/((sin(tL(ii))*cos(tR(ii))-cos(tL(ii))*sin(tR(ii)))^2));
    
    var_x_txR(ii) = (delX_delThetaQrx1^2)*var_Theta_qrx1 + (delX_delThetaQrx2^2)*var_Theta_qrx2;
    var_y_txR(ii) = (delY_delThetaQrx1^2)*var_Theta_qrx1 + (delY_delThetaQrx2^2)*var_Theta_qrx2;
   
end

    var_x = var_x_txR + var_x_txL;
    var_y = var_y_txR + var_y_txL;

end

