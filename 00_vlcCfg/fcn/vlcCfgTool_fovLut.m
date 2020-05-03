function [ lut, oob_pos_neg, oob_pos_pst ] = vlcCfgTool_fovLut( handles )
%VLCCFGTOOL_FOVLUT Summary of this function goes here
%   Detailed explanation goes here

%%% Plot fQRX/gQRX, Ray optics simulations basically..

dDL = handles.dDL;
dSL = handles.dSL;
dFL = handles.dFL;
dPX = handles.dPX;
dr = (dDL/2)*(dFL-dSL)/dFL;
dRes = handles.dRes;

theta = deg2rad(-90:handles.thetaSwpRes:90);
y_est_arr = zeros(size(theta));
y_arr = zeros(size(theta));

f = waitbar(0, 'Running FoV calculation');
for ii = 1:length(theta)
    waitbar(ii/length(theta),f, 'Running FoV calculation');
    qpd_y = tan(theta(ii))*dSL;
    qpd_vertical = 0; % Future proofing for 3D
    y_arr(ii) = qpd_y;
    
    % Generate bounds for looping
    srch_x_neg = qpd_vertical-dr;
    srch_x_pos = qpd_vertical+dr;
    x_ar = srch_x_neg:dRes:srch_x_pos;
    y_ar_pos = sqrt(abs(dr.^2-(x_ar-qpd_vertical).^2))+qpd_y;
    y_ar_neg = -sqrt(abs(dr.^2-(x_ar-qpd_vertical).^2))+qpd_y;
    
    x_filt = (x_ar <= dPX) & (x_ar >= -dPX);
    x_vals=x_ar(x_filt);
    
    y_p_vals = y_ar_pos(x_filt);
    y_p_vals(y_p_vals>dPX) = dPX;
    y_p_vals(y_p_vals<-dPX) = -dPX;
    y_n_vals = y_ar_neg(x_filt);
    y_n_vals(y_n_vals>dPX) = dPX;
    y_n_vals(y_n_vals<-dPX) = -dPX;
    
    A_ctr = 0; B_ctr = 0; C_ctr = 0; D_ctr = 0;
    for i=1:length(x_vals)
        y_vals = y_n_vals(i):dRes:y_p_vals(i);
        % Works for uniform only! But for the spot size we use, uniform is
        % a good assumption, see Manojlovic's QPD sensitivity paper for
        % more information.
        if(x_vals(i) > 0)
            s1 = y_vals > 0;
            a = sum(s1);
            B_ctr = B_ctr + a;
            D_ctr = D_ctr + length(s1) - a;
        else
            s1 = y_vals > 0;
            a = sum(s1);
            A_ctr = A_ctr + a;
            C_ctr = C_ctr + length(s1) - a;
        end
        
    end
    
    y_est = ((A_ctr+B_ctr)-(C_ctr+D_ctr))/(A_ctr+B_ctr+C_ctr+D_ctr);
    y_est_arr(ii) = y_est;
    
end
close(f)

cla(handles.thetaLutAxis)
theta_arr = -90:handles.thetaSwpRes:90;
plot(handles.thetaLutAxis,theta_arr, y_est_arr,'g','LineWidth', 1.5)
grid(handles.thetaLutAxis,'on')
hold(handles.thetaLutAxis,'on')
oob = find(abs(y_est_arr) == 1)-length(y_est_arr)/2;
oob_pos_neg = theta_arr(length(y_est_arr)/2-min(abs(oob)));
oob_pos_pst = theta_arr(length(y_est_arr)/2+min(abs(oob))+1);
plot(handles.thetaLutAxis,[oob_pos_neg oob_pos_neg], [-2 2],'c','LineWidth',1)
plot(handles.thetaLutAxis,[oob_pos_pst oob_pos_pst], [-2 2],'c','LineWidth',1)
plot(handles.thetaLutAxis,[oob_pos_neg oob_pos_pst], [-1 1],'k--','LineWidth',1)
lut = [theta_arr(length(y_est_arr)/2-min(abs(oob)) : length(y_est_arr)/2+min(abs(oob))+1);
       y_est_arr(length(y_est_arr)/2-min(abs(oob)) : length(y_est_arr)/2+min(abs(oob))+1)];
xlim(handles.thetaLutAxis,[oob_pos_neg-1 oob_pos_pst+1]);
ylim(handles.thetaLutAxis,[-1 1]);
xlabel(handles.thetaLutAxis,'Theta (deg)')
ylabel(handles.thetaLutAxis,'Spot Position (p.u.)')
hold(handles.thetaLutAxis,'off')

%%% Did have this zero point detection earlier but it's actually clutter for
%%% no real useful reason so keep it out for now:
% zero_det = find(isnan(y_est_arr))-length(y_est_arr)/2;
% zero_det_pos_neg = teta_arr(length(y_est_arr)/2-min(abs(zero_det)));
% zero_det_pos_pst = teta_arr(length(y_est_arr)/2+min(abs(zero_det))+1);
% plot([zero_det_pos_neg zero_det_pos_neg], [-2 2],'r','LineWidth',1)
% plot([zero_det_pos_pst zero_det_pos_pst], [-2 2],'r','LineWidth',1)


end

