function [ out, dRes, thetaSwpRes ] = vlcCfgTool_rxThetaLutFit( handles )
%VLCCFGTOOL_THETALUTFIT Summary of this function goes here
%   Detailed explanation goes here

%%% Recalculate the look-up table (fQRX)

prompt = {'Enter angle resolution (deg):      Recommended: 1e-2 (3 mins on i5 laptop)','Enter metric resolution for QRX sim (m)   Recommended: 1e-6 (3 mins on i5 laptop):'};
dlgtitle = 'Parameter Selection for Theta Look-up Generation';
dims = [1 45];
definput = {'1e-2','1e-6'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
handles.thetaSwpRes = str2num(answer{1});
handles.dRes = str2num(answer{2});
thetaSwpRes = handles.thetaSwpRes;
dRes = handles.dRes;
[lut,~,~] = vlcCfgTool_fovLut(handles);
a1 = lut(1,:);
a2 = lut(2,:);
lut = [a2;a1];
f = waitbar(0, 'Fitting Curve for QRX Theta LUT');
out = fit(lut(1,:)',lut(2,:)','smoothingspline');
waitbar(100, f, 'Fitting Curve for QRX Theta LUT');
close(f)
end

