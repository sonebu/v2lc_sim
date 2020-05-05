clear, clc, close all;
addpath('fcn');

% Load the polar pattern for your LED.
% NOTE: check the example for requirements. Change the values on that one
% to create new ones. The pattern has to be radially symmetric (about 0 deg) 
% and needs to span +/- 90*sqrt(2) since we use a square-symmetric trapezoidal 
% integrator (see vlcCfgTxPlrPatTool_radSymSrc3dIntegral.m). This didn't
% work as expected with quad2d, that's why we had to write our own integrator
[fn, pt] = uigetfile('data/*.mat','Select numeric polar pattern array, named "vlcCfgTxPlrPatNum_<>.mat", from file');
load(strcat(pt,fn));

% Compute volume under 3D rotated polar pattern, will be used for
% normalization later
vol = vlcCfgTxPlrPatTool_radSymSrc3dIntegral( vlcCfgTxPlrPatNum, -90, 90, -90, 90);

% Subsample pattern for faster surface fitting
answer = inputdlg('Surface fitting takes long, so if you want to subsample to make it faster (and the fit less resolute!) enter an integer >1','TxPlrPat Subsample',[1 50],{'50'});
subsmpl = str2num(answer{1});
a1 = vlcCfgTxPlrPatNum(1,1:subsmpl:end);
a2 = vlcCfgTxPlrPatNum(2,1:subsmpl:end);
vlcCfgTxPlrPatNum_dwn = [a1;a2];

% Define the initial profile
x = vlcCfgTxPlrPatNum_dwn(1,:);
a = abs(x-(-90));
m90_pos = find(a == min(a));
a = abs(x-(90));
p90_pos = find(a == min(a));
y = vlcCfgTxPlrPatNum_dwn(2,:);

% Compute normalized 3D pattern (it's a 2D matrix with "height" values for the surface)
pat = zeros(length(y));
for i = m90_pos:p90_pos
    for j = m90_pos:p90_pos
        d = abs(x-sqrt(x(i)*x(i)+x(j)*x(j)));
        a = find(d == min(d));
        pat(i,j) = y(a)/vol;
    end
end

% Surface fitting
clearvars -except pat x y vlcCfgTxPlrPatNum_dwn vlcCfgTxPlrPatNum
x_long = repmat(x,[1 length(x)]);
y_long = zeros(size(x_long));
for i=1:length(x_long)
    y_long(i) = x_long(ceil(i/sqrt(length(x_long))));
end
pat_long = reshape(pat,[1 size(pat,1)*size(pat,2)]);
vlcCfgTxPlrPatFit = fit([x_long',y_long'],pat_long','cubicinterp');

% Test the fit, tic toc to understand how much it's going to take to
% calculate a 2 deg (hrz and vert) rectangular area during the actual
% simulation. May want to increase the subsampling rate to make it faster
plot(vlcCfgTxPlrPatFit,[x_long',y_long'],pat_long')
quad2d(vlcCfgTxPlrPatFit,-1,1,-1,1)

% Save the fit
answer = inputdlg('Enter filename for the surface fit + polar pattern file','TxPlrPat Filename',[1 50],{'vlcCfgTxPlrPat_<explanation>.mat'});
save(strcat('data/',answer{1}),'vlcCfgTxPlrPatNum','vlcCfgTxPlrPatFit');
clearvars -except vlcCfgTxPlrPatFit vlcCfgTxPlrPatNum vlcCfgTxPlrPatNum_dwn
