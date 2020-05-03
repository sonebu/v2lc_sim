function [ qrx1Lum, qrx2Lum ] = vlcCfgTool_rxGainLum2PkPk( handles )
%vlcCfgTool_rxGainLum2PkPk Summary of this function goes here
%   Detailed explanation goes here

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% WE'VE MADE AN APPROXIMATION HERE                                   %%%
%%%                                                                    %%%
%%% At the angular resolution we want (max 0.0001 deg), MATLAB can't   %%%
%%% store the 3d surface representing the radiation pattern. So volume %%%
%%% under that curve was approximated by Erdem's simit method.         %%%
%%% Compute AUC for xy, revolve it 360 deg. repeat for zy              %%%
%%% You get two rings, extract the part from xy which spans angle by zy%%%
%%% repeat for zy -> xy get the average                                %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% UPDATE: There's also a rectangle approximation actually...              %%%
%%% Since our source is radially symmetric, the tilt and azimuth angles are %%%
%%% affixed to the base of the 3d plot, this defines a rectangle. We divide %%%
%%% this into as many parts as the resolution we desire (maybe we can think %%%
%%% of efficiency gains with downsampling in the future). The values these  %%%
%%% parts take in the perpendicular axis of the 3d plot (vertical) lets     %%%
%%% call them phi and ksi, sqrt(phi^2 + ksi^2) is our original rad_pat val  %%%
%%% THIS ONLY WORKS BECAUSE THE SOURCE IS RADIALLY SYMMETRIC!!              %%%
%%% Non-radially symmetric sources can maybe be represented as sum of       %%%
%%% radially symmetric sources but that should be automatized and currently %%%
%%% is not implemented. This would actually be good because automotive beam %%%
%%% patterns are pretty irregular                                           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pd_ddl = handles.dDL;
fitto = handles.txPlrPatFit;
z = 0;
beta = -deg2rad(handles.hdg); % ?

% left qrx 
x = handles.lat + 1;
y = handles.fwd;
% Geometric calc on XY
psi1_xy = atan2(y,x + pd_ddl/2);
psi2_xy = atan2(y,x - pd_ddl/2);
psi3_xy = psi2_xy - psi1_xy;
% eps1_xy = rad2deg((deg2rad(90)-psi2_xy)-atan2(dir_vect_xy(1),dir_vect_xy(2))); % this axis is not used, flat roads
eps1_xy = rad2deg((deg2rad(90)-psi2_xy)-beta);
eps2_xy = rad2deg(deg2rad(eps1_xy) + psi3_xy);

% Geometric calc on ZY
psi1_zy = atan2(y,z + pd_ddl/2);
psi2_zy = atan2(y,z - pd_ddl/2);
psi3_zy = psi2_zy - psi1_zy;
% eps1_zy = rad2deg((deg2rad(90)-psi2_zy)-atan2(dir_vect_zy(1),dir_vect_zy(2))); % this axis is not used, flat roads
eps1_zy = rad2deg((deg2rad(90)-psi2_zy)-beta);
eps2_zy = rad2deg(deg2rad(eps1_zy) + psi3_zy);

% Our own integrator, obsolete, but keep it for completeness
% vol = vlp_meas_rxGain_radSymSrc3dIntegral( rad_pat, eps1_xy, eps2_xy, eps1_zy, eps2_zy );

% New version, just use quad2d, MATLAB's integrator
vol = quad2d(fitto,eps1_xy,eps2_xy,eps1_zy,eps2_zy);
qrx1Lum = vol*handles.maxLum;


% right qrx 
x = handles.lat - 1;
y = handles.fwd;
% Geometric calc on XY
psi1_xy = atan2(y,x + pd_ddl/2);
psi2_xy = atan2(y,x - pd_ddl/2);
psi3_xy = psi2_xy - psi1_xy;
% eps1_xy = rad2deg((deg2rad(90)-psi2_xy)-atan2(dir_vect_xy(1),dir_vect_xy(2))); % this axis is not used, flat roads
eps1_xy = rad2deg((deg2rad(90)-psi2_xy)-beta);
eps2_xy = rad2deg(deg2rad(eps1_xy) + psi3_xy);

% Geometric calc on ZY
psi1_zy = atan2(y,z + pd_ddl/2);
psi2_zy = atan2(y,z - pd_ddl/2);
psi3_zy = psi2_zy - psi1_zy;
% eps1_zy = rad2deg((deg2rad(90)-psi2_zy)-atan2(dir_vect_zy(1),dir_vect_zy(2))); % this axis is not used, flat roads
eps1_zy = rad2deg((deg2rad(90)-psi2_zy)-beta);
eps2_zy = rad2deg(deg2rad(eps1_zy) + psi3_zy);

vol = quad2d(fitto,eps1_xy,eps2_xy,eps1_zy,eps2_zy);
qrx2Lum = vol*handles.maxLum;

%%% Previously we used our own integrator but this isn't actually necessary
%%% so dropping it for now... The radPat tool uses it, check that one if
%%% you need to remember
% vol = vlp_meas_rxGain_radSymSrc3dIntegral( rad_pat, eps1_xy, eps2_xy, eps1_zy, eps2_zy );

end

