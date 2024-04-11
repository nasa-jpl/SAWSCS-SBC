% Copyright 2024 by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------

function [trial] = read_trial(mp, trial_name)

ifPlot = 0;

switch trial_name
    case 'trial4'
        load([mp.falcodir 'macos/hex_maps/trial4'])
        trial = trial3;
    case 'trial6'
        load([mp.falcodir 'macos/hex_maps/trial6'])
    case 'trial7'
	load([mp.falcodir 'macos/hex_maps/trial7'])
    case 'trial8'
	load([mp.falcodir 'macos/hex_maps/trial8'])
	trial = trial(1:8);
    case 'trial9'
	load([mp.falcodir 'macos/hex_maps/trial9'])
    case 'trial10'
	load([mp.falcodir 'macos/hex_maps/trial10'])
    case 'trial11'
        load([mp.falcodir 'macos/hex_maps/trial11'])
    case 'trial12'
        load([mp.falcodir 'macos/hex_maps/trial12'])
    case 'trial13'
        load([mp.falcodir 'macos/hex_maps/trial13'])

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ifPlot
  for jj = 1:size(trial,2)
    figure(4),show_opd(1e3*trial(jj).opd0,'opd0','nm')    
    figure(5),show_opd(1e3*trial(jj).opd_wfsc,'opd__wfsc','nm')
    figure(6),show_opd(1e3*trial(jj).drift(1).opd_drift,'drift.opd__drift','nm')
    figure(7),show_opd(1e3*trial(jj).drift(1).opd_final,'drift.opd__final','nm')
    figure(8),show_opd(1e3*zernike_remove(trial(jj).delta(1).dopd_drift,[1:3]),'delta.dopd__drift','nm')
    figure(9),show_opd(1e6*zernike_remove(trial(jj).delta(1).dopd_final,[1:3]),'delta.dopd__final','pm')
    drawnow
  end
end
