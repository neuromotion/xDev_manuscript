%% fig12_invivo.m
% Author: Sam Parker
% Date: 9/30/2024
% (C) Brown Neuromotion Lab 2024
%
% PROVIDED WITHOUT WARRANTY OR GUARANTEE
clear; close all; clc
load("data\Fig12_invivo.mat")

%% Plot the EMG examples from the Intan chip
figure()
tcl = tiledlayout(2,1);
nexttile()
plot(intan.example_50Hz.emg_data.t, intan.example_50Hz.emg_data{:, 2:end});
title(sprintf("Intan Ch%d Amp: %d\\mu{}A Freq: %dHz", intan.example_50Hz.stim_data.Chan, intan.example_50Hz.stim_data.Amplitude, intan.example_50Hz.stim_data.Freq))
nexttile()
plot(intan.example_25Hz.emg_data.t, intan.example_25Hz.emg_data{:, 2:end});
title(sprintf("Intan Ch%d Amp: %d\\mu{}A Freq: %dHz", intan.example_25Hz.stim_data.Chan, intan.example_25Hz.stim_data.Amplitude, intan.example_25Hz.stim_data.Freq))
muscle_names = string([intan.example_25Hz.emg_data.Properties.VariableNames(2:end)]);
legend(muscle_names, "Location", "southoutside", "Orientation","horizontal");

%% Plot the ECAP examples from the Ripple NIP
figure()
for i = 1:size(ripple.artifact_times, 1)
    plot(ripple.artifact_times(i, :), ripple.artifact_values(i, :), '-', "Color", [187/255, 187/255, 187/255, 0.1]); hold on
    plot(ripple.ecap_times(i, :), ripple.ecap_values(i, :), '-', "Color", [51/255, 188/255, 238/255, 0.2]); hold on
end
xlim([0, 2])
ylim([-14000, -5000])
xlabel("Time (ms)")
ylabel("ECAP (\mu{}V)")
title("Stim E5 Record E58+ E53-")
