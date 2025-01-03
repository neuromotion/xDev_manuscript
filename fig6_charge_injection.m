%% fig6_charge_injection.m
% Author: Sam Parker
% Date: 8/13/2024
% (C) Brown Neuromotion Lab 2024
%
% PROVIDED WITHOUT WARRANTY OR GUARANTEE
clear; close all; clc
load("data\Fig6_charge_injection.mat")

damage_threshold_pC_per_cm2 = 30e-6 * 1e12; % Convert 30 micro-C per cm2 to pico-C (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5386002/#S9title)
dbs_electrode_area_cm2 = 0.06; % same ref as above
micro_electrode_area_cm2 = 220 * 1e-12 * 1e4; % Convert 220 um2 to cm2 (https://www.sciencedirect.com/science/article/pii/0006899394912130)
ees_electrode_area_cm2 = 3.4 * 1e-2; % Convert 3.4mm2 to cm2. Boston Scientific. 2020. Surgical Leads Directions for Use [User manual]. https://www.bostonscientific.com/content/dam/elabeling/nm/pr/92089830-02_Rev_A_Surgical_Leads_Directions_for_Use_DFU_multi_OUS_s.pdf

%% Remove amplifier gain
charge_injection.waveform.Voltage = charge_injection.waveform.Voltage ./ charge_injection.amplifier_gain;

%% Zero Crossing Detection
sign_signal = sign(charge_injection.waveform.Voltage - mean(charge_injection.waveform.Voltage));
zero_crossing_idx = find(abs(diff(sign_signal)))+1;

%% Sample waveform at high and low state
num_samples_per_level = mean(diff(zero_crossing_idx));
level_check_idxs = round(zero_crossing_idx + num_samples_per_level/2);

high_level_idxs = level_check_idxs(1:2:end-1);
low_level_idxs = level_check_idxs(2:2:end-1);

delta_v = charge_injection.waveform.Voltage(high_level_idxs) - charge_injection.waveform.Voltage(low_level_idxs);

%% Average capacitor measurements
Cap = mean(charge_injection.cap_data);

%% Calculate current injection using Q = C * V
all_cinj = Cap .* delta_v;

%% Plot
figure()
tcl = tiledlayout(1, 2);
ax(1) = nexttile();
histogram(all_cinj*1e12)
xlabel("Q_{inj} (pC)")
ylabel("Number of Transitions")
title("Charge Injection")
xlim([8.5, 9.5])
ax(2) = nexttile();
histogram(all_cinj*1e12); hold on
plot([damage_threshold_pC_per_cm2, damage_threshold_pC_per_cm2] .*  dbs_electrode_area_cm2, [0, 20], "--r");
plot([damage_threshold_pC_per_cm2, damage_threshold_pC_per_cm2] .*  micro_electrode_area_cm2, [0, 20], "--g");
plot([damage_threshold_pC_per_cm2, damage_threshold_pC_per_cm2] .*  ees_electrode_area_cm2, [0, 20], "--b");
xlabel("Q_{inj} (pC)")
ylabel("Number of Transitions")
set(gca, "XScale", "log")
title("Distribution vs Damage Thresholds")
subtitle("Areas: MEA 200\mu{}m^2, EES 3.4mm^2, DBS 0.06cm^2")
legend("Charge Injected", "DBS Damage Threshold", "MEA Damage Threshold", "EES Damage Threshold", "Location", "North")
title(tcl, "Charge Injection")
set(gcf, "Position", [100, 100, 1310, 420])

%% Plot artifact on ephys channel

figure()
yyaxis right;
plot(charge_injection.artifact.Time_s * 1e6, charge_injection.artifact.PCLK_V, '-', "Color", "#000000");
ylabel("Digital Voltage (V)")
yticks(0:2.5:5)
yyaxis left;
plot(charge_injection.artifact.Time_s * 1e6, charge_injection.artifact.X4_V * 1e6, "Color", "#33bbee");
ylabel("Artifact voltage (μV)")
ylim([-400, 400])
yticks(-400:200:400)
xlabel("Time (μs)")
title("Charge injection artifact into 1k electrode impedance")