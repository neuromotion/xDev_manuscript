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
r_vals_k = [1, 10, 100];

figure("Renderer","painters");
tcl = tiledlayout(1, 3, "TileSpacing", "tight");
for i = 1:length(r_vals_k)
    nexttile()
    yyaxis left;
    plotMeanAndStandard(charge_injection.artifact.(sprintf("R%dk", r_vals_k(i))).Time_us', charge_injection.artifact.(sprintf("R%dk", r_vals_k(i))).PCLK_V', "#262626", "#262626");
    ylabel("PCLK (V)")
    yyaxis right
    plotMeanAndStandard(charge_injection.artifact.(sprintf("R%dk", r_vals_k(i))).Time_us', charge_injection.artifact.(sprintf("R%dk", r_vals_k(i))).CH_4_uV', "#33bbee", "#33bbee");
    ylabel("Charge Injection Artifact (uV)")
    xlabel("Time (us)")
    title(sprintf("Electrode Z = %dR", r_vals_k(i)*1000))
end
set(gcf, "Position", [63, 1, 1858, 1123]);

figure("Renderer","painters");
colors = ["#06d6a0", "#b04cb9", "#33bbee"];
for i = 1:length(r_vals_k)
    charge_injection.artifact.(sprintf("R%dk", r_vals_k(i))).Current_uA = charge_injection.artifact.(sprintf("R%dk", r_vals_k(i))).CH_4_uV ./ (r_vals_k(i) * 1000);
    plotMeanAndStandard(charge_injection.artifact.(sprintf("R%dk", r_vals_k(i))).Time_us', charge_injection.artifact.(sprintf("R%dk", r_vals_k(i))).Current_uA' .* 1000, colors(i), colors(i)); hold on
end
plot([0, 0], [-10, 20], "--", "Color", "#888888");
title("All")
ylabel("Injected Current (nA)")
xlabel("Time (μs)")
yticks(-10:10:20);
xticks(0:60:120);
ylim([-10, 20]);
xlim([min(charge_injection.artifact.(sprintf("R%dk", r_vals_k(i))).Time_us, [], 'all'), max(charge_injection.artifact.(sprintf("R%dk", r_vals_k(i))).Time_us, [], 'all')]);
legend("R = 1k", "R = 10k", "R = 100k", "Connection made", "Location", "northeast")

%% Functions

function plotMeanAndStandard(time, data_mat, lineColor, fillColor)
    data_mean = nanmean(data_mat);
    data_std = nanstd(data_mat);
    time = nanmean(time);

    mask = ~isnan(time) & ~isnan(data_mean) & ~isnan(data_std);
    time = time(mask);
    data_mean = data_mean(mask);
    data_std = data_std(mask);

    curve1 = data_mean + data_std;
    curve2 = data_mean - data_std;
    x2 = [time, fliplr(time)];
    inBetween = [curve1,  fliplr(curve2)];
    fill(x2, inBetween, 'g', "FaceColor", fillColor, "LineStyle", "none", "FaceAlpha", 0.5, "HandleVisibility", "off"); hold on
    plot(time, data_mean, '-', 'Color', lineColor, 'LineWidth', 2);
end