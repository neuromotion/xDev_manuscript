%% fig10_biosignals.m
% Author: Sam Parker
% Date: 9/12/2024
% (C) Brown Neuromotion Lab 2024
%
% Note: NeuroDAC data is in mV
%
% PROVIDED WITHOUT WARRANTY OR GUARANTEE
clear; close all; clc

load("data\Fig10_biosignals.mat")

%% Plot example traces

figure()
tcl = tiledlayout(2, 1);
nexttile()
plot(sines.example_trace.time, sines.example_trace.neuroDAC, '-k'); hold on
plot(sines.example_trace.time, sines.example_trace.xDev, "Color", "#33bbee");
ylabel("Sine Example Trace (mV)")

nexttile()
plot(biosignals.example_trace.nhp_cortical.time, biosignals.example_trace.nhp_cortical.neuroDAC, '-k'); hold on
plot(biosignals.example_trace.nhp_cortical.time, biosignals.example_trace.nhp_cortical.xDev, "Color", "#33bbee");
ylabel("NHP Cortical Example Trace (mV)")
xlabel("Time (s)")
legend("NeuroDAC", "xDev", "Location", "Southoutside", "Orientation", "horizontal")


%% Compute correlations across frequencies
clearvars -except biosignals sines

figure()
semilogx(sines.frequencies, mean(sines.xDev.correlations, 2), "Color", "#06D6A0"); hold on
scatter(sines.frequencies, sines.xDev.correlations, "MarkerFaceColor", "#06D6A0", "MarkerEdgeColor", "k", "LineWidth",0.1, "HandleVisibility","off");
semilogx(sines.frequencies, mean(sines.baseline.correlations, 2), "Color", "#B04CB9");
scatter(sines.frequencies, sines.baseline.correlations, "MarkerFaceColor", "#B04CB9", "MarkerEdgeColor", "k", "LineWidth",0.1, "HandleVisibility","off");
ylabel("Correlation coefficient")
title("Input-Output correlations")
xlabel("Sine Frequency")
ylim([0, 1])
legend("xDev", "Ceiling", "Location", "southeast");

%% Compute correlations across biosignals
clearvars -except biosignals sines

figure()
data = [mean(biosignals.baseline.correlations); mean(biosignals.xDev.correlations)];
b = bar(transpose(data)); hold on
b(1).FaceColor = "#B04CB9";
b(2).FaceColor = "#06D6A0";
scatter(0.85:1:2.85, transpose(biosignals.baseline.correlations), "MarkerFaceColor", "#B04CB9", "MarkerEdgeColor", "k", "LineWidth",0.1, "HandleVisibility","off");
scatter(1.15:1:3.15, transpose(biosignals.xDev.correlations), "MarkerFaceColor", "#06D6A0", "MarkerEdgeColor", "k", "LineWidth",0.1, "HandleVisibility","off");
ylim([0, 1])
xticks(1:3)
xticklabels(strrep(biosignals.order, "_", " "))
title("Correlation between biosignals")
ylabel("Correlation Coefficient")
legend("Baseline", "xDev", "Location", "southoutside", "Orientation", "horizontal")

%% Compute power-in-band

figure()
tcl = tiledlayout(2, 1, "TileSpacing", "tight");
nexttile();
% Plot Sheep LFP
baseline_data = mean(biosignals.baseline.sheep_lfp.power_in_band.bands ./ biosignals.baseline.sheep_lfp.power_in_band.total, 2);
xdev_data = mean(biosignals.xDev.sheep_lfp.power_in_band.bands ./ biosignals.xDev.sheep_lfp.power_in_band.total, 2);
b = bar([baseline_data, xdev_data]); hold on
b(1).FaceColor = "#B04CB9";
b(2).FaceColor = "#06D6A0";
scatter(0.85:1:5.85, transpose(biosignals.baseline.sheep_lfp.power_in_band.bands ./ biosignals.baseline.sheep_lfp.power_in_band.total), "MarkerFaceColor", "#B04CB9", "MarkerEdgeColor", "k", "LineWidth",0.1, "HandleVisibility","off");
scatter(1.15:1:6.15, transpose(biosignals.xDev.sheep_lfp.power_in_band.bands ./ biosignals.xDev.sheep_lfp.power_in_band.total), "MarkerFaceColor", "#06D6A0", "MarkerEdgeColor", "k", "LineWidth",0.1, "HandleVisibility","off");
ylabel("Power in band (%)")
xlabel("Frequency Band");
xticks(1:6)
xticklabels(biosignals.baseline.sheep_lfp.power_in_band.band_names)
title("Sheep LFP")
legend("Baseline", "xDev")

baseline_data = (biosignals.baseline.sheep_lfp.power_in_band.bands ./ biosignals.baseline.sheep_lfp.power_in_band.total)';
xdev_data = (biosignals.xDev.sheep_lfp.power_in_band.bands ./ biosignals.xDev.sheep_lfp.power_in_band.total)';
p_vals = nan(size(xdev_data, 2), 1);
for i = 1:length(p_vals)
    p_vals(i) = ranksum(baseline_data(:, i), xdev_data(:, i));
end
if any(p_vals < 0.05)
    fprintf("Found statistically significant difference in %s band in Sheep LFP!\n", biosignals.baseline.sheep_emg.power_in_band.band_names(p_vals < 0.05))
else
    disp("No differences in band power for Sheep LFP")
end

% Plot Sheep EMG
nexttile();
baseline_data = mean(biosignals.baseline.sheep_emg.power_in_band.bands ./ biosignals.baseline.sheep_emg.power_in_band.total, 2);
xdev_data = mean(biosignals.xDev.sheep_emg.power_in_band.bands ./ biosignals.xDev.sheep_emg.power_in_band.total, 2);
b = bar([baseline_data, xdev_data]); hold on
b(1).FaceColor = "#B04CB9";
b(2).FaceColor = "#06D6A0";
scatter(0.85:1:5.85, transpose(biosignals.baseline.sheep_emg.power_in_band.bands ./ biosignals.baseline.sheep_emg.power_in_band.total), "MarkerFaceColor", "#B04CB9", "MarkerEdgeColor", "k", "LineWidth",0.1, "HandleVisibility","off");
scatter(1.15:1:6.15, transpose(biosignals.xDev.sheep_emg.power_in_band.bands ./ biosignals.xDev.sheep_emg.power_in_band.total), "MarkerFaceColor", "#06D6A0", "MarkerEdgeColor", "k", "LineWidth",0.1, "HandleVisibility","off");
ylabel("Power in band (%)")
xlabel("Frequency Band");
xticks(1:6)
xticklabels(biosignals.baseline.sheep_emg.power_in_band.band_names)
title("Sheep EMG")
legend("Baseline", "xDev")

baseline_data = (biosignals.baseline.sheep_emg.power_in_band.bands ./ biosignals.baseline.sheep_emg.power_in_band.total)';
xdev_data = (biosignals.xDev.sheep_emg.power_in_band.bands ./ biosignals.xDev.sheep_emg.power_in_band.total)';
p_vals = nan(size(xdev_data, 2), 1);
for i = 1:length(p_vals)
    p_vals(i) = ranksum(baseline_data(:, i), xdev_data(:, i));
end
if any(p_vals < 0.05)
    fprintf("Found statistically significant difference in %s band in Sheep EMG!\n", biosignals.baseline.sheep_emg.power_in_band.band_names(p_vals < 0.05))
else
    disp("No differences in band power for Sheep LFP")
end

%% Plot NHP cortical spike raster plot

colors = [174, 234, 219;
    6, 214, 160;
    204, 156, 210;
    176, 76, 185;
    ]./ 255;
setting = ["xDev", "baseline"];
state = ["output", "input"];
offsets = [4:-1:1];
offsets = sort([offsets, offsets+0.3, offsets-0.3], "descend");
colorIdx = 1;
offsetIdx = 1;

figure()
for settingIdx = 1:length(setting)
    for stateIdx = 1:length(state)
        for rep = 1:3
            if rep == 1
                scatter(biosignals.(setting(settingIdx)).nhp_cortical_spikes.(state(stateIdx)){rep}, offsets(offsetIdx)*ones(size(biosignals.(setting(settingIdx)).nhp_cortical_spikes.(state(stateIdx)){rep})), 'filled', "MarkerFaceColor", colors(colorIdx, :)); hold on
            else
                scatter(biosignals.(setting(settingIdx)).nhp_cortical_spikes.(state(stateIdx)){rep}, offsets(offsetIdx)*ones(size(biosignals.(setting(settingIdx)).nhp_cortical_spikes.(state(stateIdx)){rep})), 'filled', "MarkerFaceColor", colors(colorIdx, :), "HandleVisibility", "off"); hold on
            end
        offsetIdx = offsetIdx+1;
        end
        colorIdx = colorIdx+1;
    end
end
title("NHP cortical spikes")
xlabel("Time (s)")
legend("xDev Output", "xDev Input", "Baseline Output", "Baseline Input", "Location", "Southoutside", "orientation", "horizontal")

%% Use moving window to find distribution of missed spikes

start_time = 0;
window_length_s = 1;
overlap_pct = 66.6666667;

last_start_time = max([biosignals.xDev.nhp_cortical_spikes.output{:}, biosignals.xDev.nhp_cortical_spikes.input{:}, biosignals.baseline.nhp_cortical_spikes.output{:}, biosignals.baseline.nhp_cortical_spikes.input{:}]) - window_length_s;

input_output_error_xdev = [];
input_output_error_no_xdev = [];
start_times = [];
i = 1;
while(start_time < last_start_time)
    stop_time = start_time + window_length_s;
    start_times(i) = start_time;
    for rep = 1:3
        input_output_error_xdev(i, rep) = (sum(biosignals.xDev.nhp_cortical_spikes.input{rep} >= start_time & biosignals.xDev.nhp_cortical_spikes.input{rep} <= stop_time) - sum(biosignals.xDev.nhp_cortical_spikes.output{rep} >= start_time & biosignals.xDev.nhp_cortical_spikes.output{rep} <= stop_time));
        input_output_error_no_xdev(i, rep) = (sum(biosignals.baseline.nhp_cortical_spikes.input{rep} >= start_time & biosignals.baseline.nhp_cortical_spikes.input{rep} <= stop_time) - sum(biosignals.baseline.nhp_cortical_spikes.output{rep} >= start_time & biosignals.baseline.nhp_cortical_spikes.output{rep} <= stop_time));
    end
    start_time = start_time + overlap_pct / 100 * window_length_s;
    i = i+1;
end

figure()
boxplot(abs([reshape(input_output_error_no_xdev, [], 1), reshape(input_output_error_xdev, [], 1)])); hold on
scatter((rand(numel(input_output_error_no_xdev), 1)-0.5)*0.3 + 1, abs(reshape(input_output_error_no_xdev, [], 1)) + (rand(numel(input_output_error_no_xdev), 1)-0.5)*0.1, 15, "filled", "MarkerEdgeColor", "none");
scatter((rand(numel(input_output_error_xdev), 1)-0.5)*0.3 + 2, abs(reshape(input_output_error_xdev, [], 1)) + (rand(numel(input_output_error_xdev), 1)-0.5)*0.1, 15, "filled", "MarkerEdgeColor", "none");
title("Distribution of spike errors")
ylabel("Number of missed spikes")
xticklabels(["No xDev", "xDev"])
set(gcf, "Position", [345, 357, 333, 420])

p = kruskalwallis([reshape(abs(input_output_error_xdev), [], 1), reshape(abs(input_output_error_no_xdev), [], 1)], [], "off");
subtitle(sprintf("p = %0.4f | n = %d", p, numel(input_output_error_xdev)));

%% Plot average spike from neuroDAC and original data

load("data\original_spikes.mat")
load("data\neuroDAC_spikes.mat")
nd_mean = mean(nd_spike_waveforms); nd_std = std(nd_spike_waveforms);
nd_time = linspace(0, (size(nd_spike_waveforms, 2)-1)/30e3, size(nd_spike_waveforms, 2));
orig_mean = mean(orig_spike_waveforms); orig_std = std(orig_spike_waveforms);
orig_time = linspace(0, (size(orig_spike_waveforms, 2)-1)/30e3, size(orig_spike_waveforms, 2));

figure()
tcl = tiledlayout(1, 2, "TileSpacing", "tight");
ax(1) = nexttile();
c1 = nd_mean + nd_std;
c2 = nd_mean - nd_std;
t2 = [nd_time, fliplr(nd_time)];
a1 = [c1, fliplr(c2)];
h = fill(t2*1000, a1, 'k'); hold on
set(h, 'facealpha', 0.5);
set(h, 'facecolor', '#33bbee');
plot(nd_time*1000, nd_mean, 'LineWidth', 2, "Color", "#262626");
ylabel("Voltage (mV)")
xlabel("Time (ms)")
xlim([0, 5.5])
xticks(0:1:5)
legend("Standard Deviation", "Mean", "Location", "Northeast")
title("NeuroDAC output")

ax(2) = nexttile();
c1 = orig_mean + orig_std;
c2 = orig_mean - orig_std;
t2 = [orig_time, fliplr(orig_time)];
a1 = [c1, fliplr(c2)];
h = fill(t2*1000, a1, 'k'); hold on
set(h, 'facealpha', 0.5);
set(h, 'facecolor', '#33bbee');
plot(orig_time*1000, orig_mean, 'LineWidth', 2, "Color", "#262626");
legend("Standard Deviation", "Mean", "Location", "Northeast")
title("Original file")
xlabel("Time (ms)")
xlim([0, 5.5])
xticks(0:1:5)
linkaxes(ax, 'xy');
title(tcl, "NeuroDAC Performance Assessment")

set(gcf, "Position", [300, 300, 1254, 420])