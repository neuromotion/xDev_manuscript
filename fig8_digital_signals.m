%% Fig8_digital_signals.m
% Author: Sam Parker, B.E.
% Date: 7/22/2024
% (C) Brown Neuromotion Lab 2024
%
% Provided without guarantee or warranty
clear; close all; clc

%% Load the data
load("data\Fig8_ethernet.mat")

%% Summary
fprintf("Low side: (%0.2f ± %0.2f) (mean ± std_dev)\n", mean(ethernetSpeeds.lowSide.bits_per_second ./ 1e6), std(ethernetSpeeds.lowSide.bits_per_second ./ 1e6))
fprintf("High side: (%0.2f ± %0.2f) (mean ± std_dev)\n", mean(ethernetSpeeds.highSide.bits_per_second ./ 1e6), std(ethernetSpeeds.highSide.bits_per_second ./ 1e6))
fprintf("Cross over: (%0.2f ± %0.2f) (mean ± std_dev)\n", mean(ethernetSpeeds.crossOver.bits_per_second ./ 1e6), std(ethernetSpeeds.crossOver.bits_per_second ./ 1e6))
fprintf("All: (%0.2f ± %0.2f) (mean ± std_dev)\n", mean([ethernetSpeeds.crossOver.bits_per_second; ethernetSpeeds.highSide.bits_per_second; ethernetSpeeds.lowSide.bits_per_second] ./ 1e6), std([ethernetSpeeds.crossOver.bits_per_second; ethernetSpeeds.highSide.bits_per_second; ethernetSpeeds.lowSide.bits_per_second] ./ 1e6))
fprintf("Control: (%0.2f ± %0.2f) (mean ± std_dev)\n", mean(ethernetSpeeds.control.bits_per_second ./ 1e6), std(ethernetSpeeds.control.bits_per_second ./ 1e6))

%% Plot histograms

figure()
tcl = tiledlayout(1, 3, "TileSpacing","tight");
ax(1) = nexttile;
h = histogram(ethernetSpeeds.highSide.bits_per_second ./ 1e6); hold on
xlabel("Link speed (Mbps)")
ylabel("Number of intervals")
title("Lines X9:16 Y9:16")
subtitle(sprintf("%0.2f ± %0.2f", mean(ethernetSpeeds.highSide.bits_per_second ./ 1e6), std(ethernetSpeeds.highSide.bits_per_second ./ 1e6)));

ax(2) = nexttile;
histogram(ethernetSpeeds.lowSide.bits_per_second ./ 1e6, h.BinEdges); hold on
xlabel("Link speed (Mbps)")
title("Lines X0:8 Y0:8")
subtitle(sprintf("%0.2f ± %0.2f", mean(ethernetSpeeds.lowSide.bits_per_second ./ 1e6), std(ethernetSpeeds.lowSide.bits_per_second ./ 1e6)));

ax(3) = nexttile;
histogram(ethernetSpeeds.crossOver.bits_per_second ./ 1e6, h.BinEdges); hold on
xlabel("Link speed (Mbps)")
title("Lines X0:8 Y9:16")
subtitle(sprintf("%0.2f ± %0.2f", mean(ethernetSpeeds.crossOver.bits_per_second ./ 1e6), std(ethernetSpeeds.crossOver.bits_per_second ./ 1e6)));

linkaxes(ax, 'xy');
set(gcf, 'Position', [2086,  656, 1310, 420])

%% Plot boxplot
figure()
boxplot([ethernetSpeeds.highSide.bits_per_second ./ 1e6, ethernetSpeeds.lowSide.bits_per_second ./ 1e6, ethernetSpeeds.crossOver.bits_per_second ./ 1e6, ethernetSpeeds.control.bits_per_second ./ 1e6])
xticklabels(["Y9:16 Y9:16", "X0:8 Y0:8", "X0:8 Y9:16", "Control"])
ylabel("Link speed (Mbps)")
set(gcf, 'Position', [2258, 510, 382, 420])

%% Determine significance

[p, tbl, stats] = kruskalwallis([ethernetSpeeds.highSide.bits_per_second ./ 1e6, ethernetSpeeds.lowSide.bits_per_second ./ 1e6, ethernetSpeeds.crossOver.bits_per_second ./ 1e6, ethernetSpeeds.control.bits_per_second ./ 1e6], [], "off");
c = multcompare(stats, "Display","off");
group_names = ["High", "Low", "Cross", "Control"];
for i = 1:size(c, 1)
    if c(i, 6) <= 0.05
        fprintf("%s is sig diff from %s!\n", group_names(c(i, 1)), group_names(c(i, 2)))
    else
        fprintf("%s not sig diff from %s.\n", group_names(c(i, 1)), group_names(c(i, 2)))
    end
end
