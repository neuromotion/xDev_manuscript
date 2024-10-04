%% fig3_AC_DC_Impedance.m
% Author: Sam Parker
% Date: 8/21/2024
% (C) Brown Neuromotion Lab 2024
%
% PROVIDED WITHOUT WARRANTY OR GUARANTEE
clear; close all; clc
load("data\Fig3_resistance.mat")
do_dc = true;
do_ac = true;
rng(320923539); % to be consistent, set to seed used to generate figures: 320923539

%% DC
if do_dc
%% Plot histograms for each board

figure()
tcl = tiledlayout(1, 3, "TileSpacing", "tight");
for boardID = 1:3
    ax(boardID) = nexttile;
    histogram(dc_impedance.(sprintf("board%d_resistance_R", boardID)));
    xlabel("Line resistance (Ohms)")
    ylabel("Number of lines")
    title(sprintf("Board %d", boardID))
    subtitle(sprintf("(%0.2f ± %0.2f)Ω", nanmean(dc_impedance.(sprintf("board%d_resistance_R", boardID))), nanstd(dc_impedance.(sprintf("board%d_resistance_R", boardID)))));
end
linkaxes(ax, 'xy');
title(tcl, "Line DC Impedance")
subtitle(tcl, sprintf("n = %d, N = %d", length(dc_impedance.config_JSON), max(dc_impedance.config_JSON)))
set(gcf, "Position", [2383, 433, 1093, 495])

%% Plot boxplots for each board
figure()
boxplot(dc_impedance{:, ["board1_resistance_R", "board2_resistance_R", "board3_resistance_R"]})
ylabel("Line DC Impedance (Ohms)")
xlabel("Board number")
title("Line DC Impedance by Board")
ylim([99, 121])

%% Stats 
[p, tbl, stats] = kruskalwallis(dc_impedance{:, ["board1_resistance_R", "board2_resistance_R", "board3_resistance_R"]});
c = multcompare(stats, "Display","off");
for i = 1:size(c, 1)
    if c(i, 6) <= 0.05
        fprintf("Board %d is sig diff from Board %d (p=%0.4f)!\n", c(i, 1), c(i, 2), c(i, 6))
    else
        fprintf("Board %d is not sig diff from Board %d (p=%0.4f)!\n", c(i, 1), c(i, 2), c(i, 6))
    end
end

%% Impedance mismatch between random routings

numIterations = 10000;
dc_impedance_mismatch = nan(numIterations, 3);
dc_impedance_mismatch_vals = nan(numIterations, 6);

i = 1;
disp("Iterating...")
while i <= numIterations
    idxs = randi([1, height(dc_impedance)], [1, 2]);
    if dc_impedance.X_Gate(idxs(1)) == dc_impedance.X_Gate(idxs(2)) || dc_impedance.Y_Gate(idxs(1)) == dc_impedance.Y_Gate(idxs(2))
        continue
    end

    for boardID = 1:3
        dc_impedance_mismatch(i, boardID) = dc_impedance.(sprintf("board%d_resistance_R", boardID))(idxs(1)) - dc_impedance.(sprintf("board%d_resistance_R", boardID))(idxs(2));
        dc_impedance_mismatch_vals(i, (boardID-1)*2+1) = dc_impedance.(sprintf("board%d_resistance_R", boardID))(idxs(1));
        dc_impedance_mismatch_vals(i, (boardID-1)*2+2) = dc_impedance.(sprintf("board%d_resistance_R", boardID))(idxs(2));
    end
    i = i + 1;
end
disp("Done!")

mean_mismatch = nanmean(dc_impedance_mismatch);
std_mismatch = nanstd(dc_impedance_mismatch);

nbins = 17;
clear ax
figure()
tiledlayout(1, 3)
ax(1) = nexttile();
h = histogram(dc_impedance_mismatch(:, 1), nbins);
title("Board 1")
subtitle(sprintf("(%0.4f±%0.2f)Ω", mean_mismatch(1), std_mismatch(1)));
ylabel("Number of differential pairs")
xlabel("ΔR (Ω)")
ax(2) = nexttile();
histogram(dc_impedance_mismatch(:, 2), h.BinEdges);
title("Board 2")
subtitle(sprintf("(%0.4f±%0.2f)Ω", mean_mismatch(2), std_mismatch(2)));
xlabel("ΔR (Ω)")
ax(3) = nexttile();
histogram(dc_impedance_mismatch(:, 3), h.BinEdges);
title("Board 3")
subtitle(sprintf("(%0.4f±%0.2f)Ω", mean_mismatch(3), std_mismatch(3)));
xlabel("ΔR (Ω)")
linkaxes(ax, 'xy')
set(gcf, "Position", [2383, 433, 1093, 495])


[p, tbl, stats] = kruskalwallis(dc_impedance_mismatch);
c = multcompare(stats, "Display","off");
for i = 1:size(c, 1)
    if c(i, 6) <= 0.05
        fprintf("Board %d is sig diff from Board %d (p=%0.4f)!\n", c(i, 1), c(i, 2), c(i, 6))
    else
        fprintf("Board %d is not sig diff from Board %d (p=%0.4f)!\n", c(i, 1), c(i, 2), c(i, 6))
    end
end

figure()
boxplot(dc_impedance_mismatch)
ylabel("DC Impedance Mismatch (Ohms)")
xlabel("Board number")
title("Impedance Mismatch by Board")
end

%% AC Impedance
if do_ac
clearvars -except ac_impedance

% Calculate unknown impedance using voltage divider
ac_impedance.Z_xDev = ac_impedance.Output_mean_Vpp ./ (ac_impedance.Input_mean_Vpp - ac_impedance.Output_mean_Vpp) .* ac_impedance.R_Ref_Ohm;

numIterations = 10000;
ac_impedance_mismatch.freq = nan(numIterations, 1);
ac_impedance_mismatch.delta_mag = nan(numIterations, 1);
ac_impedance_mismatch.delta_phase = nan(numIterations, 1);
ac_impedance_mismatch.mag1 = nan(numIterations, 1);
ac_impedance_mismatch.mag2 = nan(numIterations, 1);
ac_impedance_mismatch.phase1 = nan(numIterations, 1);
ac_impedance_mismatch.phase2 = nan(numIterations, 1);
ac_impedance_mismatch.A_Gate1 = repmat("", numIterations, 1);
ac_impedance_mismatch.B_Gate1 = repmat("", numIterations, 1);
ac_impedance_mismatch.A_Gate2 = repmat("", numIterations, 1);
ac_impedance_mismatch.B_Gate2 = repmat("", numIterations, 1);

unique_gate_combos = unique(ac_impedance(:, 1:2), "rows");
unique_freqs = unique(ac_impedance.Freq_Hz);

i = 1;
disp("Iterating...")
while i <= numIterations
    gate_idxs = randi([1, height(unique_gate_combos)], [1, 2]);
    % check if same combo has been selected
    if all(unique_gate_combos{gate_idxs(1), :} == unique_gate_combos{gate_idxs(2), :})
        continue
    end

    ac_impedance_mismatch.A_Gate1(i) = unique_gate_combos{gate_idxs(1), 1};
    ac_impedance_mismatch.B_Gate1(i) = unique_gate_combos{gate_idxs(1), 2};
    ac_impedance_mismatch.A_Gate2(i) = unique_gate_combos{gate_idxs(2), 1};
    ac_impedance_mismatch.B_Gate2(i) = unique_gate_combos{gate_idxs(2), 2};
    ac_impedance_mismatch.freq(i) = unique_freqs(randi([1, length(unique_freqs)], 1));
    ac_impedance_mismatch.mag1(i) = mean(ac_impedance{ac_impedance.A_Gate == unique_gate_combos{gate_idxs(1), 1} & ...
                                                 ac_impedance.B_Gate == unique_gate_combos{gate_idxs(1), 2} & ...
                                                 ac_impedance.Freq_Hz == ac_impedance_mismatch.freq(i), "Z_xDev"});
    ac_impedance_mismatch.mag2(i) = mean(ac_impedance{ac_impedance.A_Gate == unique_gate_combos{gate_idxs(2), 1} & ...
                                                 ac_impedance.B_Gate == unique_gate_combos{gate_idxs(2), 2} & ...
                                                 ac_impedance.Freq_Hz == ac_impedance_mismatch.freq(i), "Z_xDev"});
    ac_impedance_mismatch.phase1(i) = mean(ac_impedance{ac_impedance.A_Gate == unique_gate_combos{gate_idxs(1), 1} & ...
                                                 ac_impedance.B_Gate == unique_gate_combos{gate_idxs(1), 2} & ...
                                                 ac_impedance.Freq_Hz == ac_impedance_mismatch.freq(i), "Phase_mean_deg"});
    ac_impedance_mismatch.phase2(i) = mean(ac_impedance{ac_impedance.A_Gate == unique_gate_combos{gate_idxs(2), 1} & ...
                                                 ac_impedance.B_Gate == unique_gate_combos{gate_idxs(2), 2} & ...
                                                 ac_impedance.Freq_Hz == ac_impedance_mismatch.freq(i), "Phase_mean_deg"});
    i = i + 1;
end
disp("Done!")

ac_impedance_mismatch = struct2table(ac_impedance_mismatch);
ac_impedance_mismatch.delta_mag = ac_impedance_mismatch.mag1 - ac_impedance_mismatch.mag2;
ac_impedance_mismatch.delta_phase = ac_impedance_mismatch.phase1 - ac_impedance_mismatch.phase2;

unique_freqs = unique(ac_impedance_mismatch.freq);
figure()
tcl = tiledlayout(2, length(unique_freqs));
for freqIdx = 1:length(unique_freqs)
    mag_ax(freqIdx) = nexttile();
    if freqIdx == 1
        mag_h = histogram(ac_impedance_mismatch.delta_mag(ac_impedance_mismatch.freq == unique_freqs(freqIdx)));
    else
        histogram(ac_impedance_mismatch.delta_mag(ac_impedance_mismatch.freq == unique_freqs(freqIdx)), mag_h.BinEdges);
    end
    title(sprintf("%dHz", unique_freqs(freqIdx)))
end
linkaxes(mag_ax);

for freqIdx = length(unique_freqs):-1:1
    phase_ax(freqIdx) = nexttile();
    if freqIdx == length(unique_freqs)
        phase_h = histogram(ac_impedance_mismatch.delta_phase(ac_impedance_mismatch.freq == unique_freqs(freqIdx)));
    else
        histogram(ac_impedance_mismatch.delta_phase(ac_impedance_mismatch.freq == unique_freqs(freqIdx)), phase_h.BinEdges);
    end
    title(sprintf("%dHz", unique_freqs(freqIdx)))
end
linkaxes(phase_ax);

mag_counts = nan(length(mag_h.BinEdges)-1, length(unique_freqs));
phase_counts = nan(length(phase_h.BinEdges)-1, length(unique_freqs));
for freqIdx = 1:length(unique_freqs)
    mag_counts(:, freqIdx) = histcounts(ac_impedance_mismatch.delta_mag(ac_impedance_mismatch.freq == unique_freqs(freqIdx)), mag_h.BinEdges);
    phase_counts(:, freqIdx) = histcounts(ac_impedance_mismatch.delta_phase(ac_impedance_mismatch.freq == unique_freqs(freqIdx)), phase_h.BinEdges);
end

mag_counts_max = max(mag_counts, [], "all"); mag_counts_min = min(mag_counts, [], "all");
mag_counts_norm = (mag_counts - mag_counts_min) / (mag_counts_max - mag_counts_min);
mag_bin_centers = mag_h.BinEdges(1:end-1) + mag_h.BinWidth/2;

phase_counts_max = max(phase_counts, [], "all"); phase_counts_min = min(phase_counts, [], "all");
phase_counts_norm = (phase_counts - phase_counts_min) / (phase_counts_max - phase_counts_min);
phase_bin_centers = phase_h.BinEdges(1:end-1) + phase_h.BinWidth/2;

figure("Renderer", "Painters")
tcl = tiledlayout(2, 1, "TileSpacing","tight", "Padding","compact");
clear ax
ax(1) = nexttile();
surf(unique_freqs, mag_bin_centers, mag_counts_norm, "EdgeColor", "none");
set(gca, "XScale", "log")
xlabel("Freq (Hz)")
ylabel("Δ|Z| (Ω)")
view([0, 90])
ylim([min(mag_bin_centers), max(mag_bin_centers)])

ax(2) = nexttile();
surf(unique_freqs, phase_bin_centers, phase_counts_norm, "EdgeColor", "none");
set(gca, "XScale", "log")
xlabel("Freq (Hz)")
ylabel("Δ∠Z (deg)")
view([0, 90])
ylim([min(phase_bin_centers), max(phase_bin_centers)])

linkaxes(ax, 'x');
xlim([min(unique_freqs), max(unique_freqs)])
title(tcl, "Impedance Mismatch Spectrogram")
set(gcf, "Position", [112, 215, 1720, 750])
colormap("abyss")
colorbar

end