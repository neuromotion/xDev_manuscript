%% fig7_freq_response.m
% Author: Sam Parker
% Date: 8/12/2024
% (C) Brown Neuromotion Lab 2024
%
% PROVIDED WITHOUT WARRANTY OR GUARANTEE
clear; close all; clc
load("data\Fig7_freq_response.mat")

%% Calculate gain in dB

single_gate_connected_freq_response.gain_mean = single_gate_connected_freq_response.output_vpp_mean ./ single_gate_connected_freq_response.input_vpp_mean;
single_gate_connected_freq_response.gain_dB_mean = 20.*log10(single_gate_connected_freq_response.gain_mean);

random_gates_connected_frequency_response.gain_mean = random_gates_connected_frequency_response.output_vpp_mean ./ random_gates_connected_frequency_response.input_vpp_mean;
random_gates_connected_frequency_response.gain_dB_mean = 20.*log10(random_gates_connected_frequency_response.gain_mean);

%% Plot the frequency response of the single gate, with the random samples shown as scatters

figure()
tcl = tiledlayout(2, 1);
nexttile
semilogx(single_gate_connected_freq_response.frequency_hz, single_gate_connected_freq_response.gain_dB_mean); hold on
scatter(random_gates_connected_frequency_response.frequency_hz, random_gates_connected_frequency_response.gain_dB_mean, 'filled');
ylim([-12, 1])
title("Magnitude")
ylabel("|H| (dB)")
nexttile
semilogx(single_gate_connected_freq_response.frequency_hz, single_gate_connected_freq_response.phase_deg_mean); hold on
scatter(random_gates_connected_frequency_response.frequency_hz, random_gates_connected_frequency_response.phase_deg_mean, 'filled');
title("Phase")
ylabel("∠H (deg)")
xlabel("Frequency (Hz)")
title(tcl, "xDev Frequency Response")


%% Plot example traces to show magnitude and phase changes

frequencyFields = ["f_12_Hz", "f_1k1_Hz", "f_5_MHz", "f_11_MHz"];
frequencyStrings = ["12Hz", "1.1kHz", "5MHz", "11MHz"];

figure()
tcl = tiledlayout(1, length(frequencyFields), "TileSpacing", "tight", "Padding", "tight");
clear ax
for i = 1:length(frequencyFields)
    ax(i) = nexttile();
    plot(example_traces.(frequencyFields(i)).Time(1:2:end), example_traces.(frequencyFields(i)).Input(1:2:end)); hold on
    plot(example_traces.(frequencyFields(i)).Time(1:2:end), example_traces.(frequencyFields(i)).Output(1:2:end));

    % Add callout points for magnitude and phase shift
    plot([min(example_traces.(frequencyFields(i)).Time), max(example_traces.(frequencyFields(i)).Time)], [min(example_traces.(frequencyFields(i)).Input), min(example_traces.(frequencyFields(i)).Input)], "--k");
    plot([min(example_traces.(frequencyFields(i)).Time), max(example_traces.(frequencyFields(i)).Time)], [max(example_traces.(frequencyFields(i)).Input), max(example_traces.(frequencyFields(i)).Input)], "--k");
    plot([min(example_traces.(frequencyFields(i)).Time), max(example_traces.(frequencyFields(i)).Time)], [min(example_traces.(frequencyFields(i)).Output), min(example_traces.(frequencyFields(i)).Output)], "--k");
    plot([min(example_traces.(frequencyFields(i)).Time), max(example_traces.(frequencyFields(i)).Time)], [max(example_traces.(frequencyFields(i)).Output), max(example_traces.(frequencyFields(i)).Output)], "--k");

    exerpt_I = example_traces.(frequencyFields(i)).Input(example_traces.(frequencyFields(i)).Time <= (0.25*max(example_traces.(frequencyFields(i)).Time)) & example_traces.(frequencyFields(i)).Time > 0);
    exerpt_O = example_traces.(frequencyFields(i)).Output(example_traces.(frequencyFields(i)).Time <= (0.25*max(example_traces.(frequencyFields(i)).Time)) & example_traces.(frequencyFields(i)).Time > 0);
    exerpt_T = example_traces.(frequencyFields(i)).Time(example_traces.(frequencyFields(i)).Time <= (0.25*max(example_traces.(frequencyFields(i)).Time)) & example_traces.(frequencyFields(i)).Time > 0);
    [~, max_i_idx] = max(exerpt_I); [~, max_o_idx] = max(exerpt_O); 

    plot([exerpt_T(max_i_idx), exerpt_T(max_i_idx)], [0.8, 1], "--b");
    plot([exerpt_T(max_o_idx), exerpt_T(max_o_idx)], [0.8, 1], "--r");

    xlabel("Time [s]");
    if i == 1
        ylabel("Voltage (V)")
    end
    title(frequencyStrings(i));
end
linkaxes(ax, "y");
title(tcl, " Example Frequency Responses")
set(gcf, "Position", [210, 490, 1570, 420])
