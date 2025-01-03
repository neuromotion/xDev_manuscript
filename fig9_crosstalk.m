%% fig9_crosstalk.m
% Author: Sam Parker
% Date: 8/16/2024
% (C) Brown Neuromotion Lab 2024
%
% PROVIDED WITHOUT WARRANTY OR GUARANTEE
clear; close all; clc
load("data\Fig9_crosstalk.mat")

%% Calculate amplifier gain in dB

amp_gain_characterization.gain = amp_gain_characterization.output_rms ./ amp_gain_characterization.input_rms_divided;
amp_gain_characterization.gain_dB = 20 * log10(amp_gain_characterization.output_rms ./ amp_gain_characterization.input_rms_divided);

crosstalk.amp_gain = nan(height(crosstalk), 1);
for i = 1:height(crosstalk)
    % amplifier only used for lower frequencies
    if crosstalk.freq(i) > 10e3
        crosstalk.amp_gain(i) = 1;
    else
        gain_entry_idx = find(amp_gain_characterization.freq == crosstalk.freq(i));
        if isempty(gain_entry_idx)
            error("Frequency for gain calculation is not in calibration curve!")
        end
        crosstalk.amp_gain(i) = amp_gain_characterization.output_rms(gain_entry_idx) / amp_gain_characterization.input_rms_divided(gain_entry_idx);
    end
end
crosstalk.output_rms_unamplified = crosstalk.output_rms_mean ./ crosstalk.amp_gain;

%% Calculate crosstalk
crosstalk.attentuation_dB = 20 * log10(crosstalk.input_rms_mean./ crosstalk.output_rms_unamplified);

%% Plot gain vs aggressor frequency

figure()
semilogx(crosstalk.freq, crosstalk.attentuation_dB); hold on
xlabel("Aggressor frequency (Hz)")
ylabel("Crosstalk Attenuation (dB)")
title("Crosstalk power vs aggressor frequency")
ylim([20, 100])

%% Plot time-domain examples

f_list = string(fieldnames(time_domain_examples));
f_list = split(f_list, "_"); f_list(:, 1) = [];
f_list = split(f_list, "Hz"); f_list(:, 2) = [];
f_list = double(f_list);

figure()
tcl = tiledlayout(ceil(length(f_list)/2), 2);
for f_idx = 1:length(f_list)
    fieldname = sprintf("f_%dHz", f_list(f_idx));
    if f_list(f_idx) <= 10e3
        time_domain_examples.(fieldname).Ch1_V = time_domain_examples.(fieldname).Ch1_V ./ amp_gain_characterization.gain(amp_gain_characterization.freq == f_list(f_idx));
    end

    % Ensure 4 cycles are plotted
    sample_rate = 1/mean(diff(time_domain_examples.(fieldname).Time_s));
    samples_per_cycle = sample_rate / f_list(f_idx);
    zero_crossing_sample = find(time_domain_examples.(fieldname).Time_s >= 0, 1)-1;
    end_sample = zero_crossing_sample + samples_per_cycle*2;
    start_sample = zero_crossing_sample - samples_per_cycle*2;


    nexttile();
    plot(time_domain_examples.(fieldname).Time_s(round(start_sample:end_sample)), time_domain_examples.(fieldname).Ch1_V(round(start_sample:end_sample))); hold on
    plot(time_domain_examples.(fieldname).Time_s(round(start_sample:end_sample)), time_domain_examples.(fieldname).Ch2_V(round(start_sample:end_sample))); hold on
    xlabel("Time (s)")
    ylabel("Voltage (V)")
    title(sprintf("%dHz", f_list(f_idx)))
    xlim([time_domain_examples.(fieldname).Time_s(round(start_sample)), time_domain_examples.(fieldname).Time_s(round(end_sample))]);
end

set(gcf, "Position", [100, 100, 560, 750]);

legend("Victim", "Aggressor", "Location", "southoutside", "Orientation","horizontal");

%% Plot single time domain example on yyaxis

% Ensure 4 cycles are plotted
sample_rate = 1/mean(diff(time_domain_examples.f_10000Hz.Time_s));
samples_per_cycle = sample_rate / 10000;
zero_crossing_sample = find(time_domain_examples.f_10000Hz.Time_s >= 0, 1)-1;
end_sample = zero_crossing_sample + samples_per_cycle*2;
start_sample = zero_crossing_sample - samples_per_cycle*2;

figure()
yyaxis left
plot((time_domain_examples.f_10000Hz.Time_s(round(start_sample:end_sample)) - time_domain_examples.f_10000Hz.Time_s(round(start_sample))).*1000, time_domain_examples.f_10000Hz.Ch1_V(round(start_sample:end_sample)));
ylabel("Victim (V)")
ylim([-5e-3, 5e-3])
yyaxis right
plot((time_domain_examples.f_10000Hz.Time_s(round(start_sample:end_sample)) - time_domain_examples.f_10000Hz.Time_s(round(start_sample))).*1000, time_domain_examples.f_10000Hz.Ch2_V(round(start_sample:end_sample)));
ylabel("Aggressor (V)")
xlabel("Time (ms)")
title("10000Hz Crosstalk")
set(gcf, "Position", [100, 100, 560, 420])

%% Square wave examples

% Single example
time_s = time_domain_examples.square_waves.a8_v7.f_10000Hz.Time_s - time_domain_examples.square_waves.a8_v7.f_10000Hz.Time_s(1);
figure()
tcl = tiledlayout(2, 1, "TileSpacing", "tight");
nexttile()
plot(time_s*1e6, time_domain_examples.square_waves.a8_v7.f_10000Hz.Ch1_V); hold on
plot(time_s*1e6, time_domain_examples.square_waves.a8_v7.f_10000Hz.Ch2_V);
ylabel("Voltage (V)")
ylim([-1.1, 1.1])
yticks([-1:0.5:1])
nexttile()
yyaxis left; plot(time_s*1e6, time_domain_examples.square_waves.a8_v7.f_10000Hz.Ch1_V .* 1000);
ylabel("Victim (mV)")
ylim([-50, 50])
yticks([-50:25:50])
yyaxis right; plot(time_s*1e6, time_domain_examples.square_waves.a8_v7.f_10000Hz.Ch2_V);
ylim([-1.5, 1.5])
yticks([-1.5:0.75:1.5])
ylabel("Aggressor (V)")
xlabel("Time (μs)")
xlim([0, 200])
xticks([0:50:200])
title(tcl, "Square wave crosstalk at 10kHz")
set(gcf, "Position", [500, 500, 560, 420])

% All frequencies
frequencies = string(fieldnames(time_domain_examples.square_waves.a8_v7));
figure()
tcl = tiledlayout(2, length(frequencies), "TileSpacing", "tight", "Padding", "tight");
for i = 1:length(frequencies)
    this_frequency = split(frequencies(i), "f_"); this_frequency = this_frequency(2);
    this_frequency = split(this_frequency, "Hz"); this_frequency = str2double(this_frequency(1));
    nexttile()
    plot(time_domain_examples.square_waves.a8_v7.(frequencies(i)).Time_s - time_domain_examples.square_waves.a8_v7.(frequencies(i)).Time_s(1), time_domain_examples.square_waves.a8_v7.(frequencies(i)).Ch1_V); hold on
    plot(time_domain_examples.square_waves.a8_v7.(frequencies(i)).Time_s - time_domain_examples.square_waves.a8_v7.(frequencies(i)).Time_s(1), time_domain_examples.square_waves.a8_v7.(frequencies(i)).Ch2_V);
    ylabel("Voltage (V)")
    title(sprintf("%d Hz", this_frequency))
end
for i = 1:length(frequencies)
    this_frequency = split(frequencies(i), "f_"); this_frequency = this_frequency(2);
    this_frequency = split(this_frequency, "Hz"); this_frequency = str2double(this_frequency(1));
    nexttile()
    yyaxis left; plot(time_domain_examples.square_waves.a8_v7.(frequencies(i)).Time_s - time_domain_examples.square_waves.a8_v7.(frequencies(i)).Time_s(1), time_domain_examples.square_waves.a8_v7.(frequencies(i)).Ch1_V .* 1000);
    ylabel("Victim (mV)")
    ylim([-50, 50])
    yticks([-50:25:50])
    yyaxis right; plot(time_domain_examples.square_waves.a8_v7.(frequencies(i)).Time_s - time_domain_examples.square_waves.a8_v7.(frequencies(i)).Time_s(1), time_domain_examples.square_waves.a8_v7.(frequencies(i)).Ch2_V);
    ylabel("Aggressor (V)")
    ylim([-1.5, 1.5])
    yticks([-1.5:0.75:1.5])
    xlabel("Time (μs)")
end
title(tcl, "Square wave crosstalk");
set(gcf, "Position", [61, 426, 1852, 420])

%% Calculate timeconstant

frequencies = string(fieldnames(time_domain_examples.square_waves.a8_v7));
freq_num = nan(length(frequencies), 1);
tc = nan(size(freq_num));

figure()
tcl = tiledlayout(length(frequencies), 1, "TileSpacing", "tight");
for i = 1:length(frequencies)-2 % not using last 2 frequencies, as half-period is shorter than time constant
    this_frequency = split(frequencies(i), "f_"); this_frequency = this_frequency(2);
    this_frequency = split(this_frequency, "Hz"); this_frequency = str2double(this_frequency(1));
    [v0, peak_idx] = max(time_domain_examples.square_waves.a8_v7.(frequencies(i)).Ch1_V(1:round(length(time_domain_examples.square_waves.a8_v7.(frequencies(i)).Ch1_V)/4*3)));
    t1_lvl = v0 / exp(1);
    tc_idx = find(time_domain_examples.square_waves.a8_v7.(frequencies(i)).Ch1_V <= t1_lvl & time_domain_examples.square_waves.a8_v7.(frequencies(i)).Time_s >= time_domain_examples.square_waves.a8_v7.(frequencies(i)).Time_s(peak_idx), 1);
    freq_num(i) = this_frequency;
    tc(i) = time_domain_examples.square_waves.a8_v7.(frequencies(i)).Time_s(tc_idx) - time_domain_examples.square_waves.a8_v7.(frequencies(i)).Time_s(peak_idx);
    
end

fprintf("Time constant for square wave transient artifact is %f microseconds\n", nanmean(tc) * 1e6)

%% Filtering victim channel for 30kHz recording

figure()
tcl = tiledlayout(2, 1, "TileSpacing", "tight");
nexttile()
plot(time_domain_examples.square_waves.a8_v7.f_10000Hz.Time_s, time_domain_examples.square_waves.a8_v7.f_10000Hz.Ch1_V); hold on
plot(time_domain_examples.square_waves.a8_v7.f_10000Hz.Time_s, lowpass(time_domain_examples.square_waves.a8_v7.f_10000Hz.Ch1_V, 15e3, 1/nanmean(diff(time_domain_examples.square_waves.a8_v7.f_10000Hz.Time_s)), "Steepness", 0.95, "StopbandAttenuation", 95));
nexttile()
yyaxis left; plot(time_domain_examples.square_waves.a8_v7.f_10000Hz.Time_s, time_domain_examples.square_waves.a8_v7.f_10000Hz.Ch1_V);
yyaxis right; plot(time_domain_examples.square_waves.a8_v7.f_10000Hz.Time_s, lowpass(time_domain_examples.square_waves.a8_v7.f_10000Hz.Ch1_V, 15e3, 1/nanmean(diff(time_domain_examples.square_waves.a8_v7.f_10000Hz.Time_s)), "Steepness", 0.95, "StopbandAttenuation", 95));