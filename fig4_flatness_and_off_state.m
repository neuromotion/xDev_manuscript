%% fig4_flatness_and_off_state.m
% Author: Xavier Lee
% Date: 9/29/2024
% (C) Brown Neuromotion Lab 2024
%
% PROVIDED WITHOUT WARRANTY OR GUARANTEE
clear; close all; clc
load("data\Fig4_Ron_Flatness.mat")

%% On resistance flatness
close all
figure()
tiledlayout(2,2);

nexttile([1 2])
plot(Tc_ba.AppliedVoltage, Tc_ba.Rmux, "Color", "#b04cb9"); hold on
plot(Tb.AppliedVoltage, Tb.Rmux, "Color", "#06D6A0");
plot(Ta.AppliedVoltage, Ta.Rmux, "Color", "#33bbee"); 

xline(-12, 'Label', 'VSS')
xline(12, 'Label', 'VDD')
ylabel('Resistance (Ω)')
xlabel('Applied Voltage (V)')
xlim([-15,15])
legend('Module 3', 'Module 2', 'Module 1', 'Location','southeast')
title('xDev Connection Resistance over Channel Voltage, B9-A6')
hold off
ylim([100, 200])

nexttile([1 2])
plot(Tc_ba.AppliedVoltage, Tc_ba.Rmux, "Color", "#b04cb9"); hold on
plot(Tc_ab.AppliedVoltage, Tc_ab.Rmux, "Color", "#FD833D");
scatter(Tr.Voltage, Tr.rmux, "filled", "MarkerFaceColor", "#bbbbbb", "MarkerEdgeColor", "none");

hold off
xline(-12, 'Label', 'VSS')
xline(12, 'Label', 'VDD')
ylabel('Resistance (Ω)')
xlabel('Applied Voltage (V)')
legend('Module 3 B-to-A B9-A6', 'Module 3 A-to-B B9-A6', 'Random Samples', 'Location','southeast')
title('Module 3 Variability, B9-A6')
ylim([100, 200])

%% Off resistance
Off_resistance.xDev_voltage = Off_resistance.input_voltage - Off_resistance.resistor_voltage;
Off_resistance.Current_amps = Off_resistance.resistor_voltage ./ Off_resistance.R_sense_ohms;
Off_resistance.xDev_resistance = Off_resistance.xDev_voltage ./ Off_resistance.Current_amps;

figure()
plot(Off_resistance.input_voltage, Off_resistance.xDev_resistance ./ 1e6)
title("Off resistance")
ylabel("Resistance (M\Omega{})")
xlabel("Applied voltage (V)")
