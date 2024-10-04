%% fig5_programming_and_switching.m
% Author: Xavier Lee
% Date: 9/30/2024
% (C) Brown Neuromotion Lab 2024
%
% PROVIDED WITHOUT WARRANTY OR GUARANTEE
clear; close all; clc
load("data\Fig5_switching_time.mat")

%% Plot SPI programming time histogram

figure()
histogram(programming_times);
ylabel("Number of configurations")
xlabel("Programming time (s)")

%% Plot channel connection trace

figure()
for i = 1:length(switching_times.t_delta_s)
    plot(switching_times.time_trace(:, i), switching_times.pclk_trace(:, i), '-k'); hold on
    plot(switching_times.time_trace(:, i), switching_times.B1_trace(:, i), "Color", "#33bbee"); hold on
end

%% Plot channel connection time histogram
figure()
histogram(switching_times.t_delta_s*1e9, 71.8:0.2:73.2);
ylabel("Number of xDev Channels")
xlabel("Connection time (ns)")