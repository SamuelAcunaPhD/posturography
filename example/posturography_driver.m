% Filename: posturography_driver.m
% Author:   Samuel Acuña
% Date:     12 Nov 2018
% Description:
% sample posturography analysis

clear; close all; clc;
addpath('..') % be able to access 'posturography.m'

%% load COP data
%filenames = {'eyesOpen_1.csv'};
filenames = {'eyesOpen_1.csv';'eyesOpen_2.csv';'eyesOpen_3.csv';'eyesClosed_1.csv';'eyesClosed_2.csv';'eyesClosed_3.csv'};

COP = posturography.load(filenames); % loads data files, filters, detrends

%% statokinesiogram
% posturography.statokinesiogram(COP(1)); % plot single statokinesiogram
% posturography.statokinesiogram_animation(COP(1)); % plot single animated statokinesiogram
posturography.statokinesiogram(COP,[1 2 3; 4 5 6]); % plot multiple statokinesiograms in specified order
%posturography.saveFigure(fig,2)

%% stabilogram
% posturography.stabilogram(COP(1)); % plot single stabilogram
% posturography.stabilogram_animation(COP(1)); % plot single animated stabilogram
posturography.stabilogram(COP,[1 2 3; 4 5 6]); % plot multiple stabilograms in specified order

%% estimate Center of Gravity
posturography.estimateCOG(COP,[1 2 3; 4 5 6])

%% posturography metrics
% [~,metrics] = posturography.calcMetrics(COP(1)) % just output metrics
COP = posturography.calcMetrics(COP) % add metrics to COP structure

%% plot posturography metrics
posturography.plotMetrics(COP,[1 2 3; 4 5 6]);
% posturography.saveFigure(fig,2)

%% plot FFTs
% posturography.plotSpectralAnalysis(COP(1)) % single FFT plot
posturography.plotSpectralAnalysis(COP,[1 2 3; 4 5 6]) % all FFTs



