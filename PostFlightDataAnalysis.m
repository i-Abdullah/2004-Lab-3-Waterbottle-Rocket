%% info:

%{

This code analyzes the data obtained from an altiemeter attached with the
water-bottle rocket. The data is stored in folder called LaunchData, it
contains temp, pressure, and time, we can use this with International
Standard Atmosphere model in matlab to see how high the rockets actually
went.


%}

%% housekeeping

clear;
clc;

close all;

%% Launch 1

Data = importdata('LaunchData/Launch1/group5_8am_launch1_processed_smoothed.dms');

Data = Data.data; % find the structure thing.
Time = Data(:,1);
Pressure = Data(:,2)*100; % convert from milliebar to pasacal
Temp = Data(:,3);


% since pressure changes linearly, we can find min, and starting
% pressure and this's our height change.


BoulderAlt = 1600:0.01:1700;

% get readings for different altitudes:
[T, a, P, rho] = atmosisa(BoulderAlt);

% find altitude corresponding to min pressure (top) and starting pressure

MinP = min(Pressure);
StartP = Pressure(1);

diffMin = abs(MinP-P);
index1 = find(min(diffMin)==diffMin);


diffStart = abs(StartP-P);
index2 = find(min(diffStart)==diffStart);

dh = BoulderAlt(index1) - BoulderAlt(index2) ; % height change



%% read data

