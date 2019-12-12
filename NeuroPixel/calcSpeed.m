function [speed,raw_speed] = calcSpeed(posx,p)
% function to calculate speed from position and time data
% Malcolm Campbell 5/21/15
%
% modified 6/6/18 MGC
% all time bins are now equal length
% has not been fully tested yet
%
% inputs:
%     posx: positions
%     p: params (spatial bin size, etc)
%     smoothingSigma: sigma for smoothing of running speed trace
% outputs:
%     speed: smoothed running speed

% smoothing parameter for running speed trace (in time bins)
smoothSigma = 10;

% calculate raw speed
speed = diff(posx)/p.TimeBin;

% throw out extreme values and interpolate
speed(speed > 150) = NaN;
speed(speed<-5) = NaN;
speed(isnan(speed)) = interp1(find(~isnan(speed)), speed(~isnan(speed)), find(isnan(speed)), 'pchip'); % interpolate NaNs
speed = [0;speed];
raw_speed = speed;
% smooth speed trace
speed = gauss_smoothing(speed,smoothSigma);

end