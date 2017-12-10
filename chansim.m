function [ output ] = chansim( input )
%CHANSIM Summary of this function goes here
%   Detailed explanation goes here

    % Load impulse response
	load('IR0.mat','impulse');
    
    % Convolve input signal with impulse response
    afterfilter = conv(impulse', input);
    
    % Add white noise
    snr = 90;
    output = awgn(afterfilter, snr);
    
    output = output(1:length(input));
end

