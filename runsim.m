M = 16;                     % Size of signal constellation
k = log2(M);                % Number of bits per symbol
n = 40000;                  % Number of bits to process
numSamplesPerSymbol = 1;    % Oversampling factor

% rng default                 % Use default random number generator
% dataIn = randi([0 1],n,1);  % Generate vector of binary data

load('bits.mat')

dataIn = bits';

% stem(dataIn(1:40),'filled');
% title('Random Bits');
% xlabel('Bit Index');
% ylabel('Binary Value');

% figure; % Create new figure window.
% stem(dataSymbolsIn(1:10));
% title('Random Symbols');
% xlabel('Symbol Index');
% ylabel('Integer Value');

[encoded_td, qamout] = qamenc(dataIn);
power = (1/length(encoded_td))*sum(encoded_td.^2)
powerperc = power/0.00125

%afterchan = encoded_td;

afterchan = chansim(encoded_td);

audiowrite('tx.wav', encoded_td, 44100, 'BitsPerSample', 24);

Fs = 44100;
%system('ccplay tx.wav rx.wav --prepause 0.27 --channel audio0 --depth 24 --rate 44100');

[afterchan, ~] = audioread('rx.wav');

% sPlotFig = scatterplot(receivedSignalG,1,0,'g.');
% hold on
% scatterplot(dataModG,1,0,'k*',sPlotFig)

[dataOutG, extracted] = qamdec(afterchan, length(encoded_td), qamout, dataIn);

[numErrorsG,berG] = biterr(dataIn,dataOutG);
fprintf('\nThe Gray coding bit error rate = %5.2e, based on %d errors\n', ...
    berG,numErrorsG)

correct = sum(dataOutG == dataIn)
R = length(bits)/(length(encoded_td)/Fs);
N = (length(bits) - correct)*200000/length(bits);
P = power;
fom = (min(R, 3000000)*(1-N/100000)^10)/max(1,800*P)