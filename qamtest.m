M = 16;                     % Size of signal constellation
k = log2(M);                % Number of bits per symbol
n = 36000;                  % Number of bits to process
numSamplesPerSymbol = 1;    % Oversampling factor

rng default                 % Use default random number generator
dataIn = randi([0 1],n,1);  % Generate vector of binary data

dataInMatrix = reshape(dataIn,length(dataIn)/k,k);   % Reshape data into binary k-tuples, k = log2(M)
dataSymbolsIn = bi2de(dataInMatrix);                 % Convert to integers

% stem(dataIn(1:40),'filled');
% title('Random Bits');
% xlabel('Bit Index');
% ylabel('Binary Value');

% figure; % Create new figure window.
% stem(dataSymbolsIn(1:10));
% title('Random Symbols');
% xlabel('Symbol Index');
% ylabel('Integer Value');

refconst = qammod(0:M-1,M);
dataModG = qammod(dataSymbolsIn, M, 0/180*pi);%.*exp(-1i*10/180*pi); % Gray coding, phase offset = 10 deg

EbNo = 20;
snr = EbNo + 10*log10(k) - 10*log10(numSamplesPerSymbol);

receivedSignalG = awgn(dataModG,snr,'measured');

sPlotFig = scatterplot(receivedSignalG,1,0,'g.');
hold on
scatterplot(dataModG,1,0,'k*',sPlotFig)
scatterplot(refconst,1,0,'r*',sPlotFig)
x = (0:15)'; 
text(real(refconst)+0.1, imag(refconst), dec2bin(x), 'fontsize', 18)
title('QAM-16, Gray Coded Mapping')
axis([-4 4 -4 4])
set(gca,'fontsize',18)
set(gcf,'color','w');

dataSymbolsOutG = qamdemod(receivedSignalG,M);

[numErrorsG,berG] = biterr(dataSymbolsIn,dataSymbolsOutG);
fprintf('\nThe Gray coding bit error rate = %5.2e, based on %d errors\n', ...
    berG,numErrorsG)