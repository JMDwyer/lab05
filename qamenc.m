function [x, qamout] = qamenc( bits )
    % Don't need special sync symbols at the beginning because commcloud is
    % high SNR. Can get the start of the signal by inspecting power
    
    M = 32;                     % Size of signal constellation
    k = log2(M);                % Number of bits per symbol (QAM symbol)
    n = length(bits);           % Number of bits to process
    numSamplesPerSymbol = 1;    % Oversampling factor
    
    % Start freq index
    freq_lo = 12;
    
    % Stop freq index
    freq_hi = 4100;
    
    % Figure out how many bits we are sending per QAM symbol
    load('QAMbits.mat');
    
    % Find where the number of bits per QAM symbol changes. This is to
    % vectorize the qammod code to speed up performance. 
    qam_bits = qam_bits(freq_lo:freq_hi); % Truncate to used freq range
    qam_bits(qam_bits > 6) = 6; % Clip the qam bits
    qam_bits_idx = [0 find(diff(qam_bits)) length(qam_bits)];
    qam_bps = qam_bits(qam_bits_idx(2:end));
%     qam_bits = ones(1, freq_hi-freq_lo+1)*5;
%     qam_bits_idx = ceil(linspace(0, length(qam_bits), 100));
%     qam_bps = qam_bits(qam_bits_idx(2:end));
    
    % Number of bits per packet
    BPP = sum(diff(qam_bits_idx).*qam_bps);

    % Number of data packets per training packet
    DPPTP = 14;
    
    % The number of zeros to append in freq domain for a cut off of 18kHz
    ignore = ceil(freq_hi/18000*(22050-18000));
    ignore = 4900 - freq_hi;

    % Number of samples to prepend in time domain
    prepend = 200;

    % power constraint and number of bits in total
    numbits = length(bits);
    P = 0.00125/((BPP + prepend/2)/(BPP+ignore + prepend/2));
    P = 0.00125*22050/18000/0.67;

    % Number of packets to send and do any required padding
    packets = ceil(numbits/BPP);
    lastsympad = packets*BPP - numbits;
    bits = [bits; zeros(lastsympad, 1)];

    % Generate the random phases
    rng(4670);
    randphase = rand([freq_hi-freq_lo+1, 1]);

    % Create the training signal
    TR = [zeros(freq_lo-1,1); sqrt(P)*exp(1i*randphase*2*pi); zeros(ignore, 1)];
    %TR = TR(2:2:end); % Half the size of the training symbol
    TR_DC = [0; TR];
    TR_full = [TR_DC; flip(conj(TR))];
    tr = sqrt(length(TR_full))*ifft(TR_full);
    tr_prepend = [tr(end-prepend+1:end); tr];

    % Create the train of packets
    x = [];
    TS = 0;
    for i = 1:packets
        dataModG = zeros(freq_hi - freq_lo + 1, 1); % QAM output (freq domain)
        currbit_idx = 1; % Keep track of which bits we are converting
        for j = 1:length(qam_bps)
            % The number of bits we are placing in the current QAM symbol
            bitsPerSymbol = qam_bps(j);
            
            % The range of freq indices that uses bitsPerSymbol
            f_lo = qam_bits_idx(j) + 1;
            f_hi = qam_bits_idx(j + 1);
            
            % Calculate the normalization factor for QAM power limit
            refconst = qammod(0:bitsPerSymbol-1,2^bitsPerSymbol);
            nf = modnorm(refconst,'avpow',P);
            
            % The number of bits we are processing this iteration
            numbits_now = (f_hi - f_lo + 1)*bitsPerSymbol;
            
            % Extract the bits that we are converting
            dataSymbolBits = bits(((i-1)*BPP + currbit_idx):...
                ((i-1)*BPP + currbit_idx + numbits_now - 1));
            
            % Reshape the bits and convert to integers for qammod
            dataBitsMatrix = reshape(dataSymbolBits,...
                length(dataSymbolBits)/bitsPerSymbol,bitsPerSymbol);
            dataIntegers = bi2de(dataBitsMatrix);
            
            % Perform QAM using Grey coding
            qamSymbol = nf*qammod(dataIntegers,2^bitsPerSymbol); % Gray coding
            
            % Add to the QAM output vector
            dataModG(f_lo:f_hi) = qamSymbol;
            
            % For debugging
            if i == 3 && j==1
                qamout = dataIntegers;
            end
            
            % Update current bit index
            currbit_idx = currbit_idx + numbits_now;
        end
        
        % For debugging
%         if i == 3
%             qamout = dataModG;
%         end

        % Append 0s in the freq domain so we don't use freq above 18khz
        % Also apply random phases
        dataModG = [zeros(freq_lo-1,1); dataModG.*exp(1i*randphase*2*pi); zeros(ignore, 1)];

        % Pad 0 for DC and append 0 if length of dataModG is even
        dataModG_DC = [0; dataModG];

        % Append flipped conjugate so that time domain is purely real
        dataModG_DC_full = [dataModG_DC; flip(conj(dataModG))];

        % iDFT to get time domain signal
        dataModG_td = sqrt(length(dataModG_DC_full))*ifft(dataModG_DC_full);

        % Prepend
        dataModG_td_prep = [dataModG_td(end-prepend+1:end); dataModG_td];

        % Whether or not to throw in a training symbol
        if mod(i-1, DPPTP) == 0
            dataModG_td_full = [tr_prepend; dataModG_td_prep];
            TS = TS + 1;
        else
            dataModG_td_full = [dataModG_td_prep];
        end

        % Append to final time domain vector
        x = [x; dataModG_td_full];
    end
    %We need to create a wav file from x. Spec'd by project.
    %audiowrite('tx.wav', x, 44100, 'BitsPerSample', 24);
end

