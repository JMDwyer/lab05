function [outbits, extracted] = qamdec(y, enclen, debugin, refbits)

    M = 32;                     % Size of signal constellation
    k = log2(M);                % Number of bits per symbol (QAM symbol)
    n = 200000;                 % Number of bits to process
    numSamplesPerSymbol = 1;    % Oversampling factor
    
    % Number of total bits
    numbits = n;
    
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
    
    % Length of the encoded signal (hardcoded)
    %enclen = 150015;

    % The number of bits per symbol
    BPP = sum(diff(qam_bits_idx).*qam_bps);
    
    % Number of data packets per training packet and # of training symbols
    DPPTP = 14;
    TS = ceil(numbits/BPP/DPPTP);
    
    % The number of zeros to append in freq domain for a cut off of 18kHz
    ignore = ceil(freq_hi/18000*(22050-18000));
    ignore = 4900 - freq_hi;
    
    % Number of samples to prepend
    prepend = 200;
    
    %Power constraint
    P = 0.00125/((BPP + prepend/2)/(BPP+ignore + prepend/2));
    P = 0.00125*22050/18000/0.67;
    
    % Calculate the number of samples per data packet
    SPP = (freq_hi + ignore)*2 + 1 + prepend;
    
    % Calculate the number of samples per training packet
    SPTP = (freq_hi + ignore)*2 + 1 + prepend;
    
    % Generate the random phases
    rng(4670);
    randphase = rand([freq_hi-freq_lo+1, 1]);
    
    %[afterchan, ~] = audioread('rx.wav');
    %y = afterchan';
    
    % Find the starting index and truncate the initial zeros
    if length(y) > enclen
        found = 0;
        start_idx = 1;
        while found == 0
            if abs(y(start_idx)) > 0.0005
                found = 1;
            else
                start_idx = start_idx + 1;
            end
        end
        y = y(start_idx:start_idx + enclen);
    end
    
    % Decode each symbol
    dataBitsOutFull = [];
    lambda_stack = []; % Lambda from the previous packet
    for t = 1:TS
        if t*(SPP*DPPTP + SPTP) > length(y)
            symsty = y(((t-1)*(SPP*DPPTP + SPTP)+1):end);
            DPPTP_adj = (length(symsty) - SPTP)/SPP;
        else
            symsty = y(((t-1)*(SPP*DPPTP + SPTP)+1):t*(SPP*DPPTP + SPTP));
            DPPTP_adj = DPPTP;
        end
        
        % Extract the training and data samples
        tr = symsty(1:SPTP);
        symsy = symsty(SPTP+1:end);

        % Remove prepends
        tr = tr(prepend + 1:end);

        % Decode the training symbols
        TR = 1/sqrt(length(tr))*fft(tr);
        % Remove top half and remove DC
        TR = TR(2:ceil(end/2));
        % Extract the frequency of interest
        TR = TR(freq_lo:freq_hi);
        % Compute channel model
        lambda = smooth(TR./(sqrt(P)*exp(1i*randphase*2*pi)), 31);
        % Add to lambda stack
        lambda_stack = [lambda_stack lambda];
        
        
        %figure(1)
        %plot(angle(lambda))
%         subplot(1,2,1)
%         plot(real(lambda))
%         subplot(1,2,2)
%         plot(imag(lambda))
        
        %figure(1)
        %plot(angle(lambda))
        %hold on
        for i = 1:DPPTP_adj
            % Extract one symbol
            symy = symsy(((i-1)*SPP+1):i*SPP);

            % Remove prepends
            symy = symy(prepend + 1:end);

            % Take DFT of the recieved signal
            Y = 1/sqrt(length(symy))*fft(symy);

            % Remove top half
            Y = Y(1:ceil(end/2));

            % Drop DC component
            Y = Y(2:end);

            % Extract the used freq range
            Y = Y(freq_lo:freq_hi);

            % Remove the random phase
            Y = Y./exp(1i*randphase*2*pi);
            
            % Remove channel effects
            Y_nochan = Y./lambda;

            % Decode QAM
            currbit_idx = 1; % Keep track of which bits we are converting
            dataBitsOutPacket = [];
            for j = 1:length(qam_bps)
                % The number of bits we have in the current QAM symbol
                bitsPerSymbol = qam_bps(j);

                % The range of freq indices that uses bitsPerSymbol
                f_lo = qam_bits_idx(j) + 1;
                f_hi = qam_bits_idx(j + 1);
                
                % Calculate the normalization factor for QAM power limit
                refconst = qammod(0:bitsPerSymbol-1,2^bitsPerSymbol);
                nf = modnorm(refconst,'avpow',P);
                
                % Number of symbols we are processing this iteration
                numsym_now = f_hi - f_lo + 1;

                % Extract the symbols that we are converting
                % Also scale by inverse the power
                dataSymbolsOut = Y_nochan(f_lo:f_hi)/nf;
                %scatterplot(dataSymbolsOut,1,0,'g.')
                
                % QAM demod the symbols
                dataIntegersOut = qamdemod(dataSymbolsOut,2^bitsPerSymbol, -1.63*j/180*pi);
                %dataIntegersOut = qamdemod(dataSymbolsOut,2^bitsPerSymbol, -0.1152525*j/180*pi);

                % Convert to bits
                dataBitsOut = de2bi(dataIntegersOut,bitsPerSymbol);
                
                % The number of bits we are processing this iteration
                numbits_now = (f_hi - f_lo + 1)*bitsPerSymbol;
            
                if (length(dataBitsOutFull) + length(dataBitsOutPacket) + length(dataBitsOut(:))) <= numbits
                    correctBits = refbits(length(dataBitsOutFull)+length(dataBitsOutPacket)+1:...
                        (length(dataBitsOutFull) + length(dataBitsOutPacket) + length(dataBitsOut(:))));
                    allcorrect = all(dataBitsOut(:) == correctBits);
                    if allcorrect == 0 && i == 1
                        [j bitsPerSymbol sum(dataBitsOut(:) ~= correctBits)]
                        %scatterplot(dataSymbolsOut,1,0,'g.')
                    end
                end
                
                % Append to current packet bits
                dataBitsOutPacket = [dataBitsOutPacket; dataBitsOut(:)];
                
                % Update current bit index
                currbit_idx = currbit_idx + numbits_now + 1;
            end
            
            % Append to current packet bits
            dataBitsOutFull = [dataBitsOutFull; dataBitsOutPacket];
            
            % Update the training using the current decoded data. This
            % assumes that the decoded data from the packet right after
            % the training packet is 100% accurate.
            % Start by performing the encoder to figure out the sent OFDM
            % signal
            sentSymbols = zeros(freq_hi - freq_lo + 1, 1); % QAM output (freq domain)
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
                dataSymbolBits = dataBitsOutPacket(currbit_idx:currbit_idx+numbits_now-1);

                % Reshape the bits and convert to integers for qammod
                dataBitsMatrix = reshape(dataSymbolBits,...
                    length(dataSymbolBits)/bitsPerSymbol,bitsPerSymbol);
                dataIntegers = bi2de(dataBitsMatrix);

                % Perform QAM using Grey coding
                qamSymbol = nf*qammod(dataIntegers,2^bitsPerSymbol); % Gray coding

                % Add to the QAM output vector
                sentSymbols(f_lo:f_hi) = qamSymbol;

                % Update current bit index
                currbit_idx = currbit_idx + numbits_now;
            end
            
            % Update the channel estimate lambda
            % Use averaging to reduce effects of any wrongly decoded bits
            lambda_stack = [lambda_stack smooth(Y./sentSymbols, 31)];
            lambda = mean(lambda_stack(:,size(lambda_stack,2):size(lambda_stack,2)),2);
            %plot(angle(lambda))
            
            %scatterplot(Y(1001:2000),1,0,'g.');
        end
    end
    extracted = 0;
    outbits = dataBitsOutFull;                 % Return data in column vector
    outbits = outbits(1:numbits);
end

