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
    freq_hi = 4000;
    
    % Figure out how many bits we are sending per QAM symbol
    load('QAMbits.mat');
    
    % Find where the number of bits per QAM symbol changes. This is to
    % vectorize the qammod code to speed up performance. 
    qam_bits = qam_bits(freq_lo:freq_hi); % Truncate to used freq range
    qam_bits(qam_bits > 6) = 6; % Clip the qam bits
    qam_bits_idx = [0 find(diff(qam_bits)) length(qam_bits)];
    qam_bps = qam_bits(qam_bits_idx(2:end));
    
    % Length of the encoded signal (hardcoded)
    %enclen = 150015;

    % The number of bits per symbol
    BPP = sum(qam_bits);
    
    % Number of data packets per training packet and # of training symbols
    SPTP = 2;
    TS = ceil(numbits/BPP/SPTP);
    
    % The number of zeros to append in freq domain for a cut off of 18kHz
    ignore = ceil(freq_hi/18000*(22050-18000));
    
    % Number of samples to prepend
    prepend = 200;
    
    %Power constraint
    P = 0.00125/((BPP + prepend/2)/(BPP+ignore + prepend/2));
    P = 0.00125;
    
    % Calculate the # of packets, amount of padding in last packet, and
    % the number of samples per packet
    SPP = (freq_hi + ignore)*2 + 1 + prepend;
    
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
            if abs(y(start_idx)) > 0.0003
                found = 1;
            else
                start_idx = start_idx + 1;
            end
        end
        y = y(start_idx:start_idx + enclen);
    end
    
    % Decode each symbol
    dataBitsOutFull = [];
    for t = 1:TS
        if t*SPP*(SPTP + 1) > length(y)
            symsty = y(((t-1)*SPP*(SPTP + 1)+1):end);
            SPTP_adj = length(symsty)/SPP - 1;
        else
            symsty = y(((t-1)*SPP*(SPTP + 1)+1):t*SPP*(SPTP + 1));
            SPTP_adj = SPTP;
        end
        
        % Extract the training and data samples
        tr = symsty(1:SPP);
        symsy = symsty(SPP+1:end);

        % Remove prepends
        tr = tr(prepend + 1:end);

        % Decode the training symbols
        TR = 1/sqrt(length(tr))*fft(tr);
        TR = TR(freq_lo+1:freq_hi+1);
        lambda = TR./(sqrt(P)*exp(1i*randphase*2*pi));
        
        for i = 1:SPTP_adj
            % Extract one symbol
            symy = symsy(((i-1)*SPP+1):i*SPP);

            % Remove prepends
            symy = symy(prepend + 1:end);

            % Take DFT of the recieved signal
            Y = 1/sqrt(length(symy))*fft(symy);

            % Remove bottom half
            Y = Y(1:ceil(end/2));

            % Drop DC component
            Y = Y(2:end);

            % Extract the used freq range
            Y = Y(freq_lo:freq_hi);

            % Remove channel effects
            Y = Y./lambda;

            % Remove the random phase
            Y = Y./exp(1i*randphase*2*pi);

            % Decode QAM
            boundary = sqrt(P/2);
            currbit_idx = 1; % Keep track of which bits we are converting
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
                dataSymbolsOut = Y(f_lo:f_hi)/nf;
                %scatterplot(dataSymbolsOut,1,0,'g.')
                
                % QAM demod the symbols
                dataIntegersOut = qamdemod(dataSymbolsOut,2^bitsPerSymbol, -1.63*((i-1)*length(qam_bps)+j)/180*pi);

                % Convert to bits
                dataBitsOut = de2bi(dataIntegersOut,bitsPerSymbol);
                
                % The number of bits we are processing this iteration
                numbits_now = (f_hi - f_lo + 1)*bitsPerSymbol;
            
                if (length(dataBitsOutFull) + length(dataBitsOut(:))) <= numbits
                    correctBits = refbits(length(dataBitsOutFull)+1:...
                        (length(dataBitsOutFull) + length(dataBitsOut(:))));
                    [all(dataBitsOut(:) == correctBits) bitsPerSymbol]
                end
                
                % Append to output bits
                dataBitsOutFull = [dataBitsOutFull; dataBitsOut(:)];
                
                % Update current bit index
                currbit_idx = currbit_idx + numbits_now + 1;
            end
            
            %scatterplot(Y(1001:2000),1,0,'g.');
        end
    end
    extracted = 0;
    outbits = dataBitsOutFull;                 % Return data in column vector
    outbits = outbits(1:numbits);
end

