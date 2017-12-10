function [outbits] = sdec()

numbits = 200000;
n_plus = 120;

%Get Waveform from wav
y = audioread('rx.wav');

%Get Untainted rand_realizations to compare against training
rand_inst = load('rand_inst.mat','rand_inst');
rand_realizations = rand_inst.rand_inst;

enclen = 497811;
%Detect delays, remove zeroes
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

%Peel off training symbol
training_cyclic = y(1:9761);
%remove prepend 
tr_time = training_cyclic(n_plus+1:end);

%Take fft
tr_freq = 1/sqrt(length(tr_time))*fft(tr_time);
   
%remove bottom half
tr_freq = tr_freq(1:ceil(end/2));

%Drop DC component
tr_freq = tr_freq(2:end);
 
%prune off zeroes for high and low frequencies
tr_freq = tr_freq(113:4112);

lambda = tr_freq./rand_realizations;

y = y(9762:end);

outbits = zeros(numbits, 1);
for i =1:50
    %Chunk into samples
    sample_i = y(((i-1)*9761+1):i*(9761));
    
    %remove prepends
    freq_data = sample_i(n_plus+1:end);
    
    %DFT
    Y = (1/sqrt(length(freq_data)))*fft(freq_data);
   
    %remove bottom half
    Y = Y(1:ceil(end/2));

    %Drop DC component
    Y = Y(2:end);
    
    %prune off zeroes for high and low frequencies
    Y = Y(113:4112);
    
    % Remove channel effects
    Y = Y./lambda;

    % Remove the random phase
    Y = Y./rand_realizations;
    
    % Decode OOK
    boundary = sqrt(1/2);
    for j = (1:4000)
                if((abs(Y(j)) >= boundary))
                    outbits(4000*(i-1)+j) = 1;
                else
                    outbits(4000*(i-1)+j) = 0;
                end
     end
end
        
end    

