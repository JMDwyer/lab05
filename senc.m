function [] = senc(input_bits)

rand_inst = load('rand_inst.mat','rand_inst');
rand_realizations = rand_inst.rand_inst;
outbits = [];

P = 0.00125;
numbits = length(input_bits);

for i = 1:50
 
    %break 4000 bits off of input bits
    data_bits = input_bits((4000*(i-1)+1:4000*(i))); 
    
    %Multiply each element within the length 4000 vector of data bits by
    %root2
    data = sqrt(2).*data_bits;
    bit_tone = data.*(rand_realizations);
    
    %Pad zeros onto bit tones
    padded_tone = [zeros(112,1);bit_tone;zeros(708,1)];

    %Flip bit_tone conjugate,prepend dc offset
    freq_symbol = [0;padded_tone;flip(conj(padded_tone))];

    %useful variables
    L = length(freq_symbol);
    n_plus = 120;
    
    %Take the IDFT
    time_symbol = sqrt(L).*ifft(freq_symbol);

    %add a cyclic prefix, multiply by sqrt L
    cyclic_time_symbol = [time_symbol((L-n_plus+1):end);time_symbol]; 

    outbits = [outbits;cyclic_time_symbol];

end

%Build the training symbol
training_data = [zeros(112,1);rand_realizations;zeros(708,1)];
training_symbol = [0;training_data;flip(conj(training_data))];

L = length(training_symbol);
n_plus = 120;

training_Ifft = sqrt(L)*ifft(training_symbol);

%Cyclic prepend
training_cyclic = [training_Ifft((L-n_plus+1):end);training_Ifft];

%Add the training symbol to the remaining bits
y = [training_cyclic; outbits];

%Write to waveform
audiowrite('tx.wav', y, 44100, 'BitsPerSample', 24);

return
