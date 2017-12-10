load('IR0.mat')

irfft = fft(impulse);
irfft = irfft(1:floor(length(irfft)/2));

irfft_interp = interp1(linspace(1,150, 150), irfft, linspace(1,150, 4900));

irfft_interp = irfft_interp/(abs(irfft_interp(1000)))*32;

qam_bits = zeros(1,length(irfft_interp));
for i = 1:length(qam_bits)
    if log2(abs(irfft_interp(i))) > 0
        qam_bits(i) = floor(log2(abs(irfft_interp(i))));
    else
        qam_bits(i) = 0;
    end
end

plot(log2(abs(irfft_interp)))
hold on
plot(qam_bits)
sum(qam_bits(1:4000))

save('QAMbits', 'qam_bits')