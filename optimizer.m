load('IR0.mat')

irfft = fft(impulse);
irfft = irfft(1:floor(length(irfft)/2));

irfft_interp = interp1(linspace(1,150, 150), irfft, linspace(1,150, 4900));

% Before: 1000
irfft_interp = irfft_interp/(abs(irfft_interp(3800)))*32;

qam_bits = zeros(1,length(irfft_interp));
for i = 1:length(qam_bits)
    if log2(abs(irfft_interp(i))) > 0
        qam_bits(i) = floor(log2(abs(irfft_interp(i))));
    else
        qam_bits(i) = 0;
    end
end

qam_bits_ideal = qam_bits;

qam_bits(qam_bits > 6) = 6; % Clip the qam bits
qam_bits(1:600) = 5;
qam_bits(2900:3145) = 5;
qam_bits(1:200) = 4;

freqs = linspace(0, 22050, 4900);
%freqs = linspace(0, 4900, 4900);
plot(freqs,log2(abs(irfft_interp)))
hold on
plot(freqs,qam_bits_ideal)
plot(freqs,qam_bits)


%title('S11 and S22 magnitudes')
legend('Log(2)', 'Quantized ideal', 'Quantized actual')
xlabel('Frequency (Hz)')
ylabel('Number of bits per QAM symbol')
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
set(gca,'fontsize',18)
set(gcf,'color','w');

sum(qam_bits(1:4000))

save('QAMbits', 'qam_bits')