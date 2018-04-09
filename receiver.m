function [pack, psd, const, eyed] = receiver(tout,fc)

% Barker sequence used for frame- and phase- synchronization
barkerSeq = [1, 1, 1, -1, -1, -1, 1, -1, -1, 1, -1];

c = [1+1i,1-1i,-1-1i,-1+1i]/sqrt(2);                % constellation
fsamp = 44e3;                                       % sampling frequency [Hz]
rb = 440;                                           % bit rate [bit/sec]
M = length(c);                                      % Number of symbols in the constellation
m = log2(M);                                        % Number of bits per symbol
fd = rb/m;                                          % Symbol rate [symb/s]
fsfd = fsamp/fd;                                    % Number of samples per symbol (choose fs such that fsfd is an integer for simplicity) [samples/symb]


% Creation of barker sequence in in-phase waveform at passband
span = 6;
[rrcpulse, t] = rtrcpuls(0.2,1/fd,fsamp,span);
barker_upsample = upsample(barkerSeq, fsfd);
barker_base = conv(rrcpulse, barker_upsample); 
N1 = length(barker_base);

% Butterworth filter used for filtering of signal
[b,a]=butter(6, fd*3/fsamp);

detectFlag = 0;                                     % Flag for detection of preamble
tic                                                 % Initialize tic to compare with tout

% Signal detection phase
while (detectFlag == 0) && (tout > toc)
    sigDet = wavrecord(N1*2,fsamp); % Record a short sample of sound
    N2 = length(sigDet);
    Bsamp1 = 1/fsamp*(0:N2-1); % Create time vector
    sigDet_base = filtfilt(b, a, sigDet*sqrt(2).*cos(2*pi*fc*Bsamp1'));
    C1_conv = conv(sigDet_base, fliplr(barker_base));
    C1 = max(abs(C1_conv)); % Find the largest value of the convolution
    
    if find(C1 > 30) % If the convolution finds the barker sequence
        detectFlag = 1;
        tic
    end
end

% Message recording phase
recMessage = wavrecord(fsfd*350,fsamp); % Record the whole message with safety margins
toc
N3 = length(recMessage);
Bsamp2 = 1/fsamp*(0:N3-1);
% Divide the recorded message into real and imaginary parts
recMessage_base1 = recMessage*sqrt(2).*cos(2*pi*fc*Bsamp2');
recMessage_base2 = recMessage*sqrt(2).*sin(2*pi*fc*Bsamp2');

% Shift to base-band
x_real = filtfilt(b, a, recMessage_base1);
x_imag = filtfilt(b, a, recMessage_base2);
x = x_real + 1j*x_imag;

% Find the start of the message by convolution, a barker sequence is
% present just before the message
C2_conv = conv(x, fliplr(barker_base));
[~,C2] = max(abs(C2_conv));

% The peak of the barker is located at the exact end of the barker sequence
% including the transient, so this is subtracted
SoM = C2-span*fsfd;
SoB = SoM-11*fsfd; % Start of barker sequence

% Synchronize the phase by using the barker sequence
phase_sync = x(SoB:fsfd:SoB+10*fsfd);
phase_synced = mean(angle(phase_sync.'.*barkerSeq));
x_synced = x.*exp(-1i*phase_synced);

%Matched filtering
matched_filter = fliplr([rrcpulse, t]);
mf_output = conv(x_synced, matched_filter)./fsfd;
save
% Truncate the signal to message size
SoM2 = SoM+12*fsfd+span*fsfd;   % First symbol value
filtered_rz = mf_output(SoM2:SoM2+215*fsfd);
% eyed = eyediagram(filtered_rz, fsfd);
eyed.r = filtered_rz;
eyed.fsfd = fsfd;

% PSD calculation
psd = {}; % Initialize psd
[p, f] = pwelch(recMessage,500,300,512,fsamp);
p = 10*log10(p/max(p));
psd.p = p;
psd.f = f-fc;

% Downsample
x_hat = downsample(filtered_rz, fsfd);

% Output constellation
const=x_hat;

% ML decoding
output = zeros(4,length(x_hat));
for i=1:4
    output(i,:) = abs(x_hat-c(i));
end
[~,idx] = min(output);

% Map symbols to bits
message = de2bi(idx-1, 'left-msb');
% Output bits
pack = reshape(message', 1,[]);
save
end