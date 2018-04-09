function [signal_pass,fsfd, fsymb, span, fsamp, Tsamp]=transmitter(packet,fc)

% Barker sequence for synchronization
barkerSeq = [1, 1, 1, -1, -1, -1, 1, -1, -1, 1, -1]; 
% Generate symbol
fsamp = 44e3;                                       % sampling frequency [Hz]
tsamp = 1/fsamp;                                    % Sampling time
rb = 440;                                           % bit rate [bit/sec]

% Constellation or bit to symbol mapping
s = [(1 + 1i) (1 - 1i) (-1 -1i) (-1 + 1i)]/sqrt(2); % Constellation 1 - QPSK/4-QAM
M = length(s);                                      % Number of symbols in the constellation
m = log2(M);                                        % Number of bits per symbol
fsymb = rb/m;                                       % Symbol rate [symb/s]
fsfd = fsamp/fsymb;                                 % Number of samples per symbol (choose fs such that fsfd is an integer for simplicity) [samples/symb]

% Map bits to symbols
message = buffer(packet, m)';                       % Message p_buffer
sym_idx = bi2de(message, 'left-msb')'+1;            % Bits to symbol index
symbol = s(sym_idx);                   % Look up symbols using the indices

% Upsample
barker_upsample = upsample(barkerSeq, fsfd); % This is for the signal detection
symbol_upsample = upsample([barkerSeq, symbol], fsfd); % Space the symbols fsfd apart, to enable pulse shaping using conv.

% Generate pulse
span = 6;
[rrcpulse, ~] = rtrcpuls(0.2,1/fsymb,fsamp,span);

% Generate baseband signal and barker sequence
barker_base = conv(rrcpulse, barker_upsample); % Barker sequence baseband
signal_base = conv(rrcpulse,symbol_upsample);  % Baseband signal

% Separate the signal
sI = real(signal_base);
sQ = imag(signal_base);

N1 = length(barker_base);
N2 = length(signal_base);
Bsamp = tsamp*(0:(N1-1));
Tsamp = tsamp*(0:(N2-1));

% Shift the signal to passband
barker_pass = sqrt(2)*barker_base.*cos(2*pi*fc*Bsamp);
signal_passI = sqrt(2)*sI.*cos(2*pi*fc*Tsamp);
signal_passQ = sqrt(2)*sQ.*sin(2*pi*fc*Tsamp);
signal_pass = signal_passI+signal_passQ;            % Passband signal

% Barker sequence on passband is used for signal detection, zeros are
% introduced to allow time for computation
preamble = [barker_pass, zeros(1,N1*2)];

% The final signal is then barker->zeros->barker+signal
final_signal = [preamble, signal_pass];
final_signal = final_signal/max(final_signal);
wavplay(final_signal, fsamp);
end