%% Vegyes Modulációs (BPSK-QPSK) Konvolúciós Kódolás Szimuláció
% Ez a szkript bemutatja az "Egyenlőtlen Hibavédelem" (UEP) módszereit.
% 3 rendszert versenyeztet RSC (Rekurzív) konvolúciós kódolóval:
% 1. Tiszta QPSK
% 2. BPSK a széleken (Anchor effektus)
% 3. Elosztott BPSK szegmensek (Periodikus horgonyok a hibaterjedés ellen)

clear; clc; close all;

%% 1. Rendszer paraméterek
K = 7; % Kényszerhossz
trellis = poly2trellis(K, [147 115], 147);

dataLen = 1000;         % Hasznos adatbitek száma
codedLen = dataLen * 2; % 2000 kódolt bit

snr_dB_array = 0:0.5:4.5;
max_keret = 1000;

% Eredménytárolók
ber_qpsk  = zeros(size(snr_dB_array));
ber_edge  = zeros(size(snr_dB_array));
ber_distr = zeros(size(snr_dB_array));

fprintf('--- Szimuláció indítása: Különböző UEP elrendezések (RSC Kóddal) ---\n');

for snr_idx = 1:length(snr_dB_array)
snr = snr_dB_array(snr_idx);
hiba_qpsk = 0; hiba_edge = 0; hiba_distr = 0;

nv = 10^(-snr/10); % Zajteljesítmény LLR-hez

for k = 1:max_keret
    %% --- ADATGENERÁLÁS ÉS KÓDOLÁS ---
    data = randi([0 1], dataLen, 1);
    codedData = convenc(data, trellis); % 2000 bit
    
    %% --- 1. RENDSZER: TISZTA QPSK (REFERENCIA) ---
    cD_qpsk = reshape(codedData, 2, [])';
    tx_qpsk = pskmod(bi2de(cD_qpsk, "left-msb"), 4, pi/4);
    rx_qpsk = awgn(tx_qpsk, snr, 'measured');
    llr_qpsk_pairs = pskdemod(rx_qpsk, 4, pi/4, 'OutputType', 'llr', 'NoiseVariance', nv);
    
    dec_qpsk_full = vitdec(llr_qpsk_pairs(:), trellis, 100, 'trunc', 'unquant');
    hiba_qpsk = hiba_qpsk + biterr(data, dec_qpsk_full);
    
    %% --- 2. RENDSZER: SZÉLEKEN BPSK ---
    % Moduláció: 200 BPSK | 1600 QPSK | 200 BPSK
    tx_e1 = pskmod(codedData(1:200), 2, 0);
    tx_e2 = pskmod(bi2de(reshape(codedData(201:1800), 2, [])', "left-msb"), 4, pi/4);
    tx_e3 = pskmod(codedData(1801:2000), 2, 0);
    
    rx_edge = awgn([tx_e1; tx_e2; tx_e3], snr, 'measured');
    
    % Demoduláció
    ll_e1 = pskdemod(rx_edge(1:200), 2, 0, 'OutputType', 'llr', 'NoiseVariance', nv);
    ll_e2 = pskdemod(rx_edge(201:1000), 4, pi/4, 'OutputType', 'llr', 'NoiseVariance', nv); % 800 QPSK szimbólum
    ll_e3 = pskdemod(rx_edge(1001:1200), 2, 0, 'OutputType', 'llr', 'NoiseVariance', nv);
    
    llr_edge_stream = [ll_e1; ll_e2(:); ll_e3];
    dec_edge_full = vitdec(llr_edge_stream, trellis, 100, 'trunc', 'unquant');
    hiba_edge = hiba_edge + biterr(data, dec_edge_full);
    
    %% --- 3. RENDSZER: ELOSZTOTT BPSK (Periodikus horgonyok) ---
    % Felosztás 9 blokkba: 5x 80 BPSK, 4x 400 QPSK (összesen 1200 szimbólum)
    tx_distr = zeros(1200, 1); % Előfoglalás
    bit_idx = 1;
    sym_idx = 1;
    
    % Modulációs Ciklus
    for b = 1:9
        if mod(b, 2) ~= 0  % Páratlan: BPSK (80 bit)
            L_bit = 80;
            L_sym = 80;
            tx_distr(sym_idx : sym_idx+L_sym-1) = pskmod(codedData(bit_idx : bit_idx+L_bit-1), 2, 0);
        else               % Páros: QPSK (400 bit -> 200 szimbólum)
            L_bit = 400;
            L_sym = 200;
            qpsk_in = bi2de(reshape(codedData(bit_idx : bit_idx+L_bit-1), 2, [])', "left-msb");
            tx_distr(sym_idx : sym_idx+L_sym-1) = pskmod(qpsk_in, 4, pi/4);
        end
        bit_idx = bit_idx + L_bit;
        sym_idx = sym_idx + L_sym;
    end
    
    rx_distr = awgn(tx_distr, snr, 'measured');
    
    % Demodulációs Ciklus
    llr_distr_stream = zeros(2000, 1); % Előfoglalás
    sym_idx = 1;
    llr_idx = 1;
    
    for b = 1:9
        if mod(b, 2) ~= 0  % Páratlan: BPSK
            L_sym = 80;
            L_llr = 80;
            ll = pskdemod(rx_distr(sym_idx : sym_idx+L_sym-1), 2, 0, 'OutputType', 'llr', 'NoiseVariance', nv);
        else               % Páros: QPSK
            L_sym = 200;
            L_llr = 400;
            ll = pskdemod(rx_distr(sym_idx : sym_idx+L_sym-1), 4, pi/4, 'OutputType', 'llr', 'NoiseVariance', nv);
        end
        llr_distr_stream(llr_idx : llr_idx+L_llr-1) = ll(:);
        sym_idx = sym_idx + L_sym;
        llr_idx = llr_idx + L_llr;
    end
    
    dec_distr_full = vitdec(llr_distr_stream, trellis, 100, 'trunc', 'unquant');
    hiba_distr = hiba_distr + biterr(data, dec_distr_full);
end

ber_qpsk(snr_idx)  = hiba_qpsk / (max_keret * dataLen);
ber_edge(snr_idx)  = hiba_edge / (max_keret * dataLen);
ber_distr(snr_idx) = hiba_distr / (max_keret * dataLen);

fprintf('SNR: %2.1f dB | Tiszta: %e | Szélek: %e | Elosztott: %e\n', ...
        snr, ber_qpsk(snr_idx), ber_edge(snr_idx), ber_distr(snr_idx));


end

%% --- Eredmények kirajzolása ---
figure('Name', 'UEP Stratégiák Versenye (RSC kódoló)', 'NumberTitle', 'off');
semilogy(snr_dB_array, ber_qpsk, 'r-x', 'LineWidth', 2, 'MarkerSize', 8); hold on;
semilogy(snr_dB_array, ber_edge, 'b-o', 'LineWidth', 2, 'MarkerSize', 8);
semilogy(snr_dB_array, ber_distr, 'g-^', 'LineWidth', 2, 'MarkerSize', 8);
grid on;
title('Waterfall Curve: UEP Elrendezések Hatása (RSC Kóddal)');
xlabel('Jel-Zaj viszony (SNR) [dB]');
ylabel('Bithibaarány (BER)');
legend('Tiszta QPSK', 'Széleken BPSK (200-1600-200)', 'Elosztott BPSK (5x80 BPSK)', 'Location', 'southwest');