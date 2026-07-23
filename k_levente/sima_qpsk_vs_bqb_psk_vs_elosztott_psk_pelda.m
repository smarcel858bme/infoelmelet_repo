%% Vegyes Modulációs (BPSK-QPSK) Konvolúciós Kódolás Szimuláció
% Ez a szkript bemutatja az "Egyenlőtlen Hibavédelem" (UEP) módszereit.
% 3 rendszert versenyeztet:
% 1. Tiszta QPSK
% 2. BPSK a széleken (Anchor effektus)
% 3. Elosztott BPSK szegmensek (Periodikus horgonyok a hibaterjedés ellen)

clear; clc; close all;

%% 1. Rendszer paraméterek
K = 7; % Kényszerhossz
poly = [171 133]; % Ipari standard polinomok
trellis = poly2trellis(K, poly);

dataLen = 994;          % Hasznos adatbitek száma
tailLen = K - 1;        % Tail (farok) bitek a regiszter kiürítéséhez
% R=1/2 miatt a kódolt bitszám pontosan a duplája:
codedLen = (dataLen + tailLen) * 2; % 1000 * 2 = 2000 kódolt bit

snr_dB_array = 0:0.5:4.5; % SNR tartomány
max_keret = 1000;         % Szimulált keretek SNR pontonként

% Eredménytárolók
ber_qpsk  = zeros(size(snr_dB_array));
ber_edge  = zeros(size(snr_dB_array));
ber_distr = zeros(size(snr_dB_array));

fprintf('--- Szimuláció indítása: Különböző UEP elrendezések ---\n');
fprintf('Kódolt kerethossz: %d bit (Mindenhol 400 BPSK és 1600 QPSK bittel gazdálkodunk!)\n\n', codedLen);

for snr_idx = 1:length(snr_dB_array)
snr = snr_dB_array(snr_idx);
hiba_qpsk = 0; hiba_edge = 0; hiba_distr = 0;

nv = 10^(-snr/10); % Zajteljesítmény

for k = 1:max_keret
    %% --- ADATGENERÁLÁS ÉS KÓDOLÁS (Közös mind a 3 rendszerhez) ---
    data = randi([0 1], dataLen, 1);
    dataWithTail = [data; zeros(tailLen, 1)]; 
    codedData = convenc(dataWithTail, trellis); % 2000 bit
    
    %% --- 1. RENDSZER: TISZTA QPSK (REFERENCIA) ---
    cD_qpsk = reshape(codedData, 2, [])';
    tx_qpsk = pskmod(bi2de(cD_qpsk, "left-msb"), 4, pi/4);
    rx_qpsk = awgn(tx_qpsk, snr, 'measured');
    llr_qpsk_pairs = pskdemod(rx_qpsk, 4, pi/4, 'OutputType', 'llr', 'NoiseVariance', nv);
    
    dec_qpsk_full = vitdec(llr_qpsk_pairs(:), trellis, 100, 'term', 'unquant');
    hiba_qpsk = hiba_qpsk + biterr(data, dec_qpsk_full(1:dataLen));
    
    %% --- 2. RENDSZER: SZÉLEKEN BPSK ---
    % Felosztás: 200 BPSK | 1600 QPSK | 200 BPSK
    tx_edge = [ pskmod(codedData(1:200), 2, 0); ...
                pskmod(bi2de(reshape(codedData(201:1800), 2, [])', "left-msb"), 4, pi/4); ...
                pskmod(codedData(1801:2000), 2, 0) ];
            
    rx_edge = awgn(tx_edge, snr, 'measured');
    
    % Demoduláció szegmensenként
    llr_edge_1 = pskdemod(rx_edge(1:200), 2, 0, 'OutputType', 'llr', 'NoiseVariance', nv);
    llr_edge_2 = pskdemod(rx_edge(201:1000), 4, pi/4, 'OutputType', 'llr', 'NoiseVariance', nv); % 800 szimbólum = 1600 bit
    llr_edge_3 = pskdemod(rx_edge(1001:1200), 2, 0, 'OutputType', 'llr', 'NoiseVariance', nv);
    
    llr_edge_stream = [llr_edge_1; llr_edge_2(:); llr_edge_3];
    dec_edge_full = vitdec(llr_edge_stream, trellis, 100, 'term', 'unquant');
    hiba_edge = hiba_edge + biterr(data, dec_edge_full(1:dataLen));
    
    %% --- 3. RENDSZER: ELOSZTOTT BPSK (A továbbfejlesztett ötlet) ---
    % Ugyanannyi (400) BPSK bitet használunk, de 9 blokkba osztva:
    % 5x 80 BPSK, és köztük 4x 400 QPSK (összesen 1200 szimbólum)
    
    tx_distr = zeros(1200, 1); % Előfoglalás a sebességért
    bit_idx = 1;
    sym_idx = 1;
    
    % Modulációs Ciklus
    for b = 1:9
        if mod(b, 2) ~= 0  % Páratlan blokk: BPSK (80 bit)
            L_bit = 80;
            L_sym = 80;
            tx_distr(sym_idx : sym_idx+L_sym-1) = pskmod(codedData(bit_idx : bit_idx+L_bit-1), 2, 0);
        else               % Páros blokk: QPSK (400 bit -> 200 szimbólum)
            L_bit = 400;
            L_sym = 200;
            qpsk_in = bi2de(reshape(codedData(bit_idx : bit_idx+L_bit-1), 2, [])', "left-msb");
            tx_distr(sym_idx : sym_idx+L_sym-1) = pskmod(qpsk_in, 4, pi/4);
        end
        bit_idx = bit_idx + L_bit;
        sym_idx = sym_idx + L_sym;
    end
    
    % Csatorna zaj
    rx_distr = awgn(tx_distr, snr, 'measured');
    
    % Demodulációs Ciklus
    llr_distr_stream = zeros(2000, 1); % 2000 kódolt bitnyi LLR fog visszajönni
    sym_idx = 1;
    llr_idx = 1;
    
    for b = 1:9
        if mod(b, 2) ~= 0  % Páratlan blokk: BPSK demoduláció
            L_sym = 80;
            L_llr = 80;
            ll = pskdemod(rx_distr(sym_idx : sym_idx+L_sym-1), 2, 0, 'OutputType', 'llr', 'NoiseVariance', nv);
        else               % Páros blokk: QPSK demoduláció
            L_sym = 200;
            L_llr = 400;
            ll = pskdemod(rx_distr(sym_idx : sym_idx+L_sym-1), 4, pi/4, 'OutputType', 'llr', 'NoiseVariance', nv);
        end
        llr_distr_stream(llr_idx : llr_idx+L_llr-1) = ll(:);
        sym_idx = sym_idx + L_sym;
        llr_idx = llr_idx + L_llr;
    end
    
    dec_distr_full = vitdec(llr_distr_stream, trellis, 100, 'term', 'unquant');
    hiba_distr = hiba_distr + biterr(data, dec_distr_full(1:dataLen));
end

% BER számítása
ber_qpsk(snr_idx)  = hiba_qpsk / (max_keret * dataLen);
ber_edge(snr_idx)  = hiba_edge / (max_keret * dataLen);
ber_distr(snr_idx) = hiba_distr / (max_keret * dataLen);

fprintf('SNR: %2.1f dB | Tiszta: %e | Szélek: %e | Elosztott: %e\n', ...
        snr, ber_qpsk(snr_idx), ber_edge(snr_idx), ber_distr(snr_idx));


end

%% --- Eredmények kirajzolása ---
figure('Name', 'UEP Stratégiák Versenye', 'NumberTitle', 'off');
semilogy(snr_dB_array, ber_qpsk, 'r-x', 'LineWidth', 2, 'MarkerSize', 8); hold on;
semilogy(snr_dB_array, ber_edge, 'b-o', 'LineWidth', 2, 'MarkerSize', 8);
semilogy(snr_dB_array, ber_distr, 'g-^', 'LineWidth', 2, 'MarkerSize', 8);
grid on;
title('Waterfall Curve: UEP Elrendezések Hatása');
xlabel('Jel-Zaj viszony (SNR) [dB]');
ylabel('Bithibaarány (BER)');
legend('Tiszta QPSK', 'Széleken BPSK (200-1600-200)', 'Elosztott BPSK (5x80 BPSK)', 'Location', 'southwest');