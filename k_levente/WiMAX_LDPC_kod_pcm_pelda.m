%% RPTU Paritásellenőrző Mátrix ($H$) beolvasása és Szimulációja
% Ez a szkript bemutatja, hogyan dolgozzuk fel az RPTU adatbázisból származó
% kicsomagolt AList formátumú fájlokat, és hogyan szimuláljuk le rajtuk a csatornakódolást.

clear; clc; close all;

hFilename = 'wimax_576_0.5.alist';

%% 2. A fájl beolvasása és Sparse Mátrix generálása (A parser)
fprintf('--- Fájl feldolgozása ---\n');
H_sparse = readAListToSparse(hFilename);

% Méretek meghatározása
[M, N] = size(H_sparse); 
K = N - M; % Hasznos bitek száma (közelítés, ha H teljes rangú)
kodarany = K / N;

fprintf('Mátrix mérete: %d oszlop (N), %d sor (M)\n', N, M);
fprintf('Kódarány: ~%.2f\n\n', kodarany);

% Opcionális: a sparse mátrix mentése MATLAB natív formátumba a jövőre nézve
save('rptu_matrix.mat', 'H_sparse'); 

%% 3. Kommunikációs Szimuláció beállítása
% Használjuk a "Csupa-Nulla" trükköt, hogy ne kelljen a Generátor mátrixszal (G) 
% vesződni. Csupa 0-t küldünk.

snr_dB_array = 0:2:8; % Szimulált Jel-Zaj viszonyok (SNR) decibelben
max_hibas_blokk = 50; % Hány hibás keretig fusson egy SNR ponton
max_keret = 10000;    % Maximum hány keretet küldjünk
ber_results = zeros(size(snr_dB_array));

% Dekódoló beállítása: A MATLAB beépített üzenettovábbító (Belief Propagation) 
% LDPC dekódoló konfigurációját használjuk a mátrixunkhoz az új ldpcDecode függvénnyel.
cfgLDPCDec = ldpcDecoderConfig(H_sparse);
max_iter = 20; % Maximum iterációk száma a dekódoláshoz

fprintf('--- Szimuláció indítása (All-Zero Codeword trükkel) ---\n');
for snr_idx = 1:length(snr_dB_array)
    snr = snr_dB_array(snr_idx);
    hiba_szam = 0;
    bit_szam = 0;
    hibas_blokkok = 0;
    
    for keret = 1:max_keret
        % 1. Adat "generálása" és Moduláció (Csupa 0 bit -> Csupa +1 szimbólum)
        % BPSK: 0 bit -> +1 feszültség, 1 bit -> -1 feszültség
        tx_szimbolumok = ones(N, 1); 
        
        % 2. Csatorna zaj (AWGN)
        % A zajos jelünk (LLR - Log Likelihood Ratio alapja)
        rx_szimbolumok = awgn(tx_szimbolumok, snr, 'measured');
        
        % Vételi valószínűségek konvertálása LLR-ré (leegyszerűsített képlet)
        % Mivel BPSK-t használunk, a vett jel önmaga használható puha döntésre, 
        % de a komparátor a 0-t szereti, ezért invertáljuk az előjelet a dekódolónak.
        rx_llr = -rx_szimbolumok; 
        
        % 3. Dekódolás a paritásmátrix (H) konfigurációja alapján
        rx_bitek = ldpcDecode(rx_llr, cfgLDPCDec, max_iter);
        
        % 4. Hibaszámítás
        % Mivel csupa 0-t küldtünk, minden 1-es bit, amit a dekódoló kiad, HIBA.
        hibak_a_keretben = sum(rx_bitek);
        hiba_szam = hiba_szam + hibak_a_keretben;
        bit_szam = bit_szam + K; % A dekódoló K bitet ad vissza
        
        if hibak_a_keretben > 0
            hibas_blokkok = hibas_blokkok + 1;
        end
        
        if hibas_blokkok >= max_hibas_blokk
            break; % Elég statisztika gyűlt össze ehhez az SNR-hez
        end
    end
    
    ber_results(snr_idx) = hiba_szam / bit_szam;
    fprintf('SNR = %d dB | BER = %e | Értékelt keretek: %d\n', snr, ber_results(snr_idx), keret);
end

%% 4. Eredmények kirajzolása (Vízesés görbe / Waterfall curve)
figure('Name', 'RPTU Kód Teljesítménye', 'NumberTitle', 'off');
semilogy(snr_dB_array, ber_results, 'b-o', 'LineWidth', 2);
grid on;
xlabel('Jel-Zaj viszony (Eb/No vagy SNR) [dB]');
ylabel('Bithibaarány (BER)');
title('Hibajavító kód teljesítménye (Waterfall Curve)');
% Ha a kód jó, a görbe egy "vízesésként" zuhan lefelé.


%% --- FÜGGVÉNYEK ---

function H = readAListToSparse(filename)
    % Ez a függvény egy szabványos RPTU AList fájlt konvertál
    % MATLAB ritka (sparse) mátrixszá.
    
    fid = fopen(filename, 'r');
    if fid == -1
        error('Nem sikerült megnyitni a fájlt!');
    end
    
    % 1. sor: N (oszlopok) és M (sorok)
    meretek = fscanf(fid, '%d', 2);
    N = meretek(1); M = meretek(2);
    
    % 2. sor: dmax_col, dmax_row (nem használjuk közvetlenül)
    fscanf(fid, '%d', 2);
    
    % 3. sor: Oszlopok súlya (hány 1-es van az egyes oszlopokban)
    col_weights = fscanf(fid, '%d', N);
    
    % 4. sor: Sorok súlya
    fscanf(fid, '%d', M);
    
    % Itt tároljuk a nem-nulla elemek koordinátáit
    row_indices = [];
    col_indices = [];
    
    % 5-től N+4-ig. sorok: Minden oszlophoz felsorolja, melyik sorokban van 1-es
    for c = 1:N
        % Beolvasunk annyi számot, amennyi az oszlop súlya
        r_idx = fscanf(fid, '%d', col_weights(c));
        % Hozzáadjuk a koordinátákat a listához
        row_indices = [row_indices; r_idx];
        col_indices = [col_indices; repmat(c, col_weights(c), 1)];
    end
    
    fclose(fid);
    
    % MATLAB beépített sparse konstruktora összerakja a ritka mátrixot a 
    % sor és oszlop koordinátákból (az értékük mindenhol 1).
    % Ensure indices are positive integers; convert from 0-based to 1-based if needed
    if isempty(row_indices) || isempty(col_indices)
        error('No nonzero entries found in file.');
    end
    if any(row_indices < 0) || any(col_indices < 0)
        error('Negative indices found in parsed AList data.');
    end
    % Robust zero-based detection: only convert if there are zeros and max matches zero-based range
    if any(row_indices == 0) || any(col_indices == 0)
        if max(row_indices) <= M-1 && max(col_indices) <= N-1
            row_indices = row_indices + 1;
            col_indices = col_indices + 1;
        else
            error('Parsed indices look like zeros but their max exceeds zero-based bounds.');
        end
    end

    % Validate final bounds
    if max(row_indices) > M || min(row_indices) < 1 || max(col_indices) > N || min(col_indices) < 1
        error('Parsed indices exceed declared matrix dimensions M=%d, N=%d.', M, N);
    end

    H = sparse(row_indices, col_indices, 1, M, N);


end