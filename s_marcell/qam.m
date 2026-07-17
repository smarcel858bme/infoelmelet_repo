%%

clear

% ============ PARAMÉTEREK ============
trellis = poly2trellis(7,[147 115],147);
numframes = 1;
K = 20;

n_sym = 400;
n     = 3*n_sym;        % 1200 adatbit, kódszó = 2400 bit = 800 durva + 1600 finom



for f0 = 0.1:0.1:2
    for c0 = 0.1:0.1:2     
        nrm     = sqrt(c0^2 + f0^2);
        cScale  = c0/nrm;
        fScale  = f0/nrm;
        
        % ============ KONSTELLÁCIÓK ÉS BIT-LEKÉPEZÉS ============
        coarseConst = qammod((0:3).',  4, 'UnitAveragePower', true);
        fineConst   = qammod((0:15).',16, 'UnitAveragePower', true);
        
        cBitsTab = zeros(4,2);
        for s = 0:3
            b = qamdemod(qammod(s,4,'UnitAveragePower',true), 4, ...
                         'UnitAveragePower',true, 'OutputType','bit');
            cBitsTab(s+1,:) = b(:).';
        end
        fBitsTab = zeros(16,4);
        for s = 0:15
            b = qamdemod(qammod(s,16,'UnitAveragePower',true), 16, ...
                         'UnitAveragePower',true, 'OutputType','bit');
            fBitsTab(s+1,:) = b(:).';
        end
        
        combPts = cScale*coarseConst + fScale*fineConst.';   % 4x16
        
        cSymFromBits = zeros(4,1);   cSymFromBits(bi2de(cBitsTab,'left-msb')+1) = (0:3).';
        fSymFromBits = zeros(16,1);  fSymFromBits(bi2de(fBitsTab,'left-msb')+1) = (0:15).';
        
        SNR_vector = 1:2:15;
        numTrials  = 100;
        
        Err_2 = zeros(1, numel(SNR_vector));
        place = zeros(1, K);
        
        for si = 1:numel(SNR_vector)
            snr = SNR_vector(si);
            err = zeros(1, numTrials);
        
            for k = 1:numTrials
                ErrCount = zeros(1,K);
        
                % ---- 1. Adatgenerálás és kódolás ----
                data      = cell(1,K);
                coarseSym = cell(1,K);
                fineSym   = cell(1,K);
                for j = 1:K
                    data{j} = randi([0 1], n, numframes, 'int8');
                    cd = Code(data{j}, trellis);
        
                    cb = reshape(cd(1:2*n_sym),        [2 n_sym])';
                    fb = reshape(cd(2*n_sym+1:6*n_sym),[4 n_sym])';
        
                    coarseSym{j} = cSymFromBits(bi2de(double(cb),'left-msb')+1).';
                    fineSym{j}   = fSymFromBits(bi2de(double(fb),'left-msb')+1).';
                end
        
                % ---- 2. Aszinkron moduláció ----
                tx = zeros(1, n_sym*(K+1));
                tx(1:n_sym) = qammod(coarseSym{1}, 4, 'UnitAveragePower', true);
                for j = 1:K-1
                    tx(n_sym*j+1 : n_sym*(j+1)) = ...
                        cScale*qammod(coarseSym{j+1}, 4,'UnitAveragePower',true) + ... % 
                        fScale*qammod(fineSym{j},    16,'UnitAveragePower',true); %   
                end
                tx(n_sym*K+1 : n_sym*(K+1)) = qammod(fineSym{K},16,'UnitAveragePower',true);
        
                % ---- 3. Csatorna ----
                Pavg = mean(abs(tx).^2);
                N0   = Pavg / 10^(snr/10);
                rx   = tx + sqrt(N0/2)*(randn(size(tx)) + 1i*randn(size(tx)));
        
                % ---- 4. Szukcesszív dekódolás ----
                decodedData = cell(1,K);
                w = rx;
        
                for j = 1:K-1
                    if j == 1
                        llrC = qamdemod(w(1:n_sym), 4, 'UnitAveragePower',true, ...
                                        'OutputType','approxllr','NoiseVariance',N0);
                    else
                        llrC = qamdemod(w(1:n_sym)/cScale, 4, 'UnitAveragePower',true, ...
                                        'OutputType','approxllr','NoiseVariance',N0/cScale^2);
                    end
                    llrC = reshape(llrC, [2 n_sym])';
        
                    llrF = jointFineLLR(w(n_sym+1:2*n_sym), combPts, fBitsTab, N0);
        
                    llrVec = [reshape(llrC.',[],1); reshape(llrF.',[],1)];
        
                    decodedData{j} = Decode(llrVec, trellis);
                    ErrCount(j) = biterr(data{j}, decodedData{j});
        
                    re  = Code(decodedData{j}, trellis);
                    rfb = reshape(re(2*n_sym+1:6*n_sym),[4 n_sym])';
                    reFineSym = fSymFromBits(bi2de(double(rfb),'left-msb')+1).';
        
                    w(n_sym+1:2*n_sym) = w(n_sym+1:2*n_sym) - ...
                        fScale*qammod(reFineSym,16,'UnitAveragePower',true);
        
                    w = w(n_sym+1:end);
                end
        
                % ---- Az utolsó (K-adik) kódszó ----
                llrC = qamdemod(w(1:n_sym)/cScale, 4, 'UnitAveragePower',true, ...
                                'OutputType','approxllr','NoiseVariance',N0/cScale^2);
                llrC = reshape(llrC, [2 n_sym])';
                llrF = qamdemod(w(n_sym+1:2*n_sym), 16, 'UnitAveragePower',true, ...
                                'OutputType','approxllr','NoiseVariance',N0);
                llrF = reshape(llrF, [4 n_sym])';
        
                llrVec = [reshape(llrC.',[],1); reshape(llrF.',[],1)];
                decodedData{K} = Decode(llrVec, trellis);
                ErrCount(K) = biterr(data{K}, decodedData{K});
        
                err(k) = sum(ErrCount)/K;
                if any(ErrCount ~= 0)
                    [~, idx] = max(ErrCount);
                    place(idx) = place(idx) + 1;
                end
            end
        
            Err_2(si) = mean(err)/n;
            disp(['SNR: ', num2str(snr), ' dB | BER: ', num2str(Err_2(si))]);
        end
        
        % ============ DIAGNOSZTIKA KIÍRÁS ============
        figure;
        semilogy(SNR_vector, Err_2, 'b.-', 'LineWidth', 1.5, 'MarkerSize', 12);
        grid on;
        xlabel(['Jel-zaj viszony - SNR (dB)', 'c0: ', num2str(c0),'f0: ', num2str(f0)]);
        ylabel('Bit hibaarány - BER');
        title('4-64 QAM aszinkron moduláció - szukcesszív dekódolás');
        legend('Szimulált 4-64 QAM (aszinkron)');
        
        figure;
        bar(place);
        xlabel('Kódszó sorszáma');
        ylabel('Hányszor volt ez a legrosszabb kódszó');
        title('Hibaterjedés a szukcesszív láncban');
        
     end
end
        % ============ SEGÉDFÜGGVÉNYEK ============
        
        function output = Code(InfBits, Params)
            output = convenc(InfBits, Params);
            output = output(:);
        end
        
        function EstimatedInfBits = Decode(LLR, Params)
            tbdepth = 60;
            % A) teszt alapján: vitdec 'unquant' konvenciója pozitív -> 0 bit,
            % ugyanaz mint a qamdemod LLR-é. NINCS negálás.
            EstimatedInfBits = vitdec(LLR, Params, tbdepth, 'trunc', 'unquant');
            EstimatedInfBits = EstimatedInfBits(:);
        end
        
        function llrFine = jointFineLLR(rxBlock, combPts, fBitsTab, N0)
            n_sym = numel(rxBlock);
            r = reshape(rxBlock, n_sym, 1, 1);
            P = reshape(combPts, 1, 4, 16);
        
            d2    = abs(r - P).^2;
            d2min = squeeze(min(d2, [], 2));   % n_sym x 16
        
            llrFine = zeros(n_sym, 4);
            for b = 1:4
                m0 = min(d2min(:, fBitsTab(:,b)==0), [], 2);
                m1 = min(d2min(:, fBitsTab(:,b)==1), [], 2);
                llrFine(:,b) = (m1 - m0) / N0;
            end
        end

%%
S = load('sync_ref.mat');

figure;
semilogy(SNR_vector, Err_2, 'b.-', 'LineWidth', 1.5, 'MarkerSize', 12);
hold on;
semilogy(S.SNR_vector, S.Err_sync, 'r.--', 'LineWidth', 1.5, 'MarkerSize', 12);
hold off;
grid on;
xlabel('SNR (dB)');
ylabel('Bit hibaarány - BER');
title('aszinkron 4-64 QAM vs. sszinkron 64-QAM');
legend('Aszinkron (4-64 QAM, szukcesszív)', 'Szinkron (64-QAM)', 'Location', 'southwest');
ylim([1e-5 1]);


%%
i = 0;
for s = 401:800:7600
    figure;
    plot(tx(s:s+800), 'or');
    title(num2str(i));
    i = i + 1;
end