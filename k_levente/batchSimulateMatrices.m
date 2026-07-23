function batchSimulateMatrices(inputDir, outputDir, snrRange, maxFrames)
% BATCHSIMULATEMATRICES Végigfuttatja a szimulációt egy mappa összes .txt fájlján.
%
% EREDMÉNY: Minden egyes bemeneti .txt fájlhoz egy KÜLÖNÁLLÓ .mat fájlt
% generál a kimeneti mappába, amely tartalmazza a TurboErr, BlockLengthHalf
% és n változókat.

%% Alapértelmezett paraméterek beállítása
if nargin < 1 || isempty(inputDir),  inputDir = pwd; end 
if nargin < 2 || isempty(outputDir), outputDir = fullfile(pwd, 'Results'); end
if nargin < 3 || isempty(snrRange),  snrRange = 1:0.1:10; end
if nargin < 4 || isempty(maxFrames), maxFrames = 1000; end 

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

files = dir(fullfile(inputDir, '*.txt'));

if isempty(files)
    fprintf('Nem találtam .txt fájlt a mappában: %s\n', inputDir);
    return;
end

%% Fő ciklus: Végigmegyünk az összes fájlon egyenként
for i = 1:length(files)
    currentFile = fullfile(files(i).folder, files(i).name);
    [~, baseFileName, ~] = fileparts(files(i).name);
    
    fprintf('--- Feldolgozás alatt (%d/%d): %s ---\n', i, length(files), files(i).name);
    
    try
        H = readTurboCode(currentFile);
        cfgLDPCEnc = ldpcEncoderConfig(H);
        cfgLDPCDec = ldpcDecoderConfig(H);
        
        n = cfgLDPCEnc.NumInformationBits;
        BlockLengthHalf = cfgLDPCDec.BlockLength / 2;
        K = 20; 
        TurboErr = []; % Eredetileg LdpcErr_2
        numframes = 1;
        place = zeros(1,20);
        
        % Innentől ugyanaz
        for snr = snrRange
            error_count = zeros(1, maxFrames);
            for k_idx = 1:maxFrames
                data = cell(1, K);
                codedData = cell(1, K);
                cD1 = cell(1, K);
                
                for j=1:K
                    data{j} = randi([0 1],n,numframes,'int8');
                end
                for j=1:K
                    codedData{j} = ldpcEncode(data{j},cfgLDPCEnc);
                    cD1{j} = reshape(codedData{j}, [2 BlockLengthHalf])';
                end
                
                input = pskmod(cD1{1}(:,1), 2, 0);
                for j=1:K-1
                    input1 = bi2de([cD1{j}(:,2), cD1{j+1}(:,1)], "left-msb");
                    input = [input; pskmod(input1, 4, pi/4)];
                end
                input = [input; pskmod(cD1{K}(:,2), 2, 0)];
                
                output = awgn(input, snr, 'measured');
                
                LLR = cell(1, K);
                decodedData = cell(1, K);
                ErrCount = zeros(1, K);
                
                for j=1:K-1
                    LLR1Half = pskdemod(output(BlockLengthHalf*(j-1)+1:BlockLengthHalf*j),2,0,OutputType="llr");
                    LLR2Half = reshape(pskdemod(output(BlockLengthHalf*j+1:BlockLengthHalf*(j+1)),4,pi/4,OutputType="llr"),[2 BlockLengthHalf])';
                    LLR{j} = [LLR1Half, LLR2Half(:,1)];
                    LlrActual = LLR{j}';
                    
                    decodedData{j} = ldpcDecode(LlrActual(:), cfgLDPCDec, 10);
                    
                    EstimatedCodeword = ldpcEncode(decodedData{j}, cfgLDPCEnc);
                    EstimatedCodeword = reshape(EstimatedCodeword, [2 BlockLengthHalf])';
                    SecondPartECWinFourier = pskmod(EstimatedCodeword(:,2), 2, pi/2);
                    output(BlockLengthHalf*j+1:BlockLengthHalf*(j+1)) = output(BlockLengthHalf*j+1:BlockLengthHalf*(j+1)) - 1/sqrt(2)*SecondPartECWinFourier;    
                end
                
                LLR{K} = reshape(pskdemod(output(BlockLengthHalf*(20-1)+1:BlockLengthHalf*(20+1)),2,0,OutputType="llr"),[BlockLengthHalf 2]);
                LlrActual = LLR{K}';
                decodedData{K} = ldpcDecode(LlrActual(:), cfgLDPCDec, 10);
                
                for j=1:K
                    ErrCount(j) = biterr(data{j}, decodedData{j});
                end
                error_count(k_idx) = sum(ErrCount)/K;
                
                if any(ErrCount ~= 0)
                    [~, j_max] = max(ErrCount);
                    place(j_max) = place(j_max) + 1;
                end
            end
            TurboErr(end+1) = mean(error_count)/n;
        end
        
        %% EGYEDI FÁJL MENTÉSE
        % Itt mentjük ki az aktuális mátrix eredményeit egy KÜLÖN .mat fájlba
        outFilePath = fullfile(outputDir, [baseFileName, '_result.mat']);
        save(outFilePath, 'TurboErr', 'BlockLengthHalf', 'n');
        fprintf('  -> Sikeresen mentve: %s\n', outFilePath);
        
    catch ME
        fprintf('  -> HIBA a %s fájl feldolgozása közben: %s\n', files(i).name, ME.message);
    end
end
fprintf('--- Összes fájl feldolgozva! ---\n');

end