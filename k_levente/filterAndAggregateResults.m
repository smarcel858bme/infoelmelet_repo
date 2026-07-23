function filterAndAggregateResults(resultsDir, finalOutputDir, filterCondition)
% FILTERANDAGGREGATERESULTS Megkeresi a különálló .mat fájlokat, és egy
% tetszőleges (szöveges) feltétel alapján szűri őket.
%
% Paraméterek:
%   resultsDir      - A különálló .mat fájlokat tartalmazó mappa
%   finalOutputDir  - Ide menti az összegző fájlt a megfelelő eredményekkel
%   filterCondition - (String) Logikai feltétel a szűréshez.
%                     Példa: 'n < 100' vagy 'BlockLengthHalf == 66'

if nargin < 1 || isempty(resultsDir), resultsDir = fullfile(pwd, 'Results'); end
if nargin < 2 || isempty(finalOutputDir), finalOutputDir = fullfile(pwd, 'Final_Output'); end
if nargin < 3 || isempty(filterCondition), filterCondition = 'true'; end % Ha nincs feltétel, mindent megtart

if ~exist(finalOutputDir, 'dir')
    mkdir(finalOutputDir);
end

files = dir(fullfile(resultsDir, '*.mat'));

if isempty(files)
    fprintf('Nem találtam .mat fájlt a mappában: %s\n', resultsDir);
    return;
end

CompiledData = struct();
validFileCount = 0;

fprintf('--- .mat fájlok szűrése a következő feltétellel: "%s" ---\n', filterCondition);
for i = 1:length(files)
    currentFile = fullfile(files(i).folder, files(i).name);

    % Betöltjük a 3 változót
    loadedData = load(currentFile, 'TurboErr', 'BlockLengthHalf', 'n');

    % Változók átemelése a lokális munkaterületre, hogy az eval() lássa őket
    if isfield(loadedData, 'n'), n = loadedData.n; else n = NaN; end
    if isfield(loadedData, 'BlockLengthHalf'), BlockLengthHalf = loadedData.BlockLengthHalf; else BlockLengthHalf = NaN; end
    if isfield(loadedData, 'TurboErr'), TurboErr = loadedData.TurboErr; else TurboErr = []; end

    % Dinamikus feltétel kiértékelése
    try
        feltetelTeljesult = eval(filterCondition);
    catch ME
        fprintf('  [HIBA] Nem sikerült kiértékelni a feltételt a %s fájlnál. Hiba: %s\n', files(i).name, ME.message);
        feltetelTeljesult = false;
    end

    if feltetelTeljesult
        validFileCount = validFileCount + 1;

        % Biztonságos név a nagy struktúrához
        [~, safeName, ~] = fileparts(files(i).name);
        safeName = matlab.lang.makeValidName(safeName); 

        CompiledData.(safeName).TurboErr = TurboErr;
        CompiledData.(safeName).BlockLengthHalf = BlockLengthHalf;
        CompiledData.(safeName).n = n;

        fprintf('  [OK] Megfelelt: %s\n', files(i).name);
    else
        fprintf('  [KIZÁRVA] Nem felelt meg: %s\n', files(i).name);
    end
end

% Ha volt sikeres találat, kimentjük őket
if validFileCount > 0
    finalFile = fullfile(finalOutputDir, 'Filtered_Aggregated_Results.mat');
    save(finalFile, 'CompiledData');
    fprintf('\nSikeresen kimentve %d fájl adata ide: %s\n', validFileCount, finalFile);
else
    fprintf('\nEgyetlen fájl sem felelt meg a megadott szűrési feltételnek!\n');
end

end