%% Központi Vezérlő Szkript
% Ez a szkript összefogja és sorban meghívja az összes korábbi modulunkat.

clear; clc;

% 1. Kipróbáljuk a mátrix beolvasót egyetlen fájlon
disp('--- 1. Lépés: Egy mátrix beolvasása teszt jelleggel ---');
% Itt meghívjuk a readTurboCode.m fájlt! Bemenet a fájlnév, kimenet a H mátrix.
H = readTurboCode('LTE_TC_N132_K40.txt'); 
disp('A mátrix sikeresen beolvasva!');

% 2. Kötegelt (Batch) szimuláció indítása
disp('--- 2. Lépés: Kötegelt szimuláció ---');
% Meghívjuk a batchSimulateMatrices.m fájlt.
% Feltételezzük, hogy van egy 'Bemeneti_Matrixok' mappád a .txt fájlokkal.
% Az eredményeket egy 'Nyers_Eredmenyek' nevű mappába kérjük.
bemenet = 'bemeneti_matrixok';
kimenet = 'nyers_eredmenyek';
snrTarto = 1:0.5:5;
keretSzam = 200;

batchSimulateMatrices(bemenet, kimenet, snrTarto, keretSzam);

% 3. Eredmények szűrése és összesítése
disp('--- 3. Lépés: Eredmények szűrése ---');
% Meghívjuk a filterAndAggregateResults.m fájlt.
% Kiszűrjük a Nyers_Eredmenyek mappából azokat a .mat fájlokat, ahol n > 10.
% A végeredményt a 'Szurt_Vegeredmeny' mappába mentjük.
feltetel = 'n > 10';
filterAndAggregateResults(kimenet, 'szurt_vegeredmeny', feltetel);

disp('--- MINDEN FOLYAMAT SIKERESEN BEFEJEZŐDÖTT! ---');