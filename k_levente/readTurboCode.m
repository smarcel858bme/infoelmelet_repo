function H_sparse = readTurboCode(filename)
% READTURBOCODE Beolvas egy sűrű bináris mátrixot tartalmazó .txt fájlt.
% Bemenet: filename - A .txt fájl útvonala és neve (string)
% Kimenet: H_sparse - Logikai, ritka (sparse) mátrix

% 1. Sűrű mátrix beolvasása a .txt fájlból
H_dense = readmatrix(filename);

% 2. Átalakítás logikai és ritka (sparse) formátumba
H_sparse = sparse(logical(H_dense));

end