% Ez csak egy function amit be kell copy-zni a Matlabba, hogy a másik kód működjön

function H = read_alist(filename)
% 1. Fájl megnyitása olvasásra
fid = fopen(filename, 'r');
if fid == -1
    error('Nem sikerült megnyitni a fájlt!');
end

% 2. Dimenziók beolvasása: N (oszlopok) és M (sorok)
dims = fscanf(fid, '%d', 2);
N = dims(1); 
M = dims(2);

% 3. Maximális fokszámok: max_dv (oszlop), max_dc (sor)
max_deg = fscanf(fid, '%d', 2);

% 4. Az egyes csomópontok pontos fokszáma 
% (Csak átugorjuk a kurzorral, de ki kell olvasni a bufferből)
dv = fscanf(fid, '%d', N);
dc = fscanf(fid, '%d', M);

% 5. A kapcsolatok (1-esek) indexeinek gyűjtése
row_idx = [];
col_idx = [];

for i = 1:N
    % Oszloponként beolvassuk a csatlakozó sorok indexeit.
    % Az alist fájl nullákkal egészíti ki az üres helyeket max_dv-ig.
    cols = fscanf(fid, '%d', max_deg(1));

    % Csak a 0-nál nagyobb, valós sor-indexeket tartjuk meg
    valid_rows = cols(cols > 0);

    % Mátrix koordináták építése
    row_idx = [row_idx; valid_rows];
    col_idx = [col_idx; repmat(i, length(valid_rows), 1)];
end

fclose(fid);

% 6. Ritka (sparse) mátrix létrehozása
% A ones() vektor garantálja, hogy minden megadott koordinátára 1-es kerül
H = sparse(row_idx, col_idx, ones(length(row_idx), 1), M, N);

% logikai GF(2) formátum kényszerítése
H = logical(H); 
end
