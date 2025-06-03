% paramètres – à adapter selon votre code Fortran
Nx       = 199;           % Nb de points en x moins 1
Ny       = 199;           % idem en y
snapshot = 20;            % pas de snapshot dans votre boucle Fortran
Nt       = 1001;          % nb d’itérations temporelles
% on calcule le nombre de snapshots réellement écrits :
n_block      = floor((Nt-1)/snapshot) + 1;
% dans chaque bloc vous écrivez Nz = (Nx/2+1)×(Ny/2+1) valeurs
nrow_sample  = floor(Nx/2) + 1;
ncol_sample  = floor(Ny/2) + 1;
n_per_block  = nrow_sample * ncol_sample;

fprintf('Expected blocks: %d   samples per block: %d\n', ...
        n_block, n_per_block);

% 1) recharge le vecteur brut
V = load('data/Hz.txt');     % doit faire n_block*n_per_block × 1
assert(numel(V)==n_block*n_per_block, ...
       'Wrong file size: got %d, expected %d', ...
        numel(V), n_block*n_per_block);

% 2) reshape en matrice [n_per_block × n_block]
M = reshape(V, n_per_block, n_block);

% 3) choisir le k-ème snapshot (ligne k)
k = 5;
b = M(:,k);

% 4) reconstruire la grille 2D
% on avait parcouru i=0:2:Nx → nrow_sample éléments
%             j=0:2:Ny → ncol_sample éléments
B = reshape(b, nrow_sample, ncol_sample);

% 5) afficher
figure;
set(gcf, 'Position', [350,200,800,600]);
imagesc(0:2:Nx, 0:2:Ny, B');  % attention au transpose selon orientation
axis equal tight;
xlabel('x');
ylabel('y');
title(sprintf('Hz à t-snapshot #%d (n = %d)', k, (k-1)*snapshot));
colorbar;