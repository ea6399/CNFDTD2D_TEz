%% Afficahge de la cartographie au temps n

% Paramètres
Nx = 501;
Ny = Nx;
Nt = 1001;
snapshot = 5;
n_block = (Nt - 1) / snapshot;


M_inter = load("data/Hz.txt");
M = reshape(M_inter, [n_block,(Nx + 1) * (Ny + 1)]);
disp('shape(M) =');
disp(size(M));

% Choix du temps à afficher
k = 5;                          % Numéro de bloc
M_t = M(k * (Nx+1),Ny + 1);     % Extraction de l'instant à afficher
figure(1);
set(gcf, 'Position', [350,250,800,600]);
imagesc(M_t);
title('Hz')
colorbar;



