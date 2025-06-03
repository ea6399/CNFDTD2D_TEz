%% Afficahge de la cartographie au temps n

% Paramètres
Nx = 199;
Nx = (Nx + 1) / 2 ;
Ny = Nx;
Nt = 1001;
snapshot = 20;
n_block = (Nt - 1) / snapshot + 1;
fprintf('Nombre de block : %d\n', n_block);


n_per_snapshot = (Nx) * (Ny);
M_inter = load("data/Hz.txt");
fprintf('Size after loading the file : %d %d\n', size(M_inter))
M = reshape(M_inter,n_per_snapshot,n_block);
disp('shape(M) =');
disp(size(M));

% Choix du temps à afficher
k = 5;                              % Numéro de bloc
M_t = M(:,k);                       % Extraction de l'instant à afficher
M_t = reshape(M_t, Nx, Ny); 
figure(1);
set(gcf, 'Position', [350,250,800,600]);
imagesc(M_t);
title('Hz')
colorbar;



