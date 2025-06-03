Nx      = 500
Ny      = 500
step    = 2
ncol    = int(Nx/step)+1
nrow    = int(Ny/step)+1
nblocks = 51             # nombre de cartes

# on garde la palette par défaut
set palette rgb 33,13,10
set view map

set terminal qt size 800,600
set grid
set xlabel "x"
set ylabel "y"
set title "Animation Hz"

# boucle sur chaque carte index k
do for [k=0:nblocks-1] {
    splot 'data/Hz.txt' index k matrix with image
    pause 0.001
}
