Nx      = 500
Ny      = 500
step    = 2
ncol    = int(Nx/step)+1
nrow    = int(Ny/step)+1
nblocks = 251             # nombre de cartes

# on garde la palette par défaut
set palette rgb 33,13,10
set view map

set terminal qt size 800,600
set grid
set xlabel "x"
set ylabel "y"
set title "Animation Hz"

block1 = 2
block2= 25
block3 = 50

splot 'data/Hz.txt' index block1 matrix with image
splot 'data/Hz.txt' index block2 matrix with image 
splot 'data/Hz.txt' index block3 matrix with image
#splot 'data/Ex.txt' index block1 matrix with image
#splot 'data/Ex.txt' index block2 matrix with image
#splot 'data/Ex.txt' index block3 matrix with image

