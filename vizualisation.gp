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

block1 = 5
block2= 15
block3 = 20

splot 'data/Ez.txt' index block1 matrix with image
splot 'data/Ez.txt' index block2 matrix with image 
splot 'data/Ez.txt' index block3 matrix with image
#splot 'data/Hx.txt' index block1 matrix with image
#splot 'data/Hx.txt' index block2 matrix with image
#splot 'data/Hx.txt' index block3 matrix with image

