MODULE fdtd

      USE numerics
      USE source

      IMPLICIT NONE
      ! Déclaration de variables

      ! Class fdtd
      TYPE :: cnfdtd
                  ! Variables locales
                  INTEGER, ALLOCATABLE :: N_d(:)
                  INTEGER, ALLOCATABLE :: S(:)
                  REAL(8), ALLOCATABLE :: dx, dy, dt
                  REAL(8), ALLOCATABLE :: Ex(:,:), Ey(:,:), Hz(:,:)
                  REAL(8), ALLOCATABLE :: B(:), rhs(:) 
                  REAL(8), ALLOCATABLE :: J(:,:)
                  REAL(8), ALLOCATABLE :: A(:,:)
                  REAL(8), ALLOCATABLE :: c_E(:,:), c_H(:,:)
                  REAL(8) :: a1, a2, bx, by
            CONTAINS
                  ! Méthodes
                  PROCEDURE :: init
                  PROCEDURE :: compute_fdtd
                  PROCEDURE :: freememory
      END TYPE cnfdtd

      CONTAINS

      SUBROUTINE init(cn)
            CLASS(cnfdtd), INTENT(inout) :: cn
            INTEGER :: i

            ! Initialisation des variables
            ALLOCATE(cn%N_d (        0:10       ) )                                      ! Grid sampling densities
            ALLOCATE(cn%S   (        0:50       ) )                                      ! Courant Number 
            ALLOCATE(cn%B   (  0 : 2 * (Nx + 1) * (Ny + 1) - 1        ) )               ! Pour matrice A entiere
            ALLOCATE(cn%rhs (  0 : 2 * (Nx - 1) * (Ny - 1) - 1        ) )               ! Pour matrice A intérieur
            ALLOCATE(cn%Ex  (                    0:Nx, 0:Ny                               ) )
            ALLOCATE(cn%Ey  (                    0:Nx, 0:Ny                               ) )
            ALLOCATE(cn%J   (                0 : Nx  ,  0 : Ny                            ) )
            ALLOCATE(cn%Hz  (                 0 : Nx , 0:Ny                               ) )
            !ALLOCATE(cn%A   (   0 : 2 * (Nx + 1) - 1 , 0: 2 * (Ny + 1) - 1                ) )              ! Précédente allocation de A
            ALLOCATE(cn%A   (     0:2 * (Nx + 1) - 1 , 0:2 * (Ny + 1) - 1                 ) )              ! Matrice A entière
            ALLOCATE(cn%c_E (                    0:Nx, 0:Ny                               ) )
            ALLOCATE(cn%c_H (                    0:Nx, 0:Ny                               ) )

            cn%N_d = (/ (10*i, i = 0,10) /)
            ! PRINT *, 'N_d = ', cn%N_d
            ! print *, 'size(N_d)' , size(cn%N_d)

            cn%S = (/ (2*i, i = 0,50) /)
            ! PRINT *, 'S = ', cn%S
            ! print *, 'size(S)' , size(cn%S)


            cn%dx = (c / fmax) / cn%N_d(3) 
            cn%dy = cn%dx
            WRITE(*, '(/,T5,A,ES17.3, /)') 'dx = ', cn%dx


            cn%dt = 0.98d0 / ( c * sqrt(  1.0d0 / (cn%dx * cn%dx)  + 1.0d0 / (cn%dy * cn%dy) ) )
            WRITE(*, '(/,T5,A,ES17.3, /)') 'dt = ', cn%dt

            WRITE(*, '(/,T5, A,F17.12, /)') 'Courant number cdt/dx = ', abs(c * cn%dt / cn%dx)

            cn%a1 = cn%dt / (2.d0 * epsilon_0) 
            cn%a2 = cn%dt / (2.d0 * mu_0)
            cn%bx = c * cn%dt / (2.d0 * cn%dx)
            cn%by = c * cn%dt / (2.d0 * cn%dy)
            WRITE(*, '(/,T5,A,ES17.3, /)') 'bx = ',cn%bx
            WRITE(*, '(/,T5,A,ES17.3, /)') 'by = ',cn%by
            WRITE(*, '(/,T5,A,ES17.3, /)') 'a1 = ', cn%a1
            WRITE(*, '(/,T5,A,ES17.3, /)') 'a1/dx = ', cn%a1 / cn%dx
            WRITE(*, '(/,T5,A,ES17.3, /)') 'a2 = ', cn%a2
            WRITE(*, '(/,T5,A,ES17.3, /)') 'a2/dx = ', cn%a2 / cn%dx

            ! Initialisation des champs
            cn%A = 0.d0
            cn%B = 0.d0
            cn%J = 0.d0
            cn%Ex = 0.d0
            cn%Ey = 0.d0
            cn%Hz = 0.d0
            cn%c_E = 1.0d0 / (epsilon_0 * cn%dx)
            cn%c_H = 1.0d0 / (mu_0 * cn%dx)


      END SUBROUTINE init

      SUBROUTINE compute_fdtd(cn)
            CLASS(cnfdtd), INTENT(inout) :: cn
            ! Variables locales
            CHARACTER(LEN=500) :: charac 
            LOGICAL :: display_it
            INTEGER :: info
            INTEGER :: n, m, nrhs, nvec, nrow, ncol, nblock
            INTEGER :: i,j
            INTEGER :: i0,j0,i1,j1
            INTEGER :: snapshot
            INTEGER :: idx_Ex, idx_Ey
            REAL(8), ALLOCATABLE :: Exx(:,:)
            REAL(8), ALLOCATABLE :: Eyy(:,:)
            REAL(8), ALLOCATABLE :: Exy(:,:)
            REAL(8), ALLOCATABLE :: Eyx(:,:)
            REAL(8), ALLOCATABLE :: B_mat(:,:)
            INTEGER :: ipiv(SIZE(cn%A,1))       ! Sert de pivot

            ALLOCATE(Exx(0:Nx, 0:Ny))
            ALLOCATE(Eyy(0:Nx, 0:Ny))
            ALLOCATE(Exy(0:Nx, 0:Ny))
            ALLOCATE(Eyx(0:Nx, 0:Ny))
            ALLOCATE(B_mat(0:2 * Nx + 1, 0:Ny))

            Exx = 0.d0; Eyy = 0.d0; Exy = 0.d0; Eyx = 0.d0; B_mat = 0.d0; ipiv = 0

            WRITE(*,'(/,T5,A,I5)') "Nx = ", Nx
            WRITE(*,'(/, T5, A, I5X, I5)') "shape(A) = "          ,  shape(cn%A)
            WRITE(*,'(/, T5, A, I5X, I5)') "shape(A_i) = "        ,  shape(Exx)
            WRITE(*,'(/, T5, A, I15X)')    "shape(B) = "          ,  shape(cn%B)
            WRITE(*,'(/, T5, A, I5X, I5)') "shape(B_mat) = "      ,  shape(B_mat)
            WRITE(*,'(/, T5, A, I5X, I5)') "shape(Ex) = "         ,  shape(cn%Ex)
            WRITE(*,'(/, T5, A, I5X, I5)') "shape(Ey) = "         ,  shape(cn%Ey)
            WRITE(*,'(/, T5, A, I5X, I5)') "shape(Hz) = "         ,  shape(cn%Hz)
            WRITE(*,'(/, T5, A, I5X, I5)') "shape(J) = "          ,  shape(cn%J)
            WRITE(*,'(/, T5, A, I5X, I5)') "shape(ipiv) = "       ,  shape(ipiv)


            m = 0
            display_it = .FALSE.  
            charac = "" 

            !-------------------------------------------------------------!
            !------------------ Ecriture de la matrice A -----------------!
            !-------------------------------------------------------------!

                        ! -------- ! -------- !
                        !   Exx    !   Exy    !
            !  A =      !----------!----------!
                        !   Eyx    !   Eyy    !
                        ! -------- ! -------- !



            ! Remplissage de la matrice A
            ! CALCUL DIRECT DE LA MATRICE A

                    !!!    


            !----------------------------------------------------!
            !------------------ Sous matrice Exx -----------------!
            !----------------------------------------------------!
                  Exx(0,0) = 1.0d0 + 2.d0 * cn%bx**2
                  Exx(0,1) = - cn%bx**2 

                  DO j = 1, Nx
                        Exx(j,j-1) = - cn%bx**2
                        Exx(j,j) = 1.0d0 + 2.d0 * cn%bx**2
                        Exx(j,j+1) = - cn%bx**2
                  END DO

                  


            ! ! Affichage de la matrice Exx
            IF (display_it) THEN
                  CALL display_matrix(Exx, "Exx")
            END IF

            



            !----------------------------------------------------!
            !------------------ Sous matrice Eyy -----------------!
            !----------------------------------------------------!     
                  Eyy = 0.d0
                  Eyy(0,0) = 1.0d0 + 2.d0 * cn%by**2
                  Eyy(0,1) = - cn%by**2

                  DO i = 1, Nx - 1
                        Eyy(i,i-1) = - cn%by**2 
                        Eyy(i,i)   = 1.0d0 + 2.d0 * cn%by**2
                        Eyy(i,i+1) = - cn%by**2
                  END DO

                  Eyy(Nx, Nx-1) =  - cn%by**2
                  Eyy(Nx, Nx) = 1.0d0 + 2.d0 * cn%by**2
                  Eyy(Nx - 1, 0) = 0.d0


            ! ! Affichage de la matrice Eyy
            IF (display_it) THEN
                  CALL display_matrix(Eyy, "Eyy")
            END IF




            !----------------------------------------------------!
            !------------------ Sous matrice Exy -----------------!
            !----------------------------------------------------!

                  Exy(0,0) = -1.d0
                  Exy(1,0) =  2.d0
                  Exy(2,0) = -1.d0

                  DO i = 1, Nx-2
                        Exy(i,i)         =  -1.d0
                        Exy(i + 1 , i)   =   2.d0
                        Exy(i + 2 , i)   =  -1.d0
                  END DO 

                  Exy(Nx-1, Nx - 1) = -1.d0
                  Exy(Nx, Nx - 1) = 2.d0
                  Exy(Nx, Nx) = -1.d0


                  Exy = cn%bx * cn%by * Exy
                  

            ! ! Affichage de la matrice Exy
            IF (display_it) THEN
                  CALL display_matrix(Exy, "Exy")
            END IF

            !----------------------------------------------------!
            !------------------ Sous matrice Eyx -----------------!
            !----------------------------------------------------!

                  Eyx = transpose(Exy)



                  !Exy = 0.0d0

            ! Affichage de la matrice Eyx
            IF (display_it) THEN
                  CALL display_matrix(Eyx, "Eyx")
            END IF
            

                  
            
             
            !---------------------------------------------------------------!
            !------------------ Assemblage de la matrice A -----------------!
            !---------------------------------------------------------------!

                  ! -------- ! -------- !
                  !   Exx    !   Exy    !
                  !----------!----------!
                  !   Eyx    !   Eyy    !
                  ! -------- ! -------- !


            ! Détermine les indices de collage
            i0 = 0;        j0 = 0
            i1 = Nx + 1;   j1 = Ny + 1
            

            WRITE(*, '(/, T5, A, I5, I5)') "i0, j0 = ", i0, j0
            WRITE(*, '(/, T5, A, I5, I5)') "i1, j1 = ", i1, j1

            ! Collage des blocs diagonaux
            cn%A(i0  :i0 + Nx, j0  :j0 + Ny)   = Exx
            cn%A(i1  :i1 + Nx, j1  :j1 + Ny)   = Eyy
            ! Collage des matrices de couplage
            cn%A(i0 : i0 + Nx, j1 : j1 + Ny) = Exy
            cn%A(i1 : i1 + Nx, j0 : j0 + Ny) = Eyx


            IF (display_it) then
                  CALL display_matrix(cn%A, " A assemblée")
                  write(*, '(/,t5,A)') " Extraction de la matrice intérieur A :"
            ENDIF




            ! Affichage de la matrice A
            ! CALL display_matrix(cn%A)
            !---------------------------------------------------!

            ! Vérification de la symétrie de la matrice A
            !CALL matrix_sym(cn%A)
            ! ! !---------------------------------------------------!

            ! Libération de mémoire 
            IF (ALLOCATED(Exx)) THEN
                  DEALLOCATE(Exx)
            END IF
            IF (ALLOCATED(Eyy)) THEN
                  DEALLOCATE(Eyy)
            END IF
            IF (ALLOCATED(Exy)) THEN
                  DEALLOCATE(Exy)
            END IF
            IF (ALLOCATED(Eyx)) THEN
                  DEALLOCATE(Eyx)
            END IF
            ! ! !---------------------------------------------------!




            ! -------------------------------------------------------------------!
            ! ------------------ Décomposition LU de A --------------------!
            ! -------------------------------------------------------------------!

            CALL DGETRF(size(cn%A,1), SIZE(cn%A,2),cn%A, size(cn%A,1),ipiv, info)
            IF (info > 0) THEN
                  WRITE(*,'(/,T5,A,I0,A,I0,A,/)') 'U(', info , ',', info ,') is exactly zero. The factorization has been completed, but the factor U is exactly singular.'
                  STOP 'LU failed'
            ELSE IF (info < 0) THEN
                  WRITE(*,'(T5,A,I0,A,/)') 'The ',info,'-th argument had an ilegal value.'
                  STOP 'LU failed'
            END IF
            

            


            ! Ouverture du fichier de sortie
            OPEN(idfile , file = "data/Ex.txt", status = "replace", action = "write", form = "formatted")
            OPEN(idfile + 1 , file = "data/Hz.txt", status = "replace", action = "write", form = "formatted")
            
            !-------------------------------------------------------------!
            !------------------- Boucle temporelle -----------------------!
            !-------------------------------------------------------------!
            WRITE(*, '(/, T5, "Injection de la source en ", I5, I5)') i_src, j_src
            WRITE(*, '(/, T5, A, /)') "Début de la boucle temporelle"
            snapshot = 20

           
            nrhs = 1
            m = 0
            nblock = (Nx + 1) * (Ny + 1) ! Nombre de blocs dans la matrice A
            DO n = 0, Nt - 1

                  IF (MOD(n,5*snapshot) == 0) THEN
                        WRITE(*, '(/, T5, "itération temporelle : ",I4)') n
                  END IF

                  cn%J(i_src, j_src) =  Esrc(n)

                  !--------------------------------------------------------------!
                  !------------------- Ecriture du vecteur B --------------------!
                  !--------------------------------------------------------------!
                  ! On enregistre les résultats du temps précédent
                  cn%Ex = B_mat(0 : Nx, 0 : Ny)
                  cn%Ey = B_mat(i1 : i1 + Nx, 0 : Ny)             ! i1 = Nx + 1

                  !print * , "pass 1"

                 
                  ! Vecteur secoind membre B
                  

                  ! Second membre Ex
                  DO i = 0,  Nx
                        !print *, "i = ", i
                        DO j = 0, Ny
                              ! Détermine le bonne indice
                              idx_Ex = i * (Nx + 1) + j
                              ! print *, "idx_Ex = ", idx_Ex, "i,j =", i , j
                              
                              if (i == 0 .OR. i == Nx .OR. j == 0 .OR. j == Ny) then
                                    cn%B(idx_Ex) = 0.d0
                              else
                                    cn%B(idx_Ex) =      (1.d0 - 2.d0 * cn%bx**2) * cn%Ex(i,j)                         & 
                                                + cn%bx**2 * ( cn%Ex(i, j - 1) + cn%Ex(i, j + 1) )                    &
                                                - cn%bx*cn%by * ( cn%Ey(i + 1, j) - cn%Ey(i, j) )         &
                                                + cn%bx*cn%by * ( cn%Ey(i + 1 , j -1) - cn%Ey(i, j - 1) )         &
                                                + 2.d0 * cn%a1 * (cn%Hz(i,j+1) - cn%Hz(i, j-1))
                              endif
                        END DO
                  END DO

                  !print *, "pass 2"


                  ! Second membre Ey
                  DO i = 0 , Nx
                        !print *, "i = ", i
                        DO j = 0,  Ny
                              ! Détermine le bonne indice
                              idx_Ey = (Nx+1)*(Ny+1) + i * (Nx + 1) + j
                              ! print *, "idx_Ey = ", idx_Ey, 'i,j =', i , j
                              if (i == 0 .OR. i == Nx .OR. j == 0 .OR. j == Ny) then 
                                    cn%B(idx_Ey) = 0.d0
                              else
                              ! Calcul du second membre Ey
                                    cn%B(idx_Ey) =      (1.d0 - 2.d0 * cn%by**2)*cn%Ey(i,j)                              & 
                                                + cn%by**2 * ( cn%Ey(i - 1, j) + cn%Ey(i + 1, j)    )                    &
                                                - cn%bx*cn%by * ( cn%Ex(i , j + 1) - cn%Ex(i , j)  )         &
                                                + cn%bx*cn%by * ( cn%Ex(i-1, j + 1)    - cn%Ex(i-1, j) )                 &
                                                - 2.d0 * cn%a1 * (cn%Hz(i+1,j) - cn%Hz(i-1, j))
                               END IF
                        END DO
                  END DO
                  !print *, "pass 3"

                  

                  ! Résolution du système linéaire
                  CALL DGETRS('N', SIZE(cn%A,1), nrhs, cn%A, SIZE(cn%A,1), ipiv, cn%B, SIZE(cn%B), info)
                 
                  !print *, "pass 4"
                  
                  !CALL display_matrix(B_mat, "B_mat")
                  
                  ! Ex dans le premier bloc
                   B_mat(    0:Nx   , 0:Ny) = RESHAPE( cn%B(0 :    nblock-1), &
                                                 SHAPE=[Nx+1,Ny+1], ORDER=[2,1] )

                  ! ! Ey dans le second bloc
                   B_mat(Nx+1:2*Nx+1, 0:Ny) = RESHAPE( cn%B(nblock : 2*nblock-1), &
                                                       SHAPE=[Nx+1,Ny+1], ORDER=[1,2] )

                  

                  


                  !print *, "pass 4"




                  ! Mise à jour explicite de Hz

                  DO i = 1, Nx-1
                        DO j = 1, Ny-1
                              cn%Hz(i,j) = cn%Hz(i,j) + cn%a2 / cn%dy * ( B_mat(i,j + 1) - B_mat(i,j - 1)                      &
                                                                        + cn%Ex(i, j + 1) - cn%Ex(i,j-1) )               &
                                                      - cn%a2 / cn%dx * ( B_mat(i1 + (i + 1),j) - B_mat(i1 + (i-1),j)          &          ! i1 = Nx + 1
                                                                        + cn%Ey(i + 1, j) - cn%Ey(i - 1,j) )
                        END DO
                  END DO

                  ! Injection de la source 
                  cn%Hz(i_src,j_src) = Esrc(n)

                  !print *, "pass 7"

                  


                  

                  ! Ecriture dans le fichier 
                  IF (MOD(n,snapshot) == 0) THEN
                        m = m + 1
                        DO i = 0, Nx, 2
                              DO j = 0, Ny, 2
                                    WRITE(idfile + 1, '(F0.15,1X)', advance='no') cn%Hz(i,j)
                              END DO
                              WRITE(idfile + 1, *)
                        END DO 
                        WRITE(idfile + 1, *)    
                  END IF
                  ! ! !---------------------------------------------------!



                  
                  
            END DO

            ! WRITE(*, '(/, T5, A, /)') "Test reshape du vecteur B :"
            ! print *, "shape(B_mat) = ", shape(B_mat)

            ! WRITE(*,'(2(AX,F16.10))') 'B(0)=',cn%B(0),' B_mat(0,0)=',B_mat(0,0)
            ! WRITE(*,'(2(AX,F16.10))') 'B(1)=',cn%B(1),' B_mat(0,1)=',B_mat(0,1)
            ! WRITE(*,'(2(AX,F16.10))') 'B(19)=',cn%B(19),' B_mat(3,3)=',B_mat(3,3)
            ! WRITE(*,'(2(AX,F16.10))') 'B(Ny+1)=',cn%B(Ny+1),' B_mat(1,0)=',B_mat(1,0)

            WRITE(*, '(/, t5, A, I5)') "Nombre de blocs : ", m
            

            CLOSE(idfile)
            CLOSE(idfile + 1)

            
            WRITE(*, '(/, T5, A, /)') "Fin de la boucle temporelle"


      END SUBROUTINE compute_fdtd



      SUBROUTINE freememory(cn)
            CLASS(cnfdtd), INTENT(inout) :: cn

            ! Libération de la mémoire
            IF (ALLOCATED(cn%N_d)) THEN
            DEALLOCATE(cn%N_d)
            END IF
            IF (ALLOCATED(cn%S)) THEN
            DEALLOCATE(cn%S)
            END IF
            IF (ALLOCATED(cn%Ex)) THEN
            DEALLOCATE(cn%Ex)
            END IF
            IF (ALLOCATED(cn%Hz)) THEN
            DEALLOCATE(cn%Hz)
            END IF
            IF (ALLOCATED(cn%c_E)) THEN
            DEALLOCATE(cn%c_E)
            END IF
            IF (ALLOCATED(cn%c_H)) THEN
            DEALLOCATE(cn%c_H)
            END IF
            IF (ALLOCATED(cn%A)) THEN
            DEALLOCATE(cn%A)
            END IF
            IF (ALLOCATED(cn%B)) THEN
            DEALLOCATE(cn%B)
            END IF
      

      END SUBROUTINE freememory


      SUBROUTINE matrix_sym(A)
            REAL(8), INTENT(in) :: A(:,:)
            INTEGER :: i,j
            LOGICAL :: is_symmetric

            is_symmetric = .TRUE.

            DO i = 1, size(A,1)
                  DO j = 1, size(A,2)
                        IF (abs(A(i,j) -  A(j,i)) > eps ) THEN
                              is_symmetric = .FALSE.
                              EXIT
                        END IF
                  END DO
            END DO
            IF (is_symmetric) THEN
                  WRITE(*, '(/, T5, A, /)') "La matrice A est symétrique."
            ELSE
                  WRITE(*, '(/, T5, A, /)') "La matrice A n'est pas symétrique."
            END IF
      ENDSUBROUTINE matrix_sym

      SUBROUTINE display_matrix(A, name)
            ! Affiche la matrice A
            REAL(8), INTENT(in) :: A(:,:)
            CHARACTER(LEN=*), INTENT(in), OPTIONAL :: name
            INTEGER :: i

            IF (PRESENT(name)) THEN
                  WRITE(*, '(/, T5, A,A,A /)', advance = 'no') "Matrice ", name, " : "
            ELSE
                  WRITE(*, '(/, T5, A, /)', advance = 'no') "Matrice :"
            END IF
            DO i = LBOUND(A,1), UBOUND(A,1)
                  WRITE(*, '(I5,500F12.2)') i-1, A(i,:)
            END DO
      ENDSUBROUTINE display_matrix

      ! FUNCTION idx_Ex(i,j)
      !       INTEGER, INTENT(in) :: i,j
      !       INTEGER :: idx_Ex

      !       idx_Ex = i * (Ny - 1) + j - 1

      ! END FUNCTION idx_Ex

      ! FUNCTION idx_Ey(i,j)
      !       INTEGER, INTENT(in) :: i,j
      !       INTEGER :: idx_Ey

      !       idx_Ey = Nx * (Ny - 1) + (i - 1) * Ny + j

      ! END FUNCTION idx_Ey

       SUBROUTINE extract_matrix_ud(A_int,A)
            ! ARGUMENTS
            REAL(8), INTENT(in), DIMENSION(:,:) :: A
            REAL(8), INTENT(inout), DIMENSION(:,:), ALLOCATABLE :: A_int

            ! VARIABLES LOCALES
            INTEGER :: idx_min
            INTEGER :: idx_max
            INTEGER :: idy_min
            INTEGER :: idy_max

            ! Initialisation 
            idx_min = LBOUND(A,1)          ! Indice inférieur de la dimension x 
            idx_max = UBOUND(A,1)          ! Indice supérieur de la dimension x

            idy_min = LBOUND(A,2)          ! Indice inférieur de la dimension y
            idy_max = UBOUND(A,2)          ! Indice supérieur de la dimension y

            ! Allocation en retirant les indices de bord en y
            ALLOCATE(A_int(idx_min : idx_max, idy_min + 1 : idy_max - 1))

            PRINT *, "idx_min, idx_max = ", idx_min, idx_max
            PRINT *, "idy_min, idy_max = ", idy_min, idy_max
            PRINT *, "shape(A_int) = ", shape(A_int)

            A_int = 0.d0

            ! Extraction de la matrice A
            A_int = A(idx_min + 1 : idx_max - 1, : ) 
      ENDSUBROUTINE extract_matrix_ud

      SUBROUTINE extract_matrix_lr(A_int,A)
            ! ARGUMENTS
            REAL(8), INTENT(in), DIMENSION(:,:) :: A
            REAL(8), INTENT(inout), DIMENSION(:,:), ALLOCATABLE :: A_int

            ! VARIABLES LOCALES
            INTEGER :: idx_min
            INTEGER :: idx_max
            INTEGER :: idy_min
            INTEGER :: idy_max

            ! Initialisation 
            idx_min = LBOUND(A,1)          ! Indice inférieur de la dimension x 
            idx_max = UBOUND(A,1)          ! Indice supérieur de la dimension x

            idy_min = LBOUND(A,2)          ! Indice inférieur de la dimension y
            idy_max = UBOUND(A,2)          ! Indice supérieur de la dimension y

            ! Allocation en retirant les indices de bord en y
            ALLOCATE(A_int(idx_min : idx_max, idy_min + 1 : idy_max - 1))

            PRINT *, "idx_min, idx_max = ", idx_min, idx_max
            PRINT *, "idy_min, idy_max = ", idy_min, idy_max
            PRINT *, "shape(A_int) = ", shape(A_int)

            A_int = 0.d0

            ! Extraction de la matrice A
            A_int = A( : , idy_min + 1 : idy_max - 1) 
      ENDSUBROUTINE extract_matrix_lr

END MODULE fdtd