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
                  REAL(8), ALLOCATABLE :: Hx(:,:), Hy(:,:), Ez(:,:)
                  REAL(8), ALLOCATABLE :: rhs(:)
                  REAL(8), ALLOCATABLE :: A(:,:)
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
            ALLOCATE(cn%Hx  (                    0:Nx, 0:Ny                               ) )
            ALLOCATE(cn%Hy  (                    0:Nx, 0:Ny                               ) )
            ALLOCATE(cn%rhs   (                      0 : 2 * (Nx + 1) * (Ny + 1) - 1        ) )               ! Pour matrice A entiere
            ALLOCATE(cn%Ez  (                 0 : Nx , 0:Ny                               ) )
            ALLOCATE(cn%A   (   0 : 2 * (Nx + 1) - 1 , 0: 2 * (Ny + 1) - 1                ) )


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
            cn%rhs = 0.d0
            cn%Hx = 0.d0
            cn%Hy = 0.d0
            cn%Ez = 0.d0


      END SUBROUTINE init

      SUBROUTINE compute_fdtd(cn)
            CLASS(cnfdtd), INTENT(inout) :: cn
            ! Variables locales
            CHARACTER(LEN=500) :: charac 
            LOGICAL :: display_it
            INTEGER :: info
            INTEGER :: n, m, nrhs, nvec, nrow, ncol
            INTEGER :: i,j
            INTEGER :: i0,j0,i1,j1
            INTEGER :: snapshot
            INTEGER :: idx_Hx, idx_Hy
            REAL(8), ALLOCATABLE :: Hxx(:,:)
            REAL(8), ALLOCATABLE :: Hyy(:,:)
            REAL(8), ALLOCATABLE :: Hxy(:,:)
            REAL(8), ALLOCATABLE :: Hyx(:,:)
            REAL(8), ALLOCATABLE :: B_pec(:,:)

            INTEGER :: ipiv(SIZE(cn%A,1))       ! Sert de pivot

            ALLOCATE(Hxx(0:Nx, 0:Ny))
            ALLOCATE(Hyy(0:Nx, 0:Ny))
            ALLOCATE(Hxy(0:Nx, 0:Ny))
            ALLOCATE(Hyx(0:Nx, 0:Ny))
            ALLOCATE(B_pec(0:2 * Nx + 1, 0:Ny))

            Hxx = 0.d0; Hyy = 0.d0; Hxy = 0.d0; Hyx = 0.d0; B_pec = 0.d0; ipiv = 0

            WRITE(*,'(/,T5,A,I5)') "Nx = ", Nx
            WRITE(*,'(/, T5, A, I5X, I5)') "shape(A) = "    ,  shape(cn%A)
            WRITE(*,'(/, T5, A, I5X, I5)') "shape(A_i) = "  ,  shape(Hxx)
            WRITE(*,'(/, T5, A, I15X)')    "shape(rhs) = "    ,  shape(cn%rhs)
            WRITE(*,'(/, T5, A, I5X, I5)') "shape(B_pec) = ",  shape(B_pec)
            WRITE(*,'(/, T5, A, I5X, I5)') "shape(Hx) = "   ,  shape(cn%Hx)
            WRITE(*,'(/, T5, A, I5X, I5)') "shape(Hy) = "   ,  shape(cn%Hy)
            WRITE(*,'(/, T5, A, I5X, I5)') "shape(Ez) = "   ,  shape(cn%Ez)
            WRITE(*,'(/, T5, A, I5X, I5)') "shape(ipiv) = " ,  shape(ipiv)

            
            m = 0
            display_it = .FALSE.  
            charac = ""   

            !!!!! Implémentation de CRANK NICOLSON standard !!!!!


            !-------------------------------------------------------------!
            !------------------ Ecriture de la matrice A -----------------!
            !-------------------------------------------------------------!

                        ! -------- ! -------- !
                        !   Hxx    !   Hxy    !
            ! A =       !----------!----------!
                        !   Hyx    !   Hyy    !
                        ! -------- ! -------- !

            !----------------------------------------------------!
            !------------------ Sous matrice Hxx -----------------!
            !----------------------------------------------------!
                  Hxx(0,0) = 1.0d0 + 2.d0 * cn%bx**2
                  Hxx(0,1) = - cn%bx**2 

                  DO j = 1, Nx
                        Hxx(j,j-1) = - cn%bx**2
                        Hxx(j,j) = 1.0d0 + 2.d0 * cn%bx**2
                        Hxx(j,j+1) = - cn%bx**2
                  END DO



            ! ! Affichage de la matrice Hxx
            IF (display_it) THEN
                  CALL display_matrix(Hxx, "Hxx")
            END IF

            



            !----------------------------------------------------!
            !------------------ Sous matrice Hyy -----------------!
            !----------------------------------------------------!     
            Hyy = Hxx

            ! ! Affichage de la matrice Hyy
            IF (display_it) THEN
                  CALL display_matrix(Hyy, "Hyy")
            END IF




            !----------------------------------------------------!
            !------------------ Sous matrice Hxy -----------------!
            !----------------------------------------------------!

                  DO i = 0, Nx-2
                        j = i
                        Hxy(i , j)     =  -1.d0
                        Hxy(i , j + 1) =   2.d0
                        Hxy(i , j + 2) =  -1.d0
                  END DO 

                  Hxy(Nx-1,Nx-1)    = - 1.d0
                  Hxy(Nx-1, Nx)     =   2.d0
                  Hxy(Nx, Nx)       = - 1.d0

                  Hxy = - cn%bx * cn%by * Hxy
                  



           

            !----------------------------------------------------!
            !------------------ Sous matrice Hyx -----------------!
            !----------------------------------------------------!

                  Hyx = transpose(Hxy)
                  Hxy = Hyx
                  Hyx = transpose(Hxy)




             ! ! Affichage de la matrice Hxy
            IF (display_it) THEN
                  CALL display_matrix(Hxy, "Hxy")
            END IF

            ! Affichage de la matrice Hyx
            IF (display_it) THEN
                  CALL display_matrix(Hyx, "Hyx")
            END IF
            

                  
            
             
            !---------------------------------------------------------------!
            !------------------ Assemblage de la matrice A -----------------!
            !---------------------------------------------------------------!

                        ! -------- ! -------- !
                        !   Hxx    !   Hxy    !
            ! A =       !----------!----------!
                        !   Hyx    !   Hyy    !
                        ! -------- ! -------- !


            ! Détermine les indices de collage
            i0 = 0;        j0 = 0
            i1 = Nx + 1;   j1 = Ny + 1
            

            WRITE(*, '(/, T5, A, I5, I5)') "i0, j0 = ", i0, j0
            WRITE(*, '(/, T5, A, I5, I5)') "i1, j1 = ", i1, j1

            ! Collage des blocs diagonaux
            cn%A(i0  :i0 + Nx, j0  :j0 + Ny)   = Hxx
            cn%A(i1  :i1 + Nx, j1  :j1 + Ny)   = Hyy
            ! Collage des matrices de couplage
            cn%A(i0 : i0 + Nx, j1 : j1 + Ny) = Hxy
            cn%A(i1 : i1 + Nx, j0 : j0 + Ny) = Hyx





            IF (display_it) then
                  CALL display_matrix(cn%A, " A assemblée")
            ENDIF

            ! Vérification de la symétrie de la matrice A
            !CALL matrix_sym(cn%A)
            ! ! !---------------------------------------------------!

            ! Libération de mémoire 
            IF (ALLOCATED(Hxx)) THEN
                  DEALLOCATE(Hxx)
            END IF
            IF (ALLOCATED(Hyy)) THEN
                  DEALLOCATE(Hyy)
            END IF
            IF (ALLOCATED(Hxy)) THEN
                  DEALLOCATE(Hxy)
            END IF
            IF (ALLOCATED(Hyx)) THEN
                  DEALLOCATE(Hyx)
            END IF
            ! ! !---------------------------------------------------!




            ! -------------------------------------------------------------------!
            ! ------------------ Décomposition LU de A --------------------!
            ! -------------------------------------------------------------------!

            CALL DGETRF(size(cn%A,1), SIZE(cn%A,2),cn%A, size(cn%A,1),ipiv, info)
            !CALL DGETRF(size(A_int,1), size(A_int,2), A_int1, size(A_int,1), ipiv, info)
            IF (info > 0) THEN
                  WRITE(*,'(/,T5,A,I0,A,I0,A,/)') 'U(', info , ',', info ,') is exactly zero. The factorization has been completed, but the factor U is exactly singular.'
                  STOP 'LU failed'
            ELSE IF (info < 0) THEN
                  WRITE(*,'(T5,A,I0,A,/)') 'The ',info,'-th argument had an ilegal value.'
                  STOP 'LU failed'
            END IF
            

            


            ! Ouverture du fichier de sortie
            OPEN(idfile , file = "data/Hx.txt", status = "replace", action = "write", form = "formatted")
            OPEN(idfile + 1 , file = "data/Ez.txt", status = "replace", action = "write", form = "formatted")
            
            !-------------------------------------------------------------!
            !------------------- Boucle temporelle -----------------------!
            !-------------------------------------------------------------!
            WRITE(*, '(/, T5, "Injection de la source en ", I5, I5)') i_src, j_src
            WRITE(*, '(/, T5, A, /)') "Début de la boucle temporelle"
            snapshot = 10

            nrow = 2 * (Nx - 1)
            ncol = Ny - 1
            nvec = nrow * ncol
            nrhs = 1
            m = 0
            DO n = 0, Nt - 1

                  IF (MOD(n,5*snapshot) == 0) THEN
                        WRITE(*, '(/, T5, "itération temporelle : ",I4)') n
                  END IF

                  !--------------------------------------------------------------!
                  !------------------- Ecriture du vecteur rhs ------------------!
                  !--------------------------------------------------------------!
                   ! On enregistre les résultats du temps n
                  cn%Hx = B_pec(0 : Nx, 0 : Ny)
                  cn%Hy = B_pec(i1 : i1 + Nx, 0 : Ny)                   !i1 = Nx + 1
                  

                  ! On parcourt l'entierté des champs Hx et Hy
                  ! Second membre Hx
                  DO i = 1,  Nx-1
                        !print *, "i = ", i)
                        DO j = 1, Ny-1
                              ! Détermine le bonne indice
                              idx_Hx = i * (Nx + 1) + j
                              !print *, "idx_Hx = ", idx_Hx, "i,j =", i , j
                              cn%rhs(idx_Hx) =      (1.d0 - 2.d0 * cn%bx**2) * cn%Hx(i,j)                       & 
                                          + cn%bx**2 * ( cn%Hx(i, j - 1) + cn%Hx(i, j + 1) )                    &
                                          - cn%bx*cn%by * ( cn%Hy(i + 1,  j)     - cn%Hy(i, j) )                &
                                          + cn%bx*cn%by * ( cn%Hy(i + 1 , j - 1) - cn%Hy(i, j - 1) )            &
                                          - 2.d0 * cn%a2 * (cn%Ez(i,j) - cn%Ez(i, j-1))
                        END DO
                  END DO

                  ! ! Conditoon aux bords du champ magnétique
                  ! DO i = 0, Nx
                  !       idx_Hx         = i * (Nx + 1)
                  !       j = 0
                  !       !print *, "idx_Hx = ", idx_Hx, "i, j = ", i, j
                  !       cn%rhs(idx_Hx) = cn%rhs(idx_Hx + 1) 
                  ! END DO




                  ! Second membre Hy
                  DO i = 1 , Nx-1
                        !print *, "i = ", i
                        DO j = 1,  Ny - 1
                              ! Détermine le bonne indice
                              idx_Hy = (Nx+1)*(Ny+1) + i * (Nx + 1) + j
                              !print *, "idx_Hy = ", idx_Hy, 'i,j =', i , j
                              ! Calcul du second membre Hy
                              cn%rhs(idx_Hy) =      (1.d0 - 2.d0 * cn%by**2)*cn%Hy(i,j)                              & 
                                          + cn%by**2 * ( cn%Hy(i - 1, j) + cn%Hy(i + 1, j)    )                    &
                                          - cn%bx*cn%by * ( cn%Hx(i , j + 1)  - cn%Hx(i , j)  )                    &
                                          + cn%bx*cn%by * ( cn%Hx(i-1, j + 1) - cn%Hx(i-1, j) )                    &
                                          + 2.d0 * cn%a2 * (cn%Ez(i,j) - cn%Ez(i-1, j))
                        END DO
                  END DO

                  ! ! Condition aux bords du champ magnétique
                  ! DO j = 0, Ny
                  !       idx_Hy = (Nx + 1) * (Ny + 1) + j
                  !       i = 0
                  !       !print *, ""
                  !       !print *, "idx_Hy = ", idx_Hy, "i, j = ", i, j
                  !       cn%rhs(idx_Hy) = cn%rhs(idx_Hy + 1) 
                  ! END DO

                                                                           ! Mise à jour explicite de Ez
                  DO i = 1, Nx-1
                        DO j = 1, Ny-1
                              cn%Ez(i,j) = cn%Ez(i,j) - cn%a1 / cn%dy * ( B_pec(i,j + 1) - B_pec(i,j)                      &
                                                                        + cn%Hx(i, j + 1) - cn%Hx(i,j) )                   &
                                                      + cn%a1 / cn%dx * ( B_pec(i1 + (i + 1),j) - B_pec(i1 + i,j)          &          ! i1 = Nx + 1
                                                                        + cn%Hy(i + 1, j) - cn%Hy(i,j) )
                        END DO
                  END DO

                  !Injection de source
                  cn%Ez(i_src,j_src) = Esrc(n)        

                                                      ! CDT DE BORD / PMC
                  cn%Ez(: ,0)  = 0.d0          ! Bord inférieur
                  cn%Ez(: ,Ny) = 0.d0          ! Bord supérieur
                  cn%Ez(0 ,:)  = 0.d0          ! Bord gauche
                  cn%Ez(Nx,:)  = 0.d0          ! Bord droit


                  ! Résolution du système linéaire pour les champs H au temps n + 1
                  CALL DGETRS('N', SIZE(cn%A,1), nrhs, cn%A, SIZE(cn%A,1), ipiv, cn%rhs, SIZE(cn%rhs), info)
                  !print *, "pass 4"

                  ! reshape du vecteur rhs / order = [2,1] fait varier j avant i
                  B_pec = reshape(cn%rhs, shape = [ 2 * (Nx + 1), Ny + 1], order = [2, 1])




                  
                  ! Ecriture dans le fichier 
                  IF (MOD(n,snapshot) == 0) THEN
                        m = m + 1
                        DO i = 0, Nx
                              DO j = 0, Ny
                                    WRITE(idfile + 1, '(F0.15,1X)', advance='no') cn%Ez(i,j)
                                    !write(idfile    , '(F0.15,1X)', advance='no') cn%Hx(i,j)
                              END DO
                              WRITE(idfile + 1, *)
                              !write(idfile    , *)
                        END DO 
                        WRITE(idfile + 1, *)    
                        WRITE(idfile    , *)
                  END IF
                  ! ! !---------------------------------------------------!



                  
                  
            END DO

            ! WRITE(*, '(/, T5, A, /)') "Test reshape du vecteur rhs :"
            ! print *, "shape(B_pec) = ", shape(B_pec)

            ! WRITE(*,'(2(AX,F16.10))') 'rhs(0)=',cn%rhs(0),' B_pec(0,0)=',B_pec(0,0)
            ! WRITE(*,'(2(AX,F16.10))') 'rhs(1)=',cn%rhs(1),' B_pec(0,1)=',B_pec(0,1)
            ! WRITE(*,'(2(AX,F16.10))') 'rhs(19)=',cn%rhs(19),' B_pec(3,3)=',B_pec(3,3)
            ! WRITE(*,'(2(AX,F16.10))') 'rhs(Ny+1)=',cn%rhs(Ny+1),' B_pec(1,0)=',B_pec(1,0)

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
            IF (ALLOCATED(cn%Hx)) THEN
            DEALLOCATE(cn%Hx)
            END IF
            IF (ALLOCATED(cn%Ez)) THEN
            DEALLOCATE(cn%Ez)
            END IF
            IF (ALLOCATED(cn%A)) THEN
            DEALLOCATE(cn%A)
            END IF
            IF (ALLOCATED(cn%rhs)) THEN
            DEALLOCATE(cn%rhs)
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