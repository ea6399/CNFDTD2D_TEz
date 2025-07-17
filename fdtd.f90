MODULE fdtd

      USE numerics
      USE source

      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'dmumps_struc.h'

      ! Déclaration de variables
      ! Class fdtd
      TYPE :: cnfdtd
                  ! Variables locales
                  INTEGER, ALLOCATABLE :: N_d(:)
                  INTEGER, ALLOCATABLE :: S(:)
                  REAL(8), ALLOCATABLE :: dx, dy, dt
                  REAL(8), ALLOCATABLE :: Ex(:,:), Ey(:,:), Hz(:,:)
                  REAL(8), ALLOCATABLE :: B(:)
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
      ! MUMPS
      TYPE(dmumps_struc) :: mumps

      CONTAINS

      SUBROUTINE init(cn)
            CLASS(cnfdtd), INTENT(inout) :: cn
            INTEGER :: i

            ! Initialisation des variables
            ALLOCATE(cn%N_d (        0:10       ) )                                      ! Grid sampling densities
            ALLOCATE(cn%S   (        0:50       ) )                                      ! Courant Number 
            ALLOCATE(cn%B   (  0 : 2 * (Nx + 1) * (Ny + 1) - 1        ) )               ! Pour matrice A entiere
            ALLOCATE(cn%Ex  (                    0:Nx, 0:Ny                               ) )
            ALLOCATE(cn%Ey  (                    0:Nx, 0:Ny                               ) )
            ALLOCATE(cn%J   (                0 : Nx  ,  0 : Ny                            ) )
            ALLOCATE(cn%Hz  (                 0 : Nx , 0:Ny                               ) )
            ALLOCATE(cn%A   (     0:2 * (Nx + 1) * (Ny + 1) - 1 , 0:2*(Ny + 1)*(Nx + 1) - 1    ) )              ! Matrice A entière
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

      INTEGER FUNCTION id_Ex(i,j)
            INTEGER, INTENT(in) :: i,j

            id_Ex = i * (Nx + 1) + j

      END FUNCTION id_Ex

      INTEGER FUNCTION id_Ey(i,j)
            INTEGER, INTENT(in) :: i,j

            id_Ey =  (Nx + 1) * (Ny + 1) + i * (Nx+1) + j

      END FUNCTION id_Ey

      SUBROUTINE compute_fdtd(cn)
            ! Class
            CLASS(cnfdtd), INTENT(inout) :: cn
      
            ! Variables locales
            CHARACTER(LEN=500) :: charac 
            LOGICAL :: display_it
            INTEGER :: info
            INTEGER :: n, m, nvec, nrow, ncol
            INTEGER :: i,j, idx, idy, k, id_nz
            INTEGER :: i1
            INTEGER :: snapshot
            REAL(8), ALLOCATABLE :: X(:,:)
            ! Mumps variables
            INTEGER(8) :: nmps, nnz, irn, jcn 
            REAL(8), ALLOCATABLE :: diag_x(:), diag_y(:), subdiag_x(:), subdiag_y(:)
            REAL(8), ALLOCATABLE :: diag_xy(:), updiag_xy(:), uupdiag_xy(:) 
            INTEGER :: counter_nnz


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

            ! Allocation et initialisation des diagonales
            ALLOCATE(diag_x(0: Nx))
            ALLOCATE(diag_y(0: Nx))
            ALLOCATE(subdiag_x(0 : Nx - 1))
            ALLOCATE(subdiag_y(0 : Nx - 1))
            diag_x = 0.d0
            diag_y = 0.d0
            subdiag_x = 0.d0
            subdiag_y = 0.d0

            ALLOCATE(diag_xy(0:Nx))
            ALLOCATE(updiag_xy(0:Nx - 1))
            ALLOCATE(uupdiag_xy(0:Nx - 2))
            diag_xy = 0.d0
            updiag_xy = 0.d0
            uupdiag_xy = 0.d0

            ! Computation des coefficients de la matrice A
            diag_x = 1.d0 + 2.d0 * cn%bx**2
            diag_y = 1.d0 + 2.d0 * cn%by**2
            subdiag_x = - cn%bx**2
            subdiag_y = - cn%by**2

            diag_xy = - cn%bx * cn%by
            updiag_xy = 2.d0 *  cn%bx * cn%by
            uupdiag_xy = - cn%bx * cn%by


            !------------------------------------------------------------------!
            !------------------------ Entering MUMPS --------------------------!
            !------------------------------------------------------------------!

            WRITE(*,'(/,/,/,/,/,"Entering MUMPS Solver",/,/,/,/,/)')

            mumps%COMM = 0                ! 0 pour séquentiel
            mumps%SYM  = 2                ! 2 pour symétrique général
            mumps%PAR  = 1                ! 1 : L'host est le seul processeur

            !-----------------------------!
            ! Initialisation d'un package !
            !-----------------------------!
            ! On résout AX = B, B ayant Ny + 1 colonnes, et 2 * (Nx + 1 ) lignes
            ! | Ex |
            ! | Ey |
            !        ______ _________
            !       |       |        |
            !       |  Exx  |  Exy   |
            !       |       |        |
            !   A = |----------------|
            !       |       |        |    
            !       |  Eyx  |  Eyy   |
            !       |_______|_______ |

            mumps%JOB = -1               ! Initialisation de MUMPS
            CALL DMUMPS(mumps)

            ! Parametrage de la matrice A
            nrow = 2 * (Nx + 1)
            ncol = 2 * (Ny + 1)
            mumps%N = nrow  ! Nombre de lignes 

                          ! Element diagonaux Exx et Eyy     +     ! Element matrice de couplage Eyx
            mumps%NNZ = 2 * (Nx + 1 + Nx)  +  (Nx + 1 + Nx + Nx - 1)  ! On ne stocke que les coefficients de la partie inférieure de la matrice
            nnz = mumps%NNZ
            PRINT *, "Nombre element non nuls : ", nnz

            ! Allocation des tableaux et initialisation
            ALLOCATE(mumps%IRN(0:nnz - 1))
            ALLOCATE(mumps%JCN(0:nnz - 1))
            ALLOCATE(mumps%A(0:nnz - 1))
            mumps%IRN = 0
            mumps%JCN = 0
            mumps%A = 0.d0
            
            ! Entrée des éléments Non Nuls
            counter_nnz = 0

            ! Matrice Exx
            DO i = 0, Nx
                  ! Éléments diagonaux
                  mumps%IRN(counter_nnz) = i 
                  mumps%JCN(counter_nnz) = i 
                  mumps%A(counter_nnz)   = 1 + 2.d0 * cn%bx**2
                  counter_nnz = counter_nnz + 1
                  ! Éléments sub-diagonaux
                  if (i < Nx) then
                        mumps%IRN(counter_nnz) = i + 1
                        mumps%JCN(counter_nnz) = i
                        mumps%A(counter_nnz)   = - cn%bx**2
                        counter_nnz = counter_nnz + 1
                  end if
            END DO

            ! Matrice Eyy
            DO i = 0, Nx
                  ! Éléments diagonaux
                  mumps%IRN(counter_nnz) = Nx + 1 + i 
                  mumps%JCN(counter_nnz) = Ny + 1 + i 
                  mumps%A(counter_nnz)   = 1 + 2.d0 * cn%by**2
                  counter_nnz = counter_nnz + 1
                  
                  ! Éléments sub diagonaux
                  IF (i < Nx) then
                        mumps%IRN(counter_nnz) = Nx + 1 + i + 1
                        mumps%JCN(counter_nnz) = Ny + 1 + i 
                        mumps%A(counter_nnz)   = - cn%by**2
                        counter_nnz = counter_nnz + 1
                  ENDIF 
            END DO


            ! Matrice Eyx
            DO i = 0 , Nx
                  ! Éléments diagonaux
                  mumps%IRN(counter_nnz) = Nx + 1 + i 
                  mumps%JCN(counter_nnz) = i
                  mumps%A(counter_nnz)   = - cn%bx * cn%by
                  counter_nnz = counter_nnz + 1

                  ! Éléments up-diagonaux
                  IF (i > 0) THEN
                  mumps%IRN(counter_nnz) = Nx + 1 + i - 1
                  mumps%JCN(counter_nnz) = i
                  mumps%A(counter_nnz)   = 2.d0 * cn%bx * cn%by
                  counter_nnz = counter_nnz + 1
                  END IF

                  ! Éléments uup-diagonaux
                  IF (i > 1) THEN
                  mumps%IRN(counter_nnz) = Nx + 1 + i - 2
                  mumps%JCN(counter_nnz) = i
                  mumps%A(counter_nnz)   = - cn%bx * cn%by
                  counter_nnz = counter_nnz + 1
                  END IF
            END DO


            OPEN(500, file = "data/mumps_elt.txt", status = "replace", action = "write", form = "formatted")
                  WRITE(500, *) mumps%IRN
                  WRITE(500, *) mumps%JCN
                  WRITE(500, *) mumps%A
            CLOSE(500)

            WRITE(*, '(/,/)')

            ! Analyse MUMPS
            mumps%JOB = 1
            mumps%ICNTL(1:3) = 0
            CALL DMUMPS(mumps)

            ! Factorisation MUMPS
            mumps%JOB = 2
            CALL DMUMPS(mumps)



            ! Ouverture du fichier de sortie
            OPEN(idfile , file = "data/Ex.txt", status = "replace", action = "write", form = "formatted")
            OPEN(idfile + 1 , file = "data/Hz.txt", status = "replace", action = "write", form = "formatted")
            
            !-------------------------------------------------------------!
            !------------------- Boucle temporelle -----------------------!
            !-------------------------------------------------------------!
            WRITE(*, '(/, T5, "Injection de la source en ", I5, I5)') i_src, j_src
            WRITE(*, '(/, T5, A, /)') "Début de la boucle temporelle"
            snapshot = 20

            m = 0

            ! Initialisation du RHS pour MUMPS
            ALLOCATE(mumps%RHS(0 : nrow - 1))
            mumps%RHS = 0.d0

             DO n = 0, Nt - 1

                  IF (MOD(n,5*snapshot) == 0) THEN
                        WRITE(*, '(/, T5, "itération temporelle : ",I4)') n
                  END IF

                  !--------------------------------------------------------------!
                  !------------------- Ecriture du vecteur rhs --------------------!
                  !--------------------------------------------------------------!

                   ! On enregistre les résultats du temps précédent
                  cn%Ex = B_pec(0 : Nx, 0 : Ny)
                  cn%Ey = B_pec(i1 : i1 + Nx, 0 : Ny)                   !i1 = Nx + 1


                  ! On parcourt l'entierté des champs Hx et Hy
                  ! Second membre Hx
                  DO i = 0,  Nx-1
                        !print *, "i = ", i
                        DO j = 1, Ny-1
                              ! Détermine le bonne indice
                              idx_Hx = i * (Nx + 1) + j
                              ! print *, "idx_Hx = ", idx_Hx, "i,j =", i , j
                              mumps%RHS(idx_Hx) =      (1.d0 - 2.d0 * cn%bx**2) * cn%Hx(i,j)                       &
                                          + cn%bx**2 * ( cn%Hx(i, j - 1) + cn%Hx(i, j + 1) )                    &
                                          - cn%bx*cn%by * ( cn%Hy(i + 1,  j)     - cn%Hy(i, j) )                &
                                          + cn%bx*cn%by * ( cn%Hy(i + 1 , j - 1) - cn%Hy(i, j - 1) )            &
                                          - 2.d0 * cn%a2 * (cn%Ez(i,j) - cn%Ez(i, j-1))
                        END DO
                  END DO

                  ! Conditoon aux bords du champ magnétique
                  DO i = 0, Nx
                        idx_Hx         = i * (Nx + 1)
                        j = 0
                        !print *, "idx_Hx = ", idx_Hx, "i, j = ", i, j
                        cn%rhs(idx_Hx) = cn%rhs(idx_Hx + 1)
                  END DO




                  ! Second membre Hy
                  DO i = 1 , Nx-1
                        !print *, "i = ", i
                        DO j = 0,  Ny - 1
                              ! Détermine le bonne indice
                              idx_Hy = (Nx+1)*(Ny+1) + i * (Nx + 1) + j
                              ! print *, "idx_Hy = ", idx_Hy, 'i,j =', i , j
                              ! Calcul du second membre Hy
                              cn%rhs(idx_Hy) =      (1.d0 - 2.d0 * cn%by**2)*cn%Hy(i,j)                              &
                                          + cn%by**2 * ( cn%Hy(i - 1, j) + cn%Hy(i + 1, j)    )                    &
                                          - cn%bx*cn%by * ( cn%Hx(i , j + 1)  - cn%Hx(i , j)  )                    &
                                          + cn%bx*cn%by * ( cn%Hx(i-1, j + 1) - cn%Hx(i-1, j) )                    &
                                          + 2.d0 * cn%a2 * (cn%Ez(i,j) - cn%Ez(i-1, j))
                        END DO
                  END DO

                  ! Condition aux bords du champ magnétique
                  DO j = 0, Ny
                        idx_Hy = (Nx + 1) * (Ny + 1) + j
                        i = 0
                        !print *, ""
                        !print *, "idx_Hy = ", idx_Hy, "i, j = ", i, j
                        cn%rhs(idx_Hy) = cn%rhs(idx_Hy + 1)
                  END DO




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


                  ! Résolution du système linéaire
                  CALL DGETRS('N', SIZE(cn%A,1), nrhs, cn%A, SIZE(cn%A,1), ipiv, cn%rhs, SIZE(cn%rhs), info)
                  !print *, "pass 4"

                  ! reshape du vecteur rhs / order = [2,1] fait varier j avant i
                  B_pec = reshape(cn%rhs, shape = [ 2 * (Nx + 1), Ny + 1], order = [2, 1])

                  ! Ecriture dans le fichier
                  IF (MOD(n,snapshot) == 0) THEN
                        m = m + 1
                        DO i = 0, Nx, 2
                              DO j = 0, Ny, 2
                                    WRITE(idfile + 1, '(F0.15,1X)', advance='no') cn%Ez(i,j)
                                    write(idfile    , '(F0.15,1X)', advance='no') cn%Hx(i,j)
                              END DO
                              WRITE(idfile + 1, *)
                              write(idfile    , *)
                        END DO
                        WRITE(idfile + 1, *)
                        WRITE(idfile    , *)
                  END IF
                  ! ! !---------------------------------------------------!


            END DO

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

      


END MODULE fdtd