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
            ALLOCATE(cn%Ex  (                    0:Nx, 0:Ny                               ) )
            ALLOCATE(cn%Ey  (                    0:Nx, 0:Ny                               ) )
            ALLOCATE(cn%rhs   (                      0 : 2 * (Nx + 1) * (Ny + 1) - 1        ) )               ! Pour matrice A entiere
            ALLOCATE(cn%Hz  (                 0 : Nx , 0:Ny                               ) )
            ALLOCATE(cn%A   (   0 : 2 * (Nx + 1) - 1 , 0: 2 * (Ny + 1) - 1                ) )


            cn%N_d = (/ (10*i, i = 0,10) /)
            ! PRINT *, 'N_d = ', cn%N_d
            ! print *, 'size(N_d)' , size(cn%N_d)



            cn%S = (/ (2*i, i = 0,50) /)
            ! PRINT *, 'S = ', cn%S
            ! print *, 'size(S)' , size(cn%S)


            cn%dx = (c / fmax) / mesh_density
            cn%dy = cn%dx
            WRITE(*, '(/,T5,A,ES17.3, /)') 'dx = ', cn%dx


            cn%dt = CFL / ( c * sqrt(  1.0d0 / (cn%dx * cn%dx)  + 1.0d0 / (cn%dy * cn%dy) ) )
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
            cn%Ex = 0.d0
            cn%Ey = 0.d0
            cn%Hz = 0.d0


      END SUBROUTINE init

      FUNCTION function_idx(i,j,n_var)
            INTEGER, INTENT(in) :: i,j,n_var
            INTEGER :: function_idx

            function_idx = (n_var - 1) * (Nx + 1) * (Ny + 1) + j * (Nx + 1) + i
            RETURN
      END FUNCTION function_idx

      SUBROUTINE compute_fdtd(cn)
            CLASS(cnfdtd), INTENT(inout) :: cn
            ! Variables locales
            CHARACTER(LEN=500) :: charac 
            LOGICAL :: display_it
            INTEGER :: info
            INTEGER :: n, m, n_var, n_elt, A_row, A_col
            INTEGER :: i,j, idx,idy,i_var
            INTEGER :: i0,j0,i1,j1
            INTEGER,ALLOCATABLE :: ipiv(:)
            REAL(8), ALLOCATABLE :: A(:, :)
            REAL(8), ALLOCATABLE :: rhs(:)
            REAL(8) :: t_start, t_end, t_elapsed

            ! Variable d'indices pour 3 variables inconnus Ex, Ey, Hz de dimension (Nx+1)*(Ny+1) chacune
            n_var = 3
            A_row = n_var * (Nx + 1)
            A_col = n_var * (Ny + 1)
            n_elt = (Nx + 1) * (Ny + 1)   ! Nombre d'éléments par variable inconnue
            n_elt = n_var * n_elt
            WRITE(*, '(/,A,I10)') "Nombre d'éléments par variable inconnue n_elt = ", n_elt
            

            ALLOCATE(A(0: A_row - 1, 0: A_col - 1))
            ALLOCATE(ipiv(0: A_row - 1))


            A = 0.d0;
            write(*, '(/,A,I5,I5,/)') "shape(A) = ", shape(A)
            WRITE(*, '(/,A,I5,I5,/)') "shape(ipiv) = ", shape(ipiv)

            OPEN(100, file = "data/params.txt", status = "replace", action = "write")
                  WRITE(100,*) Nx, Ny, Nt, cn%dx, cn%dy, cn%dt, snapshot, mesh_density, CFL, c, i_src, j_src
            CLOSE(100)
      

            
            m = 0
            display_it = .FALSE.  
            charac = ""   

            !!!!!!! Implémentation de CRANK NICOLSON standard !!!!!!!


            !-------------------------------------------------------------!
            !------------------ Ecriture de la matrice A -----------------!
            !-------------------------------------------------------------!

                        ! -------- ! -------- !-------!
                        !   A11    !   A12    !  A13  !
            ! A =       !----------!----------!-------!
                        !   A21    !   A22    !  A23  !
                        ! -------- ! -------- !-------!
                        !   A31    !   A32    !  A33  !
                        ! -------- ! -------- !-------!


            DO i = 0,  n_elt - 1
                  A(i,i) = 1.d0
            END DO

            DO i = 0, Nx
                  DO j = 0, Ny
                        !BLOC A13
                        idy = 2 * (Nx + 1) + j
                        IF (j == i) THEN
                              A(i,idy) = - cn%a1 / cn%dy
                              IF (j > 0) A(i,idy -1) = cn%a1 / cn%dy
                        END IF

                        ! Bloc A23
                        idx = (Nx +1) + i
                        IF (j == i) THEN
                              A(idx,idy) = cn%a1 / cn%dx
                              IF (i > 0) A(idx - 1,idy) = - cn%a1 / cn%dx
                        END IF

                        ! Bloc A31
                        idx = 2 * (Nx + 1) + i 
                        idy = j
                        IF (j == i) THEN
                              A(idx,idy) = - cn%a2 / cn%dy
                              IF (j < Ny) A(idx ,idy + 1) = cn%a2 / cn%dy
                        END IF

                        ! Bloc A32
                        idy = (Ny + 1) + j
                        IF (j == i) THEN
                              A(idx,idy) = - cn%a2 / cn%dx
                              IF (i < Nx) A(idx + 1,idy) = cn%a2 / cn%dx
                        END IF
                  END DO
            END DO

            !CALL display_matrix(A, "A")

            ! -------------------------------------------------------------!
            ! ------------ Décomposition LU de la matrice A ---------------!
            ! -------------------------------------------------------------!

            CALL DGETRF(A_row, A_col, A, A_row, ipiv, info)
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

            m = 0
            ALLOCATE(rhs(0 : n_var * n_elt - 1))
            rhs = 0.d0
            WRITE(*, '(/,A,I10)') "shape(rhs) = ", size(rhs)
            CALL cpu_time(t_start)

            DO n = 0, Nt - 1

                  IF (MOD(n,5 * snapshot) == 0) THEN
                        WRITE(*, '(/,A,I10,/)') "Time step n = ", n
                  END IF

                  

                  ! Récupération des champs Ex, Ey, Hz
                  DO j = 0, Ny
                        DO i = 0, Nx
                              cn%Ex(i,j) = rhs(function_idx(i,j,1))
                              cn%Ey(i,j) = rhs(function_idx(i,j,2))
                              cn%Hz(i,j) = rhs(function_idx(i,j,3))
                        END DO
                  END DO


                  DO j = 0, Ny
                        DO i = 0, Nx
                              !BLOC rhs1
                              idx = function_idx(i,j,1)
                              !print *, 'idx = ', idx
                              rhs(idx) = cn%Ex(i,j) + cn%a1 / cn%dy * (                         &
                                                                cn%Hz(i,j) - cn%Hz(i,j-1) )     
                                               

                              !BLOC rhs2
                              idx = function_idx(i,j,2)
                             ! print *, 'idx = ', idx
                              rhs(idx) = cn%Ey(i,j) - cn%a1 / cn%dx * (                         &
                                                                cn%Hz(i,j) - cn%Hz(i-1,j) )

                              !BLOC rhs3
                              idx = function_idx(i,j,3)
                              !print *, 'idx = ', idx
                              rhs(idx) = cn%Hz(i,j) + cn%a2 / cn%dy * ( cn%Ex(i,j+1) - cn%Ex(i,j) ) &
                                                    - cn%a2 / cn%dx * ( cn%Ey(i+1,j) - cn%Ey(i,j) )
                        END DO
                  END DO


                  ! Injection de la source
                  rhs(function_idx(i_src,j_src,3)) = rhs(function_idx(i_src,j_src,3)) + Esrc(n)

                  ! Résolution du système linéaire A * x = rhs
                  CALL DGETRS('N', A_row, 1, A, A_row, ipiv, rhs, A_row, info)
                  IF (info < 0) THEN
                        WRITE(*,'(T5,A,I0,A,/)') 'The ',info,'-th argument had an ilegal value.'
                        STOP 'Solve failed'
                  END IF



                   ! Ecriture dans le fichier 
                  IF (MOD(n,snapshot) == 0) THEN
                        m = m + 1
                        DO i = 0, Nx
                              DO j = 0, Ny
                                    WRITE(idfile + 1, '(F0.15,1X)', advance='no') cn%Hz(i,j)
                                    write(idfile    , '(F0.15,1X)', advance='no') cn%Ex(i,j)
                              END DO
                              WRITE(idfile + 1, *)
                              write(idfile    , *)
                        END DO 
                        WRITE(idfile + 1, *)    
                        WRITE(idfile    , *)
                  END IF


            END DO

            CALL cpu_time(t_end)
            t_elapsed = t_end - t_start
            WRITE(*, '(/, T5, A, F8.4, A, /)') "Temps de calcul total : ", t_elapsed, " secondes"

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
            IF (ALLOCATED(cn%Ey)) THEN
            DEALLOCATE(cn%Ey)
            END IF
            IF (ALLOCATED(cn%Hz)) THEN
            DEALLOCATE(cn%Hz)
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
                  WRITE(*, '(I5,500ES12.2)') i-1, A(i,:)
            END DO
      ENDSUBROUTINE display_matrix

END MODULE fdtd