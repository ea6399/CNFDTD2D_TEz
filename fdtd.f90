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

      FUNCTION function_idx(i,j,n_var)
            INTEGER, INTENT(in) :: i,j,n_var
            INTEGER :: function_idx

            function_idx = (n_var-1) * (Nx + 1) * (Ny + 1) + j * (Nx + 1) + i
            RETURN
      END FUNCTION function_idx

      SUBROUTINE compute_fdtd(cn)
            CLASS(cnfdtd), INTENT(inout) :: cn
            ! Variables locales
            CHARACTER(LEN=500) :: charac 
            LOGICAL :: display_it
            INTEGER :: info
            INTEGER :: n, m, n_var, n_node
            INTEGER :: i,j, idx,i_var
            INTEGER :: i0,j0,i1,j1
            INTEGER :: snapshot
            INTEGER, ALLOCATABLE :: A(:, :)

            ! Variable d'indices pour 3 variables inconnus Ex, Ey, Hz de dimension (Nx+1)*(Ny+1) chacune
            n_var = 3
            

            DO i_var = 1, n_var
                  DO j = 0, Ny
                        DO i = 0, Nx
                              idx = function_idx(i,j,i_var)
                              WRITE(*,'(3(AX,I5))') 'i=',i,' j=',j,' idx=',idx
                        END DO  
                  END DO
            END DO

            
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