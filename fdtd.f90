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

      CONTAINS

      SUBROUTINE init(cn)
            CLASS(cnfdtd), INTENT(inout) :: cn
            INTEGER :: i

            ! Initialisation des variables
            ALLOCATE(cn%N_d (        0:10       ) )                                      ! Grid sampling densities
            ALLOCATE(cn%S   (        0:50       ) )                                      ! Courant Number 
            ALLOCATE(cn%Ex  (                    0:Nx, 0:Ny                               ) )       
            ALLOCATE(cn%Ey  (                    0:Nx, 0:Ny                               ) )  
            ALLOCATE(cn%B   (                      0 : 2 * (Nx + 1) * (Ny + 1) - 1        ) )
            ALLOCATE(cn%J   (                0 : Nx  ,  0 : Ny                            ) )
            ALLOCATE(cn%Hz  (                 0 : Nx , 0:Ny                               ) ) 
            ALLOCATE(cn%A   (   0 : 2 * (Nx + 1) - 1 , 0: 2 * (Ny + 1) - 1                ) )
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
            INTEGER :: info
            INTEGER :: n, m, nrhs
            INTEGER :: i,j
            INTEGER :: i0,j0,i1,j1
            INTEGER :: snapshot
            INTEGER :: idx_Ex, idx_Ey
            REAL(8), ALLOCATABLE :: A1(:,:)
            REAL(8), ALLOCATABLE :: A2(:,:)
            REAL(8), ALLOCATABLE :: A3(:,:)
            REAL(8), ALLOCATABLE :: A4(:,:)
            REAL(8), ALLOCATABLE :: BB(:,:)

            ALLOCATE(A1(0:Nx, 0:Ny))
            ALLOCATE(A2(0:Nx, 0:Ny))
            ALLOCATE(A3(0:Nx, 0:Ny))
            ALLOCATE(A4(0:Nx, 0:Ny))
            ALLOCATE(BB(0:2 * Nx + 1, 0:Ny))

            A1 = 0.d0; A2 = 0.d0; A3 = 0.d0; A4 = 0.d0; BB = 0.d0

            WRITE(*,'(/,T5,A,I5)') "Nx = ", Nx
            WRITE(*,'(/, T5, A, I5X, I5)') "shape(A) = ", shape(cn%A)
            WRITE(*,'(/, T5, A, I15X)')    "shape(B) = ", shape(cn%B)
            WRITE(*,'(/, T5, A, I5X, I5)') "shape(BB) = ", shape(BB)
            WRITE(*,'(/, T5, A, I5X, I5)') "shape(Ex) = ", shape(cn%Ex)
            WRITE(*,'(/, T5, A, I5X, I5)') "shape(Ey) = ", shape(cn%Ey)
            WRITE(*,'(/, T5, A, I5X, I5)') "shape(Hz) = ", shape(cn%Hz)
            WRITE(*,'(/, T5, A, I5X, I5)') "shape(J) = ", shape(cn%J)

            m = 0       
            


            !-------------------------------------------------------------!
            !------------------ Ecriture de la matrice A -----------------!
            !-------------------------------------------------------------!

            !----------------------------------------------------!
            !------------------ Sous matrice A1 -----------------!
            !----------------------------------------------------!
                  A1(0,0) = 1.0d0 + 2.d0 * cn%bx**2
                  !A1(0,1) = - cn%bx**2 

                  DO j = 1, Nx
                        A1(j,j-1) = - cn%bx**2
                        A1(j,j) = 1.0d0 + 2.d0 * cn%bx**2
                        !A1(j,j+1) = - cn%bx**2
                  END DO

            ! ! Affichage de la matrice A1
            ! WRITE(*, '(/, T5, A, /)') "Matrice A1 :"
            ! DO i = 0, Nx
            !       WRITE(*, '(I5,500F12.2)') i, A1(i,:)
            ! END DO
            



            !----------------------------------------------------!
            !------------------ Sous matrice A2 -----------------!
            !----------------------------------------------------!     
                  A2(0,0) = 1.0d0 + 2.d0 * cn%by**2
                  !A2(0,1) = - cn%by**2

                  DO i = 1, Nx - 1
                        A2(i,i-1) = - cn%by**2 
                        A2(i,i)   = 1.0d0 + 2.d0 * cn%by**2
                        !A2(i,i+1) = - cn%by**2
                  END DO

                  A2(Nx, Nx-1) =  - cn%by**2
                  A2(Nx, Nx) = 1.0d0 + 2.d0 * cn%by**2


            ! ! Affichage de la matrice A2
            ! WRITE(*, '(/, T5, A, /)') "Matrice A2 :"
            ! DO i = 0, Nx
            !       WRITE(*, '(I5,500F12.2)') i + Nx + 1, A2(i,:)
            ! END DO




            !----------------------------------------------------!
            !------------------ Sous matrice A3 -----------------!
            !----------------------------------------------------!

                  A3(0,0) =  1.d0
                  A3(2,0) = -1.d0
                  A3(1,1) =  2.d0

                  DO i = 1, Nx-2
                        A3(i , i)   =   2.d0
                        A3(i+2 , i) =  -1.d0
                  END DO 

                  A3(Nx-1,Nx-1) =  2.d0
                  A3(Nx, Nx - 2)= -1.d0   
                  A3(Nx, Nx)    =  1.d0

                  A3 = cn%bx * cn%by * A3

            ! Affichage de la matrice A3
            WRITE(*, '(/, T5, A, /)') "Matrice A3 :"
            DO i = 0, Nx
                  WRITE(*, '(I5,500F12.2)') i, A3(i,:)
            END DO

            !----------------------------------------------------!
            !------------------ Sous matrice A4 -----------------!
            !----------------------------------------------------!

                  A4 = -transpose(A3)

                  !A3 = 0.0d0

            ! ! Affichage de la matrice A4
            ! WRITE(*, '(/, T5, A, /)') "Matrice A4 :"
            ! DO i = 0, Nx
            !       WRITE(*, '(I5,500F12.2)') i + 3 * (Nx + 1), A4(i,:)
            ! END DO
            ! ! !---------------------------------------------------!

                  
            
             
            !---------------------------------------------------------------!
            !------------------ Assemblage de la matrice A -----------------!
            !---------------------------------------------------------------!


            ! Détermine les indices de collage
            i0 = 0;        j0 = 0
            i1 = Nx + 1;   j1 = Ny + 1
            

            WRITE(*, '(/, T5, A, I5, I5)') "i0, j0 = ", i0, j0
            WRITE(*, '(/, T5, A, I5, I5)') "i1, j1 = ", i1, j1

            ! Collage des blocs diagonaux
            cn%A(i0  :i0 + Nx, j0  :j0 + Ny)   = A1
            cn%A(i1  :i1 + Nx, j1  :j1 + Ny)   = A2
            ! Collage des matrices de couplage
            cn%A(i0 : i0 + Nx, j1 : j1 + Ny) = A3
            cn%A(i1 : i1 + Nx, j0 : j0 + Ny) = A4




            ! ! Affichage de la matrice A
            ! WRITE(*, '(/, T5, A, /)') "Matrice A :"
            ! DO i = 0, 2 * Nx + 1
            !       WRITE(*, '(I5,500F12.2)') i, cn%A(i,:)
            ! END DO
            ! !---------------------------------------------------!

            ! Vérification de la symétrie de la matrice A
            !CALL matrix_sym(cn%A)
            ! ! !---------------------------------------------------!

            ! Libération de mémoire 
            IF (ALLOCATED(A1)) THEN
                  DEALLOCATE(A1)
            END IF
            IF (ALLOCATED(A2)) THEN
                  DEALLOCATE(A2)
            END IF
            IF (ALLOCATED(A3)) THEN
                  DEALLOCATE(A3)
            END IF
            IF (ALLOCATED(A4)) THEN
                  DEALLOCATE(A4)
            END IF
            ! ! !---------------------------------------------------!




            ! -------------------------------------------------------------------!
            ! ------------------ Décomposition Cholesky de A --------------------!
            ! -------------------------------------------------------------------!

            CALL DPOTRF('L',size(cn%A,1),cn%A, size(cn%A,1), info)
            IF (info > 0) THEN
                  WRITE(*,'(/,T5,A,I0,A,/)') ' the leading proincipal minor of order ', info, 'is not positive.'
                  STOP 'Choletsky failed'
            ELSE IF (info < 0) THEN
                  WRITE(*,'(T5,A,I0,A,/)') 'The ',info,'-th argument had an ilegal value.'
                  STOP 'Choletsky failed'
            END IF
            

            


            ! Ouverture du fichier de sortie
            OPEN(idfile , file = "data/Ex.txt", status = "replace", action = "write", form = "formatted")
            OPEN(idfile + 1 , file = "data/Hz.txt", status = "replace", action = "write", form = "formatted")
            
            !-------------------------------------------------------------!
            !------------------- Boucle temporelle -----------------------!
            !-------------------------------------------------------------!
            WRITE(*, '(/, T5, "Injection de la source en ", I5)') i_src, j_src
            WRITE(*, '(/, T5, A, /)') "Début de la boucle temporelle"
            snapshot = 20

            nrhs = 1
            m = 0
            DO n = 0, Nt - 1

                  IF (MOD(n,5*snapshot) == 0) THEN
                        WRITE(*, '(/, T5, "itération temporelle : ",I4)') n
                  END IF

                  cn%J(i_src, j_src) =  Esrc(n)

                  !--------------------------------------------------------------!
                  !------------------- Ecriture du vecteur B --------------------!
                  !--------------------------------------------------------------!
                  ! On enregistre les résultats du temps précédent
                  cn%Ex = BB(0 : Nx, 0 : Ny)
                  cn%Ey = BB(i1 : i1 + Nx, 0 : Ny)

                 
                  

                  ! Second membre Ex
                  DO i = 0,  Nx
                        !print *, "i = ", i
                        DO j = 0, Ny
                              ! Détermine le bonne indice
                              idx_Ex = (i) * (Nx + 1) + (j)
                              !print *, "idx_Ex = ", idx_Ex
                              if (i == 0 .OR. i == Nx .OR. j == 0 .OR. j == Ny) then 
                                    cn%B(idx_Ex) = 0.d0
                              else
                                    cn%B(idx_Ex) =      (1.d0 - 2.d0 * cn%bx**2) * cn%Ex(i,j)                         & 
                                                + cn%bx**2 * ( cn%Ex(i, j - 1) + cn%Ex(i, j + 1) )                    &
                                                - cn%bx*cn%by * ( cn%Ey(i + 1, j + 1) - cn%Ey(i - 1, j + 1) )         &
                                                + cn%bx*cn%by * ( cn%Ey(i + 1 , j -1) - cn%Ey(i - 1, j - 1) )         &
                                                + 2.d0 * cn%a1 * (cn%Hz(i,j+1) - cn%Hz(i, j-1))
                              endif
                        END DO
                  END DO



                  ! Second membre Ey
                  DO i = 0 , Nx 
                        !print *, "i = ", i
                        DO j = 0,  Ny 
                              ! Détermine le bonne indice
                              idx_Ey = (Nx+1)*(Ny+1) + (i) * (Nx + 1) + j
                              !print *, "idx_Ey = ", idx_Ey
                              if (i == 0 .OR. i == Nx .OR. j == 0 .OR. j == Ny) then 
                                    cn%B(idx_Ey) = 0.d0
                              else
                              ! Calcul du second membre Ey
                                    cn%B(idx_Ey) =      (1.d0 - 2.d0 * cn%by**2)*cn%Ey(i,j)                              & 
                                                + cn%by**2 * ( cn%Ey(i - 1, j) + cn%Ey(i + 1, j)    )                    &
                                                - cn%bx*cn%by * ( cn%Ex(i + 1 , j + 1) - cn%Ex(i + 1 , j - 1)  )         &
                                                + cn%bx*cn%by * ( cn%Ex(i-1, j + 1)    - cn%Ex(i-1, j) )                 &
                                                - 2.d0 * cn%a1 * (cn%Hz(i+1,j) - cn%Hz(i-1, j))
                              END IF
                        END DO
                  END DO

                  

                  ! Résolution du système linéaire
                  CALL DPOTRS('L',size(cn%A,1),nrhs,cn%A,size(cn%A,1),cn%B,size(cn%B,1),info)

                  ! reshape du vecteur B / order = [2,1] fait varier j avant i
                  BB = reshape(cn%B, shape = [ 2 * Nx + 2, Ny + 1], order = [2, 1])

                  !  print *, "shape(BB) = ", shape(BB)

                  ! WRITE(*,'(2(AX,F16.10))') 'B(0)=',cn%B(0),' BB(0,0)=',BB(0,0)
                  ! WRITE(*,'(2(AX,F16.10))') 'B(1)=',cn%B(1),' BB(0,1)=',BB(0,1)
                  ! WRITE(*,'(2(AX,F16.10))') 'B(Ny+1)=',cn%B(Ny+1),' BB(1,0)=',BB(1,0)

                  

                  ! Injection de la source
                  BB = BB + cn%J




                  ! Mise à jour explicite de Hz
                  DO i = 0, Nx-1
                        DO j = 0, Ny-1
                              cn%Hz(i,j) = cn%Hz(i,j) + cn%a2 / cn%dy * ( BB(i,j + 1) - BB(i,j)               &
                                                                        + cn%Ex(i, j+ 1) - cn%Ex(i,j) )       &
                                                      - cn%a2 / cn%dx * ( BB(i1 + i + 1,j) - BB(i1 + i,j)     &
                                                                        + cn%Ey(i + 1, j) - cn%Ey(i,j) )
                        END DO
                  END DO

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


END MODULE fdtd