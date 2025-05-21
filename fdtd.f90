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
                  REAL(8), ALLOCATABLE :: B(:,:)
                  REAL(8), ALLOCATABLE :: A(:,:)
                  REAL(8), ALLOCATABLE :: c_E(:,:), c_H(:,:)
                  REAL(8) :: a1, a2
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
            ALLOCATE(cn%Ex  (                    0:Nx, 0:Ny                   ) )       
            ALLOCATE(cn%Ey  (                    0:Nx, 0:Ny                   ) )  
            ALLOCATE(cn%B   (   0 : 3 * (Nx + 1) - 1 , 0 : 3 * (Ny + 1) - 1   ) )
            ALLOCATE(cn%Hz  (                 0 : Nx , 0:Ny                   ) ) 
            ALLOCATE(cn%A   (   0 : 3 * (Nx + 1) - 1 , 0: 3 * (Ny + 1) - 1    ) )
            ALLOCATE(cn%c_E (                    0:Nx, 0:Ny                   ) )
            ALLOCATE(cn%c_H (                    0:Nx, 0:Ny                   ) )

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
            WRITE(*, '(/,T5,A,ES17.3, /)') 'a1 = ', cn%a1
            WRITE(*, '(/,T5,A,ES17.3, /)') 'a1/dx = ', cn%a1 / cn%dx
            WRITE(*, '(/,T5,A,ES17.3, /)') 'a2 = ', cn%a2
            WRITE(*, '(/,T5,A,ES17.3, /)') 'a2/dx = ', cn%a2 / cn%dx

            ! Initialisation des champs
            cn%A = 0.d0
            cn%B = 0.d0
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
            INTEGER :: ipiv(0:3*Nx+2)
            INTEGER :: n, m, nrhs
            INTEGER :: i,j
            INTEGER :: i0,j0,i1,j1,i2,j2
            INTEGER :: snapshot
            REAL(8), ALLOCATABLE :: A1(:,:)
            REAL(8), ALLOCATABLE :: A2(:,:)
            REAL(8), ALLOCATABLE :: A3(:,:)

            ALLOCATE(A1(0:Nx, 0:Ny))
            ALLOCATE(A2(0:Nx, 0:Ny))
            ALLOCATE(A3(0:Nx, 0:Ny))

            A1 = 0.d0; A2 = 0.d0; A3 = 0.d0

            WRITE (*, '(/,T5,A,I5)') "Nx = ", Nx
            WRITE(*,'(/, T5, A, I5X, I5)') "shape(A) = ", shape(cn%A)


            m = 0       
            


            !-------------------------------------------------------------!
            !------------------ Ecriture de la matrice A -----------------!
            !-------------------------------------------------------------!

            !---------------------------------------------------!
            !------------------ Sous matrice 1 -----------------!
            !---------------------------------------------------!
                  A1(0,0) = 1.0d0
                  A1(0,1) = - cn%a1 / cn%dy

                  DO j = 1, Nx - 1
                        A1(j,j-1) = cn%a1 / cn%dy
                        A1(j,j) = 1.0d0
                        A1(j,j+1) = - cn%a1 / cn%dy
                  END DO

                  A1(Nx, Nx - 1) = cn%a1 / cn%dy
                  A1(Nx, Nx) = 1.0d0

            ! ! Affichage de la matrice A1
            ! WRITE(*, '(/, T5, A, /)') "Matrice A1 :"
            ! DO i = 0, Nx
            !       WRITE(*, '(I5,500F12.2)') i, A1(i,:)
            ! END DO
            



            !---------------------------------------------------!
            !------------------ Sous matrice 2 -----------------!
            !---------------------------------------------------!     
                  A2(0,0) = 1.0d0
                  A2(0,1) = -cn%a1 / cn%dx

                  DO i = 1, Nx - 1
                        A2(i,i-1) =   cn%a1 / cn%dy
                        A2(i,i) = 1.0d0
                        A2(i,i+1) = - cn%a1 / cn%dy
                  END DO

                  A2(Nx, Nx - 1) =  cn%a1 / cn%dy
                  A2(Nx, Nx) = 1.0d0


            ! ! Affichage de la matrice A2
            ! WRITE(*, '(/, T5, A, /)') "Matrice A2 :"
            ! DO i = 0, Nx
            !       WRITE(*, '(I5,500F12.2)') i + Nx + 1, A2(i,:)
            ! END DO




            !---------------------------------------------------!
            !------------------ Sous matrice 3 -----------------!
            !---------------------------------------------------!

                  A3(0,0) = 1.d0
                  A3(0,1) = - cn%a2 / cn%dx
                  A3(1,0) =   cn%a2 / cn%dy
                  A3(1,1) = 0.d0

                  DO i = 3, Nx-3, 3
                        j = i
                        A3(i-1 ,j-1) = 0.d0
                        A3(i-1 ,j)   = - cn%a2/cn%dx 
                        A3(i-1 ,j+1) = 0.d0

                        A3(i,j-1)   = cn%a2/cn%dy
                        A3(i,j)     = 1.d0
                        A3(i,j+1)   = - cn%a2/cn%dx

                        A3(i+1,j-1) = 0.d0
                        A3(i+1,j)   = cn%a2/cn%dx
                        A3(i+1,j+1) = 0.d0
                  END DO

                  A3(Nx, Nx) = 1.d0
                  A3(Nx, Nx - 1) =   cn%a2 / cn%dx
                  A3(Nx - 1, Nx) = - cn%a2 / cn%dy

            ! ! Affichage de la matrice A3
            ! WRITE(*, '(/, T5, A, /)') "Matrice A3 :"
            ! DO i = 0, Nx
            !       WRITE(*, '(I5,500F12.5)') i + 2 * (Nx + 1), A3(i,:)
            ! END DO
            
             
            !---------------------------------------------------------------!
            !------------------ Assemblage de la matrice A -----------------!
            !---------------------------------------------------------------!

            i0 = 0;        j0 = 0
            i1 = Nx + 1;   j1 = Ny + 1
            i2 = 2*(Nx+1); j2 = 2*(Ny+1)

            ! Collage des blocs
            cn%A(i0  :i0 + Nx, j0  :j0 + Ny)   = A1
            cn%A(i1  :i1 + Nx, j1  :j1 + Ny)   = A2
            cn%A(i2  :i2 + Nx, j2  :j2 + Ny)   = A3

            ! ! Affichage de la matrice A
            ! WRITE(*, '(/, T5, A, /)') "Matrice A :"
            ! DO i = 0, 3 * Nx + 2
            !       WRITE(*, '(I5,500F12.5)') i, cn%A(i,:)
            ! END DO
            ! !---------------------------------------------------!

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
            ! ! !---------------------------------------------------!




            ! -------------------------------------------------------------!
            ! ------------------ Décomposition LU de A --------------------!
            ! -------------------------------------------------------------!

            CALL DGETRF(size(cn%A,1),size(cn%A,2),cn%A, size(cn%A,1),ipiv, info)
            IF (info > 0) THEN
                  WRITE(*,'(/,T5,A,I0,A,I0,A,/)') 'A(',info,',',info,') is exactly zero. '
            ELSE
                  WRITE(*,'(T5,A,I0,A,/)') 'The ',info,'-th argument had an ilegal value.'
            END IF
            !STOP 'LU decomposition error'


            ! ! ! Affichage de la matrice A
            ! WRITE(*, '(/, T5, A, /)') "Matrice A :"
            ! DO i = 0, 3 * Nx + 2
            !       WRITE(*, '(I5,500F12.5)') i, cn%A(i,:)
            ! END DO
            ! ! !---------------------------------------------------!


            


            ! Ouverture du fichier de sortie
            OPEN(idfile , file = "data/Ex.txt", status = "replace", action = "write", form = "formatted")
            OPEN(idfile + 1 , file = "data/Hz.txt", status = "replace", action = "write", form = "formatted")
            
            !-------------------------------------------------------------!
            !------------------- Boucle temporelle -----------------------!
            !-------------------------------------------------------------!
            WRITE(*, '(/, T5, "Injection de la source en ", I5)') i_src
            WRITE(*, '(/, T5, A, /)') "Début de la boucle temporelle"
            snapshot = 5

            nrhs = size(cn%B,2)
            DO n = 0, Nt - 1

                  m = m + 1

                  IF (MOD(n,20*snapshot) == 0) THEN
                        WRITE(*, '(/, T5, "itération temporelle : ",I4)') n
                  END IF

                  cn%Hz(i_src,j_src) =  Esrc(n)

                  !-------------------------------------------------------------!
                  !------------------- Ecriture du vecteur B -------------------!
                  !-------------------------------------------------------------!
                  DO i = 0, Nx 
                        DO j =1, Ny-1
                              cn%B(i,j)             = cn%Ex(i,j) + cn%a1 /cn%dy * (cn%Hz(i,j+1) - cn%Hz(i,j-1))
                        END DO
                  END DO

                  DO j = 0, Ny
                        DO i = 1, Nx-1
                              cn%B(i + i1, j + j1) = cn%Ey(i,j) - cn%a1 /cn%dx * (cn%Hz(i+1,j) - cn%Hz(i-1,j))
                        END DO
                  END DO

                  DO i = 1, Nx - 1
                        DO j = 1, Ny - 1
                              cn%B(i + i2, j + j2) = cn%Hz(i,j) - cn%a2 / cn%dy * (cn%Ex(i,j+1) - cn%Ex(i,j-1)) &
                                                                + cn%a2 / cn%dx * (cn%Ey(i+1,j) - cn%Ey(i-1,j))
                        END DO
                  END DO

                  CALL DGETRS('N',size(cn%A,1),nrhs,cn%A,size(cn%A,1),ipiv,cn%B,size(cn%B,1),info)

                  ! Sauvegarde des champs 
                  cn%Ex = cn%B( 0:Nx      ,       0:Ny)
                  cn%Ey = cn%B( i1:i1 + Nx, j1:j1 + Ny)
                  cn%Hz = cn%B( i2:i2 + Nx, j2:j2 + Ny)

                  ! Ecriture dans le fichier 
                  IF (MOD(n,snapshot) == 0) THEN
                        WRITE(idfile + 1, *) cn%Hz
                        WRITE(idfile    , *) cn%Ex
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


END MODULE fdtd