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
            ALLOCATE(cn%Ex  (     0:Nx, 0:Ny    ) )       
            ALLOCATE(cn%Ey  (     0:Nx, 0:Ny    ) )  
            ALLOCATE(cn%B   (        0:Nx       ) )
            ALLOCATE(cn%Hz  (       0:Nx, 0:Ny        ) ) 
            ALLOCATE(cn%A   ( 0 : 3 * Nx + 2, 0: 3 * Ny + 2  ))
            ALLOCATE(cn%c_E (     0:Nx, 0:Ny    ) )
            ALLOCATE(cn%c_H (     0:Nx, 0:Ny    ) )

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
            WRITE(*, '(/,T5,A,ES17.3, /)') 'a2 = ', cn%a2

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
            INTEGER :: n, m
            INTEGER :: i,j
            INTEGER :: i0,j0,i1,j1,i2,j2
            INTEGER :: snapshot
            REAL(8) :: A1(0:Nx, 0:Ny)
            REAL(8) :: A2(0:Nx, 0:Ny)
            REAL(8) :: A3(0:Nx, 0:Ny)

            A1 = 0.d0; A2 = 0.d0; A3 = 0.d0

            WRITE (*, '(/,T5,A,I5)') "Nx = ", Nx
            WRITE(*,'(/, T5, A, I3X, I3)') "shape(A) = ", shape(cn%A)


            m = 0       


            !-------------------------------------------------------------!
            !------------------ Ecriture de la matrice A -----------------!
            !-------------------------------------------------------------!

            !---------------------------------------------------!
            !------------------ Sous matrice 1 -----------------!
            !---------------------------------------------------!
                  A1(0,0) = 1.0d0
                  A1(0,1) = - cn%a1 / cn%dy

                  DO i = 1, Nx - 1
                        A1(i,i-1) = cn%a1 / cn%dy
                        A1(i,i) = 1.0d0
                        A1(i,i+1) = - cn%a1 / cn%dy
                  END DO

                  A1(Nx, Nx - 1) = cn%a1 / cn%dy
                  A1(Nx, Nx) = 1.0d0

            ! Affichage de la matrice A1
            WRITE(*, '(/, T5, A, /)') "Matrice A1 :"
            DO i = 0, Nx
                  WRITE(*, '(I5,500F12.2)') i, A1(i,:)
            END DO
            
            !---------------------------------------------------!
            !------------------ Sous matrice 2 -----------------!
            !---------------------------------------------------!     
                  
                  A2 = A1

            ! Affichage de la matrice A2
            WRITE(*, '(/, T5, A, /)') "Matrice A2 :"
            DO i = 0, Nx
                  WRITE(*, '(I5,500F12.2)') i, A2(i,:)
            END DO

            !---------------------------------------------------!
            !------------------ Sous matrice 3 -----------------!
            !---------------------------------------------------!

                  A3(0,0) = 1.d0
                  A3(0,1) = - cn%a2 / cn%dx
                  A3(1,0) = - cn%a2 / cn%dy
                  A3(1,1) = 0.d0

                  DO i = 3, Nx-3, 3
                        j = i
                        A3(i-1 ,j-1) = 0.d0
                        A3(i-1 ,j)   = cn%a2/cn%dx 
                        A3(i-1 ,j+1) = 0.d0

                        A3(i,j-1)   = - cn%a2/cn%dy
                        A3(i,j)     = 1.d0
                        A3(i,j+1)   = - cn%a2/cn%dx

                        A3(i+1,j-1) = 0.d0
                        A3(i+1,j)   = cn%a2/cn%dx
                        A3(i+1,j+1) = 0.d0
                  END DO

                  A3(Nx, Nx) = 1.d0
                  A3(Nx, Nx - 1) = - cn%a2 / cn%dx
                  A3(Nx - 1, Nx) =   cn%a2 / cn%dy

            ! Affichage de la matrice A3
            WRITE(*, '(/, T5, A, /)') "Matrice A3 :"
            DO i = 0, Nx
                  WRITE(*, '(I5,500F12.5)') i, A3(i,:)
            END DO
            !---------------------------------------------------!

            !------------------ Assemblage de la matrice A -----------------!
            !---------------------------------------------------------------!

            i0 = 0;        j0 = 0
            i1 = Nx + 1;   j1 = Ny + 1
            i2 = 2*(Nx+1); j2 = 2*(Ny+1)

            ! Collage des blocs
            cn%A(i0  :i0 + Nx, j0  :j0 + Ny)   = A1
            cn%A(i1  :i1 + Nx, j1  :j1 + Ny)   = A2
            cn%A(i2  :i2 + Nx, j2  :j2 + Ny)   = A3

            ! Affichage de la matrice A
            WRITE(*, '(/, T5, A, /)') "Matrice A :"
            DO i = 0, 3 * Nx + 2
                  WRITE(*, '(I5,500F12.5)') i, cn%A(i,:)
            END DO
            ! !---------------------------------------------------!


            


            ! Ouverture du fichier de sortie
            OPEN(idfile , file = "data/E.txt", status = "replace", action = "write", form = "formatted")
            OPEN(idfile + 1 , file = "data/H.txt", status = "replace", action = "write", form = "formatted")
            
            !-------------------------------------------------------------!
            !------------------- Boucle temporelle -----------------------!
            !-------------------------------------------------------------!
            WRITE(*, '(/, T5, "Injection de la source en ", I5)') i_src
            WRITE(*, '(/, T5, A, /)') "Début de la boucle temporelle"
            snapshot = 5
            DO n = 0, Nt - 1

                  IF (MOD(n,20*snapshot) == 0) THEN
                        WRITE(*, '(/, T5, "itération temporelle : ",I4)') n
                  END IF

                  !-----------------------------------------------------------------------------------------!
                  !------------------ Ecriture du vecteur second membre B_E et résolution ------------------!
                  !-----------------------------------------------------------------------------------------!

                  
                  !-----------------------------------------------------------------------------------------!
                  !------------------ Ecriture du vecteur second membre B_H et résolution ------------------!
                  !-----------------------------------------------------------------------------------------!

                  
                  
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