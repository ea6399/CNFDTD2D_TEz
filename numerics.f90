MODULE numerics
      ! This module contains numerical methods and constants
      IMPLICIT NONE
      ! Constant Parameters
      REAL(8), PARAMETER :: pi = atan(1.0d0) * 4.0d0
      REAL(8), PARAMETER :: EPS = 1.0d-10                     ! Epsilon pour les comparaisons
      REAL(8), PARAMETER :: epsilon_0 = 8.854187817d-12                 ! Permittivité du vide
      REAL(8), PARAMETER :: epsilon_r = 4.0d0 * epsilon_0               ! Permittivité relative
      REAL(8), PARAMETER :: mu_0 = 1.256637061d-6                       ! Perméabilité du vide
      REAL(8), PARAMETER :: c = 1.0d0 / sqrt(epsilon_0 * mu_0)          ! Vitesse de la lumière

      ! Constants for the Gaussian pulse
      REAL(8), PARAMETER :: fmax = 1.0d9                                ! Fréquence max
      REAL(8), PARAMETER :: attfmax = 10.0d0                            ! Atténuation de la fréquence max
      REAL(8), PARAMETER :: att0 = 1000.0d0                             ! Atténuation de la fréquence 0
      REAL(8), PARAMETER :: a0 = 1.0d0
      REAL(8), PARAMETER :: T = sqrt( log(attfmax) ) / (PI * fmax)      ! Largeur de la gaussienne
      REAL(8), PARAMETER :: t0 = T * sqrt( log(att0) )                  ! Retard de la gaussienne

      ! Global Parameters
      INTEGER, PARAMETER :: Nt = 1001                                   ! Nombre d'échantillons temps
      INTEGER, PARAMETER :: Nx = 99                                    ! Nombre d'échantillons espace suivant x
      INTEGER, PARAMETER :: Ny = Nx                                    ! Nombre d'échantillons espace suivant y
      INTEGER, PARAMETER :: i_src = int(Nx / 2) + 1                    ! Injection de la source suivant l'axe x
      INTEGER, PARAMETER :: j_src = i_src                                 ! Injection de la source suivant l'axe y
      INTEGER, PARAMETER :: idfile = 50                                 ! idfile
      


      REAL(8), ALLOCATABLE :: Esrc(:), base_Esrc(:)

      CONTAINS

      SUBROUTINE init_source()
            IMPLICIT NONE 
            
            if( .NOT. ALLOCATED(Esrc) )          ALLOCATE(Esrc(0 : Nt - 1))
            if( .NOT. ALLOCATED(base_Esrc) )     ALLOCATE(base_Esrc(0 : Nt - 1))

            Esrc = 0.d0
            base_Esrc = 0.d0

      ENDSUBROUTINE init_source
            
END MODULE numerics