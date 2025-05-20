module source

      use numerics

      IMPLICIT NONE

      CONTAINS

            ! Fonction gaussienne
            FUNCTION gauss_t(n,dt)
                  REAL(8) :: gauss_t, dt

                  ! Arguments
                  INTEGER, intent(in) :: n

                  gauss_t = a0 * exp( - ( ( n * dt - t0 ) / T )**2 )
            ENDFUNCTION gauss_t

            SUBROUTINE compute_gauss(E,base, dt)
                  REAL(8), intent(inout) :: base(0:Nt - 1)
                  REAL(8), intent(inout) :: E(0:Nt - 1)
                  REAL(8), intent(in) :: dt

                  ! Variables locales
                  INTEGER :: n

                  ! Intervalles de définition de la gaussienne temporelle
                  DO n = 0, Nt - 1
                        base(n) = n * dt
                  end do

                  ! Calcul de la gaussienne
                  DO n = 0, Nt - 1
                        E(n) = gauss_t(n, dt)
                  END DO

                  OPEN(15, file = "data/Esrc.txt", status = "replace", action = "write", form = "formatted")
                        DO n = 0, Nt - 1
                              WRITE(15, *) base_Esrc(n), Esrc(n)
                        END DO
                  CLOSE(15)
            ENDSUBROUTINE compute_gauss


            SUBROUTINE free_source()
                  ! Libération de la mémoire
                  IF (ALLOCATED(Esrc)) THEN
                        DEALLOCATE(Esrc)
                  END IF
                  IF (ALLOCATED(base_Esrc)) THEN
                        DEALLOCATE(base_Esrc)
                  END IF
            END SUBROUTINE free_source

end module source