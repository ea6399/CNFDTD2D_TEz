PROGRAM main

      ! Module
      USE numerics
      USE source
      USE fdtd
      USE test

      ! Variables
      IMPLICIT NONE
            type(cnfdtd) :: cn


      ! Début du programme
      CALL cn%init()

      

      ! Initialisation de la source
      CALL init_source()
      CALL compute_gauss(Esrc, base_Esrc, cn%dt)

      ! Calcul de la FDTD par CN
      CALL cn%compute_fdtd()

      ! Libération de la mémoire
      CALL cn%freememory()
      CALL free_source()


            ! Module de test
      !CALL init_test()





END PROGRAM main