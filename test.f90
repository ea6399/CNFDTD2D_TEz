! Moduke de test 
MODULE test
    USE numerics
    USE source
    USE fdtd
    
    IMPLICIT NONE

    ! Variables de test
    REAL(8) :: Atest(0:3,0:3), Btest(0:5,0:5), Ctest(0:10,0:10)
    REAL(8), ALLOCATABLE, DIMENSION(:,:) :: Aint_test, Bint_test, Cint_test
    REAL(8), ALLOCATABLE, dimension(:) :: rhs_test
    REAL(8), ALLOCATABLE, dimension(:,:) :: rhs_mat
    INTEGER :: i, nrow, ncol, nvec

    contains

    SUBROUTINE init_test()





    ENDSUBROUTINE init_test

END MODULE test