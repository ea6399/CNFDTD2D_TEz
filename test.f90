! Moduke de test 
MODULE test
    USE numerics
    USE source
    USE fdtd
    
    IMPLICIT NONE

    ! Variables de test
    REAL(8) :: Atest(0:3,0:3), Btest(0:5,0:5), Ctest(0:10,0:10)
    REAL(8), ALLOCATABLE, DIMENSION(:,:) :: Aint_test, Bint_test, Cint_test

    contains

    SUBROUTINE init_test()


        Atest = 1.d0
        Atest(0,:) = 0.0d0
        Atest(:,0) = 0.0d0
        Atest(3,:) = 0.0d0
        Atest(:,3) = 0.0d0

        CALL display_matrix(Atest, "Atest")
        CALL extract_matrix(Aint_test,Atest)
        CALL display_matrix(Aint_test)

        Btest = 1.d0
        Btest(0,:) = 0.0d0
        Btest(:,0) = 0.0d0
        Btest(5,:) = 0.0d0
        Btest(:,5) = 0.0d0

        CALL display_matrix(Btest, "Btest")
        CALL extract_matrix(Bint_test,Btest)
        CALL display_matrix(Bint_test, "B extracted")

        Ctest = 1.d0
        Ctest(0,:) = 0.0d0
        Ctest(:,0) = 0.0d0
        Ctest(10,:) = 0.0d0
        Ctest(:,10) = 0.0d0

        CALL display_matrix(Ctest, "Ctest")
        CALL extract_matrix(Cint_test, Ctest)
        CALL display_matrix(Cint_test, "B extracted")



    ENDSUBROUTINE init_test

END MODULE test