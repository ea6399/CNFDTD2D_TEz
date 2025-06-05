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

        ! Dimension intérieur
        nrow = 2*(Nx-1)
        ncol = (Ny-1)
        nvec = nrow*ncol

        allocate(rhs_test(0 : nvec - 1))
        rhs_test = 0.d0
        
        DO i = 0, nvec - 1 
            rhs_test(i) = i
        END DO
        DO i = 0, nvec - 1
            print*,i,rhs_test(i)
        END DO

        ALLOCATE(rhs_mat(0:nrow-1,0:ncol-1))

        rhs_mat = reshape(rhs_test, shape = [2 * (Nx - 1),Ny-1], order = [2,1])

        CALL display_matrix(rhs_mat,"rhs side reshaped")

        DO i = 0, nrow
            PRINT *, i, rhs_mat(i,:)
        END DO




    ENDSUBROUTINE init_test

END MODULE test