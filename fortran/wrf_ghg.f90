SUBROUTINE DCOMPUTEXGHG(sfc_p, pres, nx, ny, nz, ant, bio, bck, xghg)
    IMPLICIT NONE

    !f2py threadsafe
    !f2py intent(in,out) :: xghg

    INTEGER, INTENT(IN) :: nx, ny, nz
    REAL(KIND=8), DIMENSION(nx, ny), INTENT(IN) :: sfc_p
    REAL(KIND=8), DIMENSION(nx, ny, nz), INTENT(IN) :: pres
    REAL(KIND=8), DIMENSION(nx, ny, nz), INTENT(IN) :: ant
    REAL(KIND=8), DIMENSION(nx, ny, nz), INTENT(IN) :: bio
    REAL(KIND=8), DIMENSION(nx, ny, nz), INTENT(IN) :: bck
    REAL(KIND=8), DIMENSION(nx, ny), INTENT(IN,OUT) :: xghg

    INTEGER :: i, j, k
    REAL(KIND=8), DIMENSION(nx, ny, nz) :: ghg
    REAL(KIND=8), DIMENSION(nx, ny, nz) :: pres_bound
    REAL(KIND=8), DIMENSION(nx, ny, nz) :: p_layer_diff
    REAL(KIND=8), DIMENSION(nx, ny) :: p_diff



    !$OMP PARALLEL
    !$OMP DO COLLAPSE(3) SCHEDULE(runtime)
    DO k = 1,nz
        DO j = 1,ny
            DO i = 1,nx
                ghg(i , j, k) = ant(i, j, k) + bio(i, j, k) - bck(i, j, k)
                IF (k .eq. 1) THEN
                    pres_bound(i, j, k) = sfc_p(i, j)
                ELSE 
                    pres_bound(i, j, k) = pres_bound(i, j, k-1) + (2*(pres(i, j, k-1)-pres_bound(i, j, k-1)))
                END IF
            END DO
        END DO
    END DO
    !$OMP END DO
    
    !$OMP DO COLLAPSE(2) SCHEDULE(runtime)
    DO j = i,ny
        DO i = 1,nx
        END DO
    END DO

