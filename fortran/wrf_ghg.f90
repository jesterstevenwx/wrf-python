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
    REAL(KIND=8), DIMENSION(nx, ny), INTENT(OUT) :: xghg

    INTEGER :: i, j, k
    REAL(KIND=8), DIMENSION(nx, ny, nz) :: ghg
    REAL(KIND=8), DIMENSION(nx, ny, nz) :: pres_bound
    REAL(KIND=8), DIMENSION(nx, ny, nz) :: p_layer_diff
    REAL(KIND=8), DIMENSION(nx, ny) :: p_diff
    REAL(KIND=8), DIMENSION(nx, ny) :: weighted_ghg



    !$OMP PARALLEL
    !$OMP DO COLLAPSE(3) SCHEDULE(runtime)
    DO k = 1,nz
        DO j = 1,ny
            DO i = 1,nx
                ghg(i , j, k) = ant(i, j, k) + bio(i, j, k) - bck(i, j, k)
                IF (k .eq. 1) THEN
                    pres_bound(i, j, k) = sfc_p(i, j)
                    p_layer_diff(i, j, k) = pres(i, j, k) - sfc_p(i, j)
                ELSE 
                    pres_bound(i, j, k) = pres_bound(i, j, k-1) + (2*(pres(i, j, k-1)-pres_bound(i, j, k-1)))
                    p_layer_diff(i, j, k) = pres_bound(i, j, k-1) - pres_bound(i, j, k)
                END IF
            END DO
        END DO
    END DO
    !$OMP END DO

    weighted_ghg = SUM(ghg * p_layer_diff, 3)

    !$OMP DO COLLAPSE(2) SCHEDULE(runtime)
    DO j = i,ny
        DO i = 1,nx
            p_diff(i, j) = pres_bound(i, j, 1) - pres_bound(i, j, nz)
            xghg(i, j) = weighted_ghg(i, j) / p_diff(i, j)
        END DO
    END DO
    !$OMP END DO
    !$OMP END PARALLEL

    RETURN
END SUBROUTINE DCOMPUTEXGHG