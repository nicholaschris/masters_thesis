module FortranModule
implicit none
contains



!*******************************************************************************
subroutine RBF_Function_Matrix(r,nKnown,KnownDispl,KnownCoords,nMovingDOF,MovingCoords,r_vector,mat)
    !------------------------------------------------------------------------------
    implicit none
    integer                                 :: i, j, k
    integer,intent(in)                      :: nKnown, nMovingDOF
    double precision,intent(in)                 :: KnownDispl(nKnown,1), KnownCoords(nKnown,3), MovingCoords(nMovingDOF,3)
    double precision                            :: xb(nKnown,3),Mbb(nKnown,nKnown), r, x_val, Pb(nKnown,4)
    double precision                            :: db(nKnown,1), zeros(4,4), zeros2(4,3)
    double precision                            :: result_m(nKnown+4,3)
    double precision, intent(out)               ::  r_vector(nKnown+4),mat(nKnown+4,nKnown+4)
    double precision,allocatable                :: work(:)
    integer,allocatable                     :: iwork(:)
    integer                                 :: lda,n, ind
    double precision                            :: alpha(nKnown,3),beta(4,3), x,y,z, sumx, sumy, sumz, px(4)
!     double precision,intent(out)                :: s_out(nMovingDOF,1)
    double precision                            :: s(nMovingDOF,1)
    double precision                            :: Kdispl(nknown,1), kcoords(nknown,3), mcoords(nMovingDOF,3)
!     real(kind=8)                            :: mat2D(nKnown+3,nKnown+3), r_vector2D(nKnown+3,2)
!     real(kind=8)                            :: temp_mat(nKnown+4,nKnown+4), temp_mat2(nKnown+4,nKnown+4)
    
    
!     real(kind=8)                            :: RCOND, z_solver(nKnown+4)
    integer                                 :: ipvt(nKnown+4)
    !------------------------------------------------------------------------------
Kdispl = knowndispl
Kcoords = knowncoords
Mcoords = movingcoords
allocate(work(nKnown+4),iwork(nknown+4))
! s_out = 0.0


xb = KCoords 
    !setup Mbb matrix
    Mbb = 0d0 
!$OMP PARALLEL DEFAULT(SHARED) private(i,j,x_val)
!$OMP DO
    do i=1,nKnown
        do j=1,nKnown
            x_val = dsqrt((xb(i,1)-xb(j,1))**2 + (xb(i,2)-xb(j,2))**2 + (xb(i,3)-xb(j,3))**2)
!             x_val = dsqrt((xb(i,1)-xb(j,1))**2 + (xb(i,2)-xb(j,2))**2 )
!             write(*,*) x_val,'x_val'
            if (x_val .le. r) then
                Mbb(i,j) = ((1.d0-x_val/r)**4)*(4.d0*(x_val/r)+1.d0)
            end if
            
        end do
    end do
!$omp end do
!$omp end parallel

    
    !setup Pb matrix
    do i = 1,nKnown
        Pb(i,1:4) = [1d0, xb(i,1), xb(i,2), xb(i,3)]
!         Pb(i,1:2) = [1d0,
    end do 
    
    
    
    !define known boundary movement
    db = 0d0
    db = KDispl
    
    !NOTE: this is done because it is actually 2D
!     db(:,3) = db(:,3)
    
    !--------Define solution matrix---------------------------------------------
!     allocate(mat(nbound+3,nbound+3),zeros(4,3))
    
    mat = 0d0
    zeros = 0d0
    mat(1:nKnown,1:nKnown) = Mbb(:,:)
    mat(1:nKnown,nKnown+1:nKnown+4) = Pb(1:nKnown,1:4)
    
    mat(nKnown+1:nKnown+4,1:nKnown) = transpose(Pb)
    mat(nKnown+1:nKnown+4,nKnown+1:nKnown+4) = zeros
    
!     mat2D = mat(1:nKnown+3,1:nKnown+3)
    
    
    !-------Define solution vector--------------------------------------------
    zeros2 = 0d0

    r_vector(1:nKnown) = db(:,1)
    r_vector(nKnown+1:nKnown+4) = 0d0

end subroutine RBF_Function_Matrix
!*******************************************************************************




!*******************************************************************************
subroutine RBF_Evaluate(r,nKnown,KnownDispl,KnownCoords,nMovingDOF,MovingCoords,s_out,result_m)
    !------------------------------------------------------------------------------
    implicit none
    integer                                 :: i, j, k, n_increments
    integer,intent(in)                      :: nKnown, nMovingDOF
    double precision,intent(in)                 :: KnownDispl(nKnown,1), KnownCoords(nKnown,3), MovingCoords(nMovingDOF,3)
    double precision,intent(in)                 :: result_m(nKnown+4)
    double precision                            :: xb(nKnown,3),Mbb(nKnown,nKnown), r, x_val, Pb(nKnown,4)
    double precision                            :: db(nKnown,1),mat(nKnown+4,nKnown+4), zeros(4,4), zeros2(4,3)
    double precision                            :: r_vector(nKnown+4,3)
    double precision,allocatable                :: work(:)
    integer,allocatable                     :: iwork(:)
    integer                                 :: lda,n, ind
    double precision                            :: alpha(nKnown,3),beta(4,3), x,y,z, sumx, sumy, sumz, px(4)
    double precision,intent(out)                :: s_out(nMovingDOF,1)
    double precision                            :: s(nMovingDOF,1)
    double precision                            :: Kdispl(nknown,1), kcoords(nknown,3), mcoords(nMovingDOF,3)
!     real(kind=8)                            :: mat2D(nKnown+3,nKnown+3), r_vector2D(nKnown+3,2)
!     real(kind=8)                            :: temp_mat(nKnown+4,nKnown+4), temp_mat2(nKnown+4,nKnown+4)
    
    
!     real(kind=8)                            :: RCOND, z_solver(nKnown+4)
    integer                                 :: ipvt(nKnown+4)
    !------------------------------------------------------------------------------
    Kdispl = knowndispl
Kcoords = knowncoords
Mcoords = movingcoords
allocate(work(nKnown+4),iwork(nknown+4))
s_out = 0.0
xb = KCoords 
  !--------Compute alpha and beta---------------------------------------------
!     allocate(alpha(nbound,2),beta(3,2))
    alpha(1:nKnown,1) = result_m(1:nKnown)
    beta(1:4,1) = result_m(nKnown+1:nKnown+4)
    s = 0d0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I,J,SUMX,SUMY,sumz,X,Y,z,X_VAL,PX)
!$OMP DO
    do i = 1,nMovingDOF
        x = MCoords(i,1)
        y = MCoords(i,2)
        z = MCoords(i,3)
        
        sumx=0d0; sumy = 0d0; sumz = 0d0;
        do j = 1,nKnown
            x_val = dsqrt((x-xb(j,1))**2 + (y-xb(j,2))**2 + (z-xb(j,3))**2)
!             x_val = dsqrt((x-xb(j,1))**2 + (y-xb(j,2))**2)
!             write(*,*) x_val
            if (x_val .le. r) then
!                 sumx = sumx + alpha(j,1)*(1+x_val**2)
!                 sumy = sumy + alpha(j,2)*(1+x_val**2)
!                 sumz = sumz + alpha(j,3)*(1+x_val**2)
                sumx = sumx + alpha(j,1)*((1.d0-x_val/r)**4)*(4.d0*(x_val/r)+1.d0)
!                 sumy = sumy + alpha(j,2)*((1.d0-x_val/r)**4)*(4.d0*(x_val/r)+1.d0)
!                 sumz = sumz + alpha(j,3)*((1.d0-x_val/r)**4)*(4.d0*(x_val/r)+1.d0)
            end if
        end do
        
        px = [1d0, x, y,z]

        s(i,1) = sumx
!         s(i,2) = sumy
!         s(i,3) = sumz
        
        do j = 1,4
            s(i,1) = s(i,1) + beta(j,1)*px(j)
!             s(i,2) = s(i,2) + beta(j,2)*px(j)
!             s(i,3) = s(i,3) + beta(j,3)*px(j)
        end do
        
        
    end do
!$OMP END DO NOWAIT
!$OMP END PARALLEL 

s_out = s
! write(*,*) 's from rbf pressure: ',s

end subroutine RBF_Evaluate
!*******************************************************************************

end module FortranModule