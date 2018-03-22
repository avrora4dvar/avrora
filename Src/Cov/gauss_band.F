      SUBROUTINE GAUSS_BAND(A,b,N,p,q)

! A. Kurapov, 1/8/2009: Solve matix eqn. (of size N) A x=b, where
! A is the banded matrix, p/q is is the lower/upper band width.
! b is overwritten with x

      implicit none 

      integer :: N, p, q

      real(8) :: A(p+q+1,N)
      real(8) :: b(N)

      integer :: k,ii,q1,ij,ik,kj,i,j
      
!===================================
! BAND GAUSS ELIMINATION: A=LU
!===================================
      q1=q+1
      DO k=1,N-1

       DO i=k+1,MIN(k+p,N)
        ii=i-k+q1
        A(ii,k)=A(ii,k)/A(q1,k)
       END DO
       
       DO j=k+1,MIN(k+q,N)
        DO i=k+1,MIN(k+p,N)     
         ij=i-j+q1
         ik=i-k+q1
         kj=k-j+q1
         A(ij,j)=A(ij,j)-A(ik,k)*A(kj,j)
        END DO
       END DO

      END DO      


!===================================
! BAND FWD SUBSTITUTION
!===================================
      DO j=1,N
       DO i=j+1,MIN(j+p,N)
        b(i)=b(i)-A(i-j+q1,j)*b(j)
       END DO
      END DO

!===================================
! BAND BACK SUBSTITUTION
!===================================
      DO j=N,1,-1
       b(j)=b(j)/A(q1,j)
       DO i=MAX(1,j-q),j-1
        b(i)=b(i)-A(i-j+q1,j)*b(j)
       END DO
      END DO

      RETURN
      END
