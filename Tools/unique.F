      integer FUNCTION nunique(typeo,K)

      implicit none
      integer :: K,k1,k2,n
      integer :: typeo(K)
      integer, allocatable :: tmp(:)
 
      n=0
      allocate(tmp(K))

      DO k1=1,K
       tmp(k1)=typeo(k1)
      END DO

      DO k1=1,K
       if (tmp(k1).ne.0) then
         n=n+1
         DO k2=k1+1,K
          if (tmp(k2).eq.tmp(k1)) then
           tmp(k2)=0      
          end if
         END DO
        end if
       END DO 

      deallocate(tmp)

      nunique=n

      RETURN
      END

      SUBROUTINE unique(typeo,K,uni,n)

      implicit none
      integer :: K,k1,k2,n,nn
      integer :: typeo(K)
      integer :: uni(n)
      integer, allocatable :: tmp(:)
 
      allocate(tmp(K))

      DO k1=1,K
       tmp(k1)=typeo(k1)
      END DO

      nn=0

      DO k1=1,K
       if (tmp(k1).ne.0) then
         nn=nn+1
         uni(nn)=tmp(k1)
         DO k2=k1+1,K
          if (tmp(k2).eq.tmp(k1)) then
           tmp(k2)=0      
          end if
         END DO
        end if
       END DO 

      deallocate(tmp)

      RETURN
      END
