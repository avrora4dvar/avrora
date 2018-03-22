      real(8) function interp1(x,f,x1,n)

      implicit none
      integer i,n
      real(8) :: x1
      real(8) :: x(n),f(n) 

      IF (x1.lt.x(1)) THEN
       STOP 'interp1: x1<x(1)'
      END IF
      
      IF (x1.gt.x(n)) THEN
       STOP 'interp1: x1>x(n)'
      END IF
      
      i=1
      DO WHILE (x1.gt.x(i+1))
       i=i+1
      END DO

      interp1=f(i)+(f(i+1)-f(i))/(x(i+1)-x(i))*(x1-x(i))

      return
      end
   
