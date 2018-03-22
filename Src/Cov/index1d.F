      integer FUNCTION index1d(i,j,L)

      implicit none
      integer i,j,L

      index1d=j*(L+1)+i+1

      RETURN
      END
