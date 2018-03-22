MODULE mod_constants

 IMPLICIT NONE

 REAL(8)::R_earth=DBLE(6371000)
 REAL(8)::pi=DBLE(3.141592653589793)

 INTERFACE array2mat
  MODULE PROCEDURE array2mat_2D,array2mat_3D,array2mat_4D,&
  &array2mat_2B,array2mat_3B,array2mat_4B
 END INTERFACE

 INTERFACE mat2array
  MODULE PROCEDURE mat2array_2D,mat2array_3D,mat2array_4D,&
  &mat2array_2B,mat2array_3B,mat2array_4B
 END INTERFACE

CONTAINS

 SUBROUTINE array2mat_2D(n,val_in,val_out)
 !ARRAY2MAT_2D Convert a array to a 2D matrix with array dimension n
 !as its first dimension
 !
 !Syntax:
 ! CALL array2mat(n,val_in,val_out)
 !Input:
 ! n: INTEGER with the dimension of the input array that has to become
 !    the 1st dimension of the output array
 ! val_in: input DOUBLE ARRAY
 !Output:
 ! val_out: 2D output DOUBLE ARRAY

  INTEGER,INTENT(in)::n
  REAL(8),INTENT(in)::val_in(:,:)
  REAL(8),ALLOCATABLE,INTENT(inout)::val_out(:,:)
  INTEGER::m_size(3)
  REAL(8)::pad(2)

  pad=DBLE(0)
  IF(n.GT.2) THEN
   WRITE(*,*) 'mod_constants.array2mat: n may not exceede dimension input&
   & array.'
   STOP
  END IF

  m_size(1)=MAXVAL([1,PRODUCT(SHAPE(val_in),[1,2].LT.n)])
  m_size(2)=MAXVAL([1,PRODUCT(SHAPE(val_in),[1,2].EQ.n)])
  m_size(3)=MAXVAL([1,PRODUCT(SHAPE(val_in),[1,2].GT.n)])

  IF(ALLOCATED(val_out)) THEN
   IF(ANY(SHAPE(val_out).NE.[m_size(2),m_size(1)*m_size(3)])) THEN
    DEALLOCATE(val_out)
   END IF
  END IF
  IF(.NOT.ALLOCATED(val_out)) THEN
    ALLOCATE(val_out(m_size(2),m_size(1)*m_size(3)))
  END IF  

  val_out=&
  &RESHAPE(&
  &RESHAPE(val_in,[m_size(2)*m_size(3),m_size(1)],pad,[2,1])&
  &,[m_size(2),m_size(1)*m_size(3)]) 

 END SUBROUTINE array2mat_2D

!---------------------------------------------------------------------------

 SUBROUTINE array2mat_3D(n,val_in,val_out)
 !ARRAY2MAT_3D As array2mat_2 but now for 3D input arrays.  

  INTEGER,INTENT(in)::n
  REAL(8),INTENT(in)::val_in(:,:,:)
  REAL(8),ALLOCATABLE,INTENT(inout)::val_out(:,:)
  INTEGER::m_size(3)
  REAL(8)::pad(2)

  pad=DBLE(0)
  IF(n.GT.3) THEN
   WRITE(*,*) 'mod_constants.array2mat: n may not exceede dimension input&
   & array.'
   STOP
  END IF

  m_size(1)=MAXVAL([1,PRODUCT(SHAPE(val_in),[1,2,3].LT.n)])
  m_size(2)=MAXVAL([1,PRODUCT(SHAPE(val_in),[1,2,3].EQ.n)])
  m_size(3)=MAXVAL([1,PRODUCT(SHAPE(val_in),[1,2,3].GT.n)])

  IF(ALLOCATED(val_out)) THEN
   IF(ANY(SHAPE(val_out).NE.[m_size(2),m_size(1)*m_size(3)])) THEN
    DEALLOCATE(val_out)
   END IF
  END IF
  IF(.NOT.ALLOCATED(val_out)) THEN
    ALLOCATE(val_out(m_size(2),m_size(1)*m_size(3)))
  END IF  

  val_out=&
  &RESHAPE(&
  &RESHAPE(val_in,[m_size(2)*m_size(3),m_size(1)],pad,[2,1])&
  &,[m_size(2),m_size(1)*m_size(3)]) 

 END SUBROUTINE array2mat_3D

!----------------------------------------------------------------------------

 SUBROUTINE array2mat_4D(n,val_in,val_out)
 !ARRAY2MAT_4D As array2mat_2D but now for 4D input arrays 

  INTEGER,INTENT(in)::n
  REAL(8),INTENT(in)::val_in(:,:,:,:)
  REAL(8),ALLOCATABLE,INTENT(inout)::val_out(:,:)
  INTEGER::m_size(3)
  REAL(8)::pad(2)

  pad=DBLE(0)
  IF(n.GT.4) THEN
   WRITE(*,*) 'mod_constants.array2mat: n may not exceede dimension input&
   & array.'
   STOP
  END IF

  m_size(1)=MAXVAL([1,PRODUCT(SHAPE(val_in),[1,2,3,4].LT.n)])
  m_size(2)=MAXVAL([1,PRODUCT(SHAPE(val_in),[1,2,3,4].EQ.n)])
  m_size(3)=MAXVAL([1,PRODUCT(SHAPE(val_in),[1,2,3,4].GT.n)])

  IF(ALLOCATED(val_out)) THEN
   IF(ANY(SHAPE(val_out).NE.[m_size(2),m_size(1)*m_size(3)])) THEN
    DEALLOCATE(val_out)
   END IF
  END IF
  IF(.NOT.ALLOCATED(val_out)) THEN
    ALLOCATE(val_out(m_size(2),m_size(1)*m_size(3)))
  END IF  

  val_out=&
  &RESHAPE(&
  &RESHAPE(val_in,[m_size(2)*m_size(3),m_size(1)],pad,[2,1])&
  &,[m_size(2),m_size(1)*m_size(3)]) 

 END SUBROUTINE array2mat_4D

!----------------------------------------------------------------------------

 SUBROUTINE mat2array_2D(n,val_in,val_out)
 !MAT2ARRAY_2D Convert a 2D array created with array2mat back to a multi-
 !dimensional array.
 !
 !Syntax:
 ! CALL mat2array(n,val_in,val_out)
 !Input:
 ! n: INTEGER. The 1st dimension of val_in is going to become the nth
 !    dimension of val_out
 ! val_in: 2D DOUBLE ARRAY with input values
 !Output:
 ! val_out: DOUBLE ARRAY with the same number of entries as val_in

  INTEGER,INTENT(in)::n
  REAL(8),ALLOCATABLE,INTENT(in)::val_in(:,:)
  REAL(8),ALLOCATABLE,INTENT(inout)::val_out(:,:)
  INTEGER::m_size(3)
  REAL(8)::pad(2)

  pad=DBLE(0)
  IF(n.GT.2) THEN
   WRITE(*,*) 'mod_constants.array2mat: n may not exceede dimension input&
   & array.'
   STOP
  END IF
  IF(SIZE(val_in).NE.SIZE(val_out)) THEN
   WRITE(*,*) 'mod_constants.mat2array: dimensions val_in and val_out are&
   & incongruent.'
   STOP
  END IF

  m_size(1)=MAXVAL([1,PRODUCT(SHAPE(val_out),[1,2].LT.n)])
  m_size(2)=MAXVAL([1,PRODUCT(SHAPE(val_out),[1,2].EQ.n)])
  m_size(3)=MAXVAL([1,PRODUCT(SHAPE(val_out),[1,2].GT.n)])

  val_out=RESHAPE(&
  &RESHAPE(val_in,[m_size(1),m_size(2)*m_size(3)],pad,[2,1])&
  &,SHAPE(val_out))  

 END SUBROUTINE mat2array_2D

!----------------------------------------------------------------------------

 SUBROUTINE mat2array_3D(n,val_in,val_out)
 !MAT2ARRAY_3D as mat2array_2D but now for 3D output array

  INTEGER,INTENT(in)::n
  REAL(8),ALLOCATABLE,INTENT(in)::val_in(:,:)
  REAL(8),ALLOCATABLE,INTENT(inout)::val_out(:,:,:)
  INTEGER::m_size(3)
  REAL(8)::pad(2)

  pad=DBLE(0)
  IF(n.GT.3) THEN
   WRITE(*,*) 'mod_constants.array2mat: n may not exceede dimension input&
   & array.'
   STOP
  END IF
  IF(SIZE(val_in).NE.SIZE(val_out)) THEN
   WRITE(*,*) 'mod_constants.mat2array: dimensions val_in and val_out are&
   & incongruent.'
   STOP
  END IF

  m_size(1)=MAXVAL([1,PRODUCT(SHAPE(val_out),[1,2,3].LT.n)])
  m_size(2)=MAXVAL([1,PRODUCT(SHAPE(val_out),[1,2,3].EQ.n)])
  m_size(3)=MAXVAL([1,PRODUCT(SHAPE(val_out),[1,2,3].GT.n)])

  val_out=RESHAPE(&
  &RESHAPE(val_in,[m_size(1),m_size(2)*m_size(3)],pad,[2,1])&
  &,SHAPE(val_out))  

 END SUBROUTINE mat2array_3D

!----------------------------------------------------------------------------

 SUBROUTINE mat2array_4D(n,val_in,val_out)
 !MAT2ARRAY_4D As mat2array_2D but now for 4D output array

  INTEGER,INTENT(in)::n
  REAL(8),ALLOCATABLE,INTENT(in)::val_in(:,:)
  REAL(8),ALLOCATABLE,INTENT(inout)::val_out(:,:,:,:)
  INTEGER::m_size(3)
  REAL(8)::pad(2)

  pad=DBLE(0)
  IF(n.GT.4) THEN
   WRITE(*,*) 'mod_constants.array2mat: n may not exceede dimension input&
   & array.'
   STOP
  END IF
  IF(SIZE(val_in).NE.SIZE(val_out)) THEN
   WRITE(*,*) 'mod_constants.mat2array: dimensions val_in and val_out are&
   & incongruent.'
   STOP
  END IF

  m_size(1)=MAXVAL([1,PRODUCT(SHAPE(val_out),[1,2,3,4].LT.n)])
  m_size(2)=MAXVAL([1,PRODUCT(SHAPE(val_out),[1,2,3,4].EQ.n)])
  m_size(3)=MAXVAL([1,PRODUCT(SHAPE(val_out),[1,2,3,4].GT.n)])

  val_out=RESHAPE(&
  &RESHAPE(val_in,[m_size(1),m_size(2)*m_size(3)],pad,[2,1])&
  &,SHAPE(val_out))  

 END SUBROUTINE mat2array_4D

!----------------------------------------------------------------------------

SUBROUTINE array2mat_2B(n,val_in,val_out)
 !ARRAY2MAT_2B As array2mat_2D but now for LOGICALS

  INTEGER,INTENT(in)::n
  LOGICAL,INTENT(in)::val_in(:,:)
  LOGICAL,ALLOCATABLE,INTENT(inout)::val_out(:,:)
  INTEGER::m_size(3)
  LOGICAL::pad(2)

  pad=.FALSE.
  IF(n.GT.2) THEN
   WRITE(*,*) 'mod_constants.array2mat: n may not exceede dimension input&
   & array.'
   STOP
  END IF

  m_size(1)=MAXVAL([1,PRODUCT(SHAPE(val_in),[1,2].LT.n)])
  m_size(2)=MAXVAL([1,PRODUCT(SHAPE(val_in),[1,2].EQ.n)])
  m_size(3)=MAXVAL([1,PRODUCT(SHAPE(val_in),[1,2].GT.n)])

  IF(ALLOCATED(val_out)) THEN
   IF(ANY(SHAPE(val_out).NE.[m_size(2),m_size(1)*m_size(3)])) THEN
    DEALLOCATE(val_out)
   END IF
  END IF
  IF(.NOT.ALLOCATED(val_out)) THEN
    ALLOCATE(val_out(m_size(2),m_size(1)*m_size(3)))
  END IF  

  val_out=&
  &RESHAPE(&
  &RESHAPE(val_in,[m_size(2)*m_size(3),m_size(1)],pad,[2,1])&
  &,[m_size(2),m_size(1)*m_size(3)]) 

 END SUBROUTINE array2mat_2B

!---------------------------------------------------------------------------

 SUBROUTINE array2mat_3B(n,val_in,val_out)
 !ARRAY2MAT_3B As array2mat_2 but now for 3D LOGICAL arrays.  

  INTEGER,INTENT(in)::n
  LOGICAL,INTENT(in)::val_in(:,:,:)
  LOGICAL,ALLOCATABLE,INTENT(inout)::val_out(:,:)
  INTEGER::m_size(3)
  LOGICAL::pad(2)

  pad=.FALSE.
  IF(n.GT.3) THEN
   WRITE(*,*) 'mod_constants.array2mat: n may not exceede dimension input&
   & array.'
   STOP
  END IF

  m_size(1)=MAXVAL([1,PRODUCT(SHAPE(val_in),[1,2,3].LT.n)])
  m_size(2)=MAXVAL([1,PRODUCT(SHAPE(val_in),[1,2,3].EQ.n)])
  m_size(3)=MAXVAL([1,PRODUCT(SHAPE(val_in),[1,2,3].GT.n)])

  IF(ALLOCATED(val_out)) THEN
   IF(ANY(SHAPE(val_out).NE.[m_size(2),m_size(1)*m_size(3)])) THEN
    DEALLOCATE(val_out)
   END IF
  END IF
  IF(.NOT.ALLOCATED(val_out)) THEN
    ALLOCATE(val_out(m_size(2),m_size(1)*m_size(3)))
  END IF  

  val_out=&
  &RESHAPE(&
  &RESHAPE(val_in,[m_size(2)*m_size(3),m_size(1)],pad,[2,1])&
  &,[m_size(2),m_size(1)*m_size(3)]) 

 END SUBROUTINE array2mat_3B

!----------------------------------------------------------------------------

 SUBROUTINE array2mat_4B(n,val_in,val_out)
 !ARRAY2MAT_4D As array2mat_2D but now for 4D LOGICAL arrays 

  INTEGER,INTENT(in)::n
  LOGICAL,INTENT(in)::val_in(:,:,:,:)
  LOGICAL,ALLOCATABLE,INTENT(inout)::val_out(:,:)
  INTEGER::m_size(3)
  LOGICAL::pad(2)

  pad=.FALSE.
  IF(n.GT.4) THEN
   WRITE(*,*) 'mod_constants.array2mat: n may not exceede dimension input&
   & array.'
   STOP
  END IF

  m_size(1)=MAXVAL([1,PRODUCT(SHAPE(val_in),[1,2,3,4].LT.n)])
  m_size(2)=MAXVAL([1,PRODUCT(SHAPE(val_in),[1,2,3,4].EQ.n)])
  m_size(3)=MAXVAL([1,PRODUCT(SHAPE(val_in),[1,2,3,4].GT.n)])

  IF(ALLOCATED(val_out)) THEN
   IF(ANY(SHAPE(val_out).NE.[m_size(2),m_size(1)*m_size(3)])) THEN
    DEALLOCATE(val_out)
   END IF
  END IF
  IF(.NOT.ALLOCATED(val_out)) THEN
    ALLOCATE(val_out(m_size(2),m_size(1)*m_size(3)))
  END IF  

  val_out=&
  &RESHAPE(&
  &RESHAPE(val_in,[m_size(2)*m_size(3),m_size(1)],pad,[2,1])&
  &,[m_size(2),m_size(1)*m_size(3)]) 

 END SUBROUTINE array2mat_4B

!----------------------------------------------------------------------------

 SUBROUTINE mat2array_2B(n,val_in,val_out)
 !MAT2ARRAY_2B as mat2array_2d but now for logicals

  INTEGER,INTENT(in)::n
  LOGICAL,ALLOCATABLE,INTENT(in)::val_in(:,:)
  LOGICAL,ALLOCATABLE,INTENT(inout)::val_out(:,:)
  INTEGER::m_size(3)
  LOGICAL::pad(2)

  pad=.FALSE.
  IF(n.GT.2) THEN
   WRITE(*,*) 'mod_constants.array2mat: n may not exceede dimension input&
   & array.'
   STOP
  END IF
  IF(SIZE(val_in).NE.SIZE(val_out)) THEN
   WRITE(*,*) 'mod_constants.mat2array: dimensions val_in and val_out are&
   & incongruent.'
   STOP
  END IF

  m_size(1)=MAXVAL([1,PRODUCT(SHAPE(val_out),[1,2].LT.n)])
  m_size(2)=MAXVAL([1,PRODUCT(SHAPE(val_out),[1,2].EQ.n)])
  m_size(3)=MAXVAL([1,PRODUCT(SHAPE(val_out),[1,2].GT.n)])

  val_out=RESHAPE(&
  &RESHAPE(val_in,[m_size(1),m_size(2)*m_size(3)],pad,[2,1])&
  &,SHAPE(val_out))  

 END SUBROUTINE mat2array_2B

!----------------------------------------------------------------------------

 SUBROUTINE mat2array_3B(n,val_in,val_out)
 !MAT2ARRAY_3D as mat2array_2D but now for 3D LOGICAL array

  INTEGER,INTENT(in)::n
  LOGICAL,ALLOCATABLE,INTENT(in)::val_in(:,:)
  LOGICAL,ALLOCATABLE,INTENT(inout)::val_out(:,:,:)
  INTEGER::m_size(3)
  LOGICAL::pad(2)

  pad=.FALSE.
  IF(n.GT.3) THEN
   WRITE(*,*) 'mod_constants.array2mat: n may not exceede dimension input&
   & array.'
   STOP
  END IF
  IF(SIZE(val_in).NE.SIZE(val_out)) THEN
   WRITE(*,*) 'mod_constants.mat2array: dimensions val_in and val_out are&
   & incongruent.'
   STOP
  END IF

  m_size(1)=MAXVAL([1,PRODUCT(SHAPE(val_out),[1,2,3].LT.n)])
  m_size(2)=MAXVAL([1,PRODUCT(SHAPE(val_out),[1,2,3].EQ.n)])
  m_size(3)=MAXVAL([1,PRODUCT(SHAPE(val_out),[1,2,3].GT.n)])

  val_out=RESHAPE(&
  &RESHAPE(val_in,[m_size(1),m_size(2)*m_size(3)],pad,[2,1])&
  &,SHAPE(val_out))  

 END SUBROUTINE mat2array_3B

!----------------------------------------------------------------------------

 SUBROUTINE mat2array_4B(n,val_in,val_out)
 !MAT2ARRAY_4B As mat2array_2D but now for LOGICAL array

  INTEGER,INTENT(in)::n
  LOGICAL,ALLOCATABLE,INTENT(in)::val_in(:,:)
  LOGICAL,ALLOCATABLE,INTENT(inout)::val_out(:,:,:,:)
  INTEGER::m_size(3)
  LOGICAL::pad(2)

  pad=.FALSE.
  IF(n.GT.4) THEN
   WRITE(*,*)'mod_constants.array2mat: n may not exceede dimension input&
   & array.'
   STOP
  END IF
  IF(SIZE(val_in).NE.SIZE(val_out)) THEN
   WRITE(*,*) 'mod_constants.mat2array: dimensions val_in and val_out are&
   & incongruent.'
   STOP
  END IF

  m_size(1)=MAXVAL([1,PRODUCT(SHAPE(val_out),[1,2,3,4].LT.n)])
  m_size(2)=MAXVAL([1,PRODUCT(SHAPE(val_out),[1,2,3,4].EQ.n)])
  m_size(3)=MAXVAL([1,PRODUCT(SHAPE(val_out),[1,2,3,4].GT.n)])

  val_out=RESHAPE(&
  &RESHAPE(val_in,[m_size(1),m_size(2)*m_size(3)],pad,[2,1])&
  &,SHAPE(val_out))  

 END SUBROUTINE mat2array_4B

!----------------------------------------------------------------------------
 
END MODULE
