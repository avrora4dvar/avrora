MODULE mod_time

CONTAINS

 FUNCTION unix2cal(utime) RESULT(idate)
  !utime  input  Unix system time, seconds since 1970.0
  !idate  output Array: 1=year, 2=month, 3=date, 4=hour, 5=minute, 6=secs
  !Author  Clive Page, Leicester University, UK.   1995-MAY-2
  implicit none
  integer,intent(in)::utime
  integer::idate(6)
  integer::mjday, nsecs,tutime
  real::day

  !Note the MJD algorithm only works from years 1901 to 2099.
  tutime = utime
  !use tutime = utime - 3600 for daylight savings time!
  mjday    = int(tutime/86400 + 40587)
  idate(1) = 1858 + int( (mjday + 321.51) / 365.25)
  day      = aint( mod(mjday + 262.25, 365.25) ) + 0.5
  idate(2) = 1 + int(mod(day / 30.6 + 2.0, 12.0) )
  idate(3) = 1 + int(mod(day,30.6))
  nsecs    = mod(tutime, 86400)
  idate(6) = mod(nsecs, 60)
  nsecs    = nsecs / 60
  idate(5) = mod(nsecs, 60)
  idate(4) = nsecs / 60
 END FUNCTION unix2cal

!----------------------------------------------------------------------------

 FUNCTION cal2unix(idate) RESULT(utime)
  IMPLICIT NONE
  INTEGER,intent(in)::idate(6)
  INTEGER::utime
  INTEGER::gdate(6)

  !First guess
  utime=(idate(1)-1970)*31557600+(idate(2)-1)*2643840+&
  &(idate(3)-1)*86400+idate(4)*3600+idate(5)*60+idate(6)
  gdate=unix2cal(utime)

  !Get correct year and month
  DO WHILE(gdate(1).NE.idate(1))
   IF(gdate(1).LT.idate(1)) utime=utime+28*86400
   IF(gdate(1).GT.idate(1)) utime=utime-28*86400
   gdate=unix2cal(utime)
  END DO
  !Get correct month
  DO WHILE(gdate(2).NE.idate(2))
   IF(gdate(2).LT.idate(2)) utime=utime+28*86400
   IF(gdate(2).GT.idate(2)) utime=utime-28*86400
   gdate=unix2cal(utime)
  END DO    
  !Get correct day,hour,minute,second
  utime=utime+(idate(3)-gdate(3))*86400+&
  &(idate(4)-gdate(4))*3600+(idate(5)-gdate(5))*60+(idate(6)-gdate(6))

 END FUNCTION cal2unix

!-------------------------------------------------------------------------
 FUNCTION ref2unix(ref,iref) RESULT(utime)
 !Calculate the UNIX time for the time in seconds given with respect to the 
 !reference time ref

  INTEGER,intent(in)::iref(6),ref
  INTEGER::utime
 
  utime=ref+cal2unix(iref)

 END FUNCTION ref2unix


END MODULE
