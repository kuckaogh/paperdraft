SUBROUTINE BC_wall(uu,vv,ww,NXm,NYm,NZm,wall)          
!___________________________________________

USE COMDAT

INTEGER(4)  ,INTENT(IN)           :: NXm,NYm,NZm,wall   ! 1:left 2:right else:both
REAL(8)   ,DIMENSION(0:NXm  ,0:NYm+1,0:NZm+1) ,INTENT(INOUT) :: uu
REAL(8)   ,DIMENSION(0:NXm+1,0:NYm  ,0:NZm+1) ,INTENT(INOUT) :: vv
REAL(8)   ,DIMENSION(0:NXm+1,0:NYm+1,0:NZm  ) ,INTENT(INOUT) :: ww

!_____________________________________________________________


!Assign zero for no-flux wall
!-------------------------------

SELECTCASE (wall)

CASE(1) !left

 uu(0  , : ,:) = 0
 vv(:  , 0 ,:) = 0
 vv(:  ,NYm,:) = 0
 ww(:  , : ,0) = 0

CASE(2) !right

 uu(NXm, : ,:) = 0
 vv( : , 0 ,:) = 0
 vv( : ,NYm,:) = 0 
 ww( : , : ,0) = 0
 
CASE DEFAULT !all

 uu(0  , : ,:) = 0
 uu(NXm, : ,:) = 0
 vv( : , 0 ,:) = 0
 vv( : ,NYm,:) = 0 
 ww( : , : ,0) = 0

ENDSELECT

ENDSUBROUTINE BC_wall