SUBROUTINE BC_wall_colocated(uu,vv,ww,NXm,NYm,NZm,wall)          
!___________________________________________

USE COMDAT

INTEGER(4)  ,INTENT(IN)           :: NXm,NYm,NZm,wall   ! 1:left 2:right else:both
REAL(8),DIMENSION(0:NXm+1,0:NYm+1,0:NZm+1),INTENT(INOUT) :: uu,vv,ww
!_____________________________________________________________

!Assign zero for no-flux wall
!-------------------------------

SELECTCASE (wall)

CASE(1) !left

 uu(0  , :   ,:) = -uu(1  , : ,:) !right side is not assigned
 vv(:  , 0   ,:) = -vv(:  , 1 ,:)
 vv(:  ,NYm+1,:) = -vv(:  ,NYm,:)
 ww(:  , :   ,0) = -ww(:  , : ,1) !top is not assigned

CASE(2) !right

 uu(NXm+1, :   ,:) = -uu(NXm, : ,:) !left side is not assigned
 vv( :   , 0   ,:) = -vv( : , 1 ,:)
 vv( :   ,NYm+1,:) = -vv( : ,NYm,:) 
 ww( :   , :   ,0) = -ww( : , : ,1) !top is not assigned
 
CASE DEFAULT !both

 uu(0    , :   ,:) = -uu(1  , : ,:)
 uu(NXm+1, :   ,:) = -uu(NXm, : ,:)
 vv( :   , 0   ,:) = -vv( : , 1 ,:)
 vv( :   ,NYm+1,:) = -vv( : ,NYm,:) 
 ww( :   , :   ,0) = -ww( : , : ,1) !top is not assigned

ENDSELECT

ENDSUBROUTINE BC_wall_colocated