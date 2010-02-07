SUBROUTINE INIT_KP_All(kp,NXk,NYk,NZk)

USE COMDAT

INTEGER(4)                               ,INTENT(IN)  :: NXk,NYk,NZk
REAL(8)    ,DIMENSION(1:NXk,1:NYk,1:NZk) ,INTENT(OUT) :: kp
!REAL(8) ::dx2,dz2,dy2

dx2 = (dx*dx)
dz2 = (dz*dz)
dy2 = (dy*dy)

kp      = 0

IF(NYk==1) THEN
 IF(OPENZ) THEN
 !***** kp *********************************************************	

 kp    = 2./dx2 + 2./dz2   ! for interior points

 DO i=1, NXk
  kp(i,1,1) = 2./dx2 + 1./dz2 !+2./dy2 !for boundary area
 ! kp(i,1,NZ_high)= 2./dx2 + 1./dz2 !+2./dy2 !for boundary area
 ENDDO

 DO j=1, NZk
  kp(1,1,j) = 1./dx2 + 2./dz2 !+2./dy2 !for boundary area
  kp(NXk,1,j)= 1./dx2 + 2./dz2 !+2./dy2 !for boundary area
 ENDDO

  kp(1 , 1, 1)= 1./dx2 + 1./dz2  !+1./dy2     !for corner points
  kp(NXk, 1, 1)= 1./dx2 + 1./dz2  !+1./dy2
 ! kp(NX, 1,NZ_high)= 1./dx2 + 1./dz2 !+1./dy2
 ! kp(1 , 1,NZ_high)= 1./dx2 + 1./dz2 !+1./dy2

 !*********************************************************************************
 ELSE
  WRITE(*,*) 'OPENZ in INIT_KP_All Error!';STOP
 ENDIF

ELSE
 WRITE(*,*) '2D in INIT_KP_All Error!';STOP
ENDIF

ENDSUBROUTINE INIT_KP_All