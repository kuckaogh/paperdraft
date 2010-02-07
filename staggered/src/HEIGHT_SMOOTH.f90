SUBROUTINE HEIGHT_SMOOTH(volsum)   

USE COMDAT

REAL(8)    ,DIMENSION( 0:NX+1, 0:NY+1)        ,INTENT(INOUT)    :: volsum
REAL(8)    ,DIMENSION(-2:NX+3,-2:NY+3)              :: h_sm,h_sm_star 

h_sm          = 0
h_sm_star     = 0

! smoothing height
!-------------------------------------------------------------------------------------------
DO i=1,NX
DO m=1,NY
h_sm(i,m)   = volsum(i,m)
ENDDO
ENDDO

!h_sm(-1)  = volsumOUT(1)
!h_sm( 0)  = volsumOUT(1)
!h_sm(NX+1)= volsumOUT(NX)
!h_sm(NX+2)= volsumOUT(NX)

DO m=1,NY
h_sm(-1,m)  = volsum(1,m)
h_sm(-0,m)  = volsum(1,m)
h_sm(NX+1,m)= volsum(NX,m)
h_sm(NX+2,m)= volsum(NX,m)
ENDDO

DO i=1,NX
h_sm(i,-1)  = volsum(i,1)
h_sm(i,-0)  = volsum(i,1)
h_sm(i,NY+1)= volsum(i,NY)
h_sm(i,NY+2)= volsum(i,NY)
ENDDO

h_sm(0,0)      = 0.5* (h_sm(0,1)     + h_sm(1,0) )
h_sm(NX+1,0)   = 0.5* (h_sm(NX,   0) + h_sm(NX+1, 1) )
h_sm(NX+1,NY+1)= 0.5* (h_sm(NX,NY+1) + h_sm(NX+1,NY) )
h_sm(0,NY+1)   = 0.5* (h_sm( 1,NY+1) + h_sm(   0,NY) )

h_sm(-1,-1)     =   h_sm(0,0) 
h_sm(NX+2,-1)   = h_sm(NX+1,0) 
h_sm(NX+2,NY+2) = h_sm(NX+1,NY+1)
h_sm(-1,NY+2)   = h_sm(0,NY+1) 

!DO i=1,NX
!h_sm_star(i)= -0.0625*h_sm(i-2)+0.25*h_sm(i-1)+0.625*h_sm(i)+0.25*h_sm(i+1)-0.0625*h_sm(i+2)
!ENDDO

DO i=1,NX
DO m=NY,1,-1
!h_sm_star(i,m)= 0.375*h_sm(i,m)+0.125*(h_sm(i-1,m)+h_sm(i+1,m)+h_sm(i,m-1)+h_sm(i,m+1)) &
!               +1./64.*(h_sm(i-2,m)+h_sm(i+2,m)+h_sm(i,m-2)+h_sm(i,m+2))      &
!               +1./64.*(h_sm(i-1,m-1)+h_sm(i+1,m+1)+h_sm(i+1,m-1)+h_sm(i-1,m+1))
IF (NY >= 3) THEN

h_sm_star(i,m)=    4./9.*h_sm(i,m)                                           &
                +0.75/9.*(h_sm(i-1,m)+h_sm(i+1,m)+h_sm(i,m-1)+h_sm(i,m+1))   &
                + 0.5/9.*(h_sm(i-1,m-1)+h_sm(i+1,m-1)+h_sm(i-1,m+1)+h_sm(i+1,m+1))  
!h_sm_star(i,m)= (9.5-4*sqrt(2.0))/16.*h_sm(i,m)    &
!               +0.125*(h_sm(i-1,m)+h_sm(i+1,m)+h_sm(i,m-1)+h_sm(i,m+1)) &
!               +1./(8.*sqrt(2.0))*(h_sm(i-1,m-1)+h_sm(i+1,m+1)+h_sm(i+1,m-1)+h_sm(i-1,m+1))  &
!               -1./64.*(h_sm(i-2,m)+h_sm(i+2,m)+h_sm(i,m-2)+h_sm(i,m+2))                     &
!               -1./(64.*2.)*(h_sm(i-2,m)+h_sm(i+2,m)+h_sm(i,m-2)+h_sm(i,m+2))
               

!h_sm_star(i,m)= 0.625*h_sm(i,m)+0.125*(h_sm(i-1,m)+h_sm(i+1,m)+h_sm(i,m-1)+h_sm(i,m+1)) &
!               -1./64.*(h_sm(i-2,m)+h_sm(i+2,m)+h_sm(i,m-2)+h_sm(i,m+2))      &
!               -1./64.*(h_sm(i-1,m-1)+h_sm(i+1,m+1)+h_sm(i+1,m-1)+h_sm(i-1,m+1))

ELSEIF (NY>=2) THEN

h_sm_star(i,m)= 0.5*h_sm(i,m)+0.125*(h_sm(i-1,m)+h_sm(i+1,m)+h_sm(i,m-1)+h_sm(i,m+1))

ELSE

h_sm_star(i,m)= -0.0625*h_sm(i-2,m)+0.25*h_sm(i-1,m)+0.625*h_sm(i,m)+0.25*h_sm(i+1,m)-0.0625*h_sm(i+2,m)

ENDIF
!h_sm_star(i,m)= 0.2*h_sm(i,m)+0.2*(h_sm(i-1,m)+h_sm(i+1,m)+h_sm(i,m-1)+h_sm(i,m+1))
ENDDO
ENDDO

DO i=1,NX
DO m=1,NY
 volsum(i,m) = h_sm_star(i,m)
ENDDO
ENDDO

DO m=1,NY
volsum(NX+1,m)=volsum(NX,m)
volsum(0,m)=volsum(1,m)
ENDDO

DO i=1,NX
volsum(i,NY+1)=volsum(i,NY)
volsum(i,0)=volsum(i,1)
ENDDO

!  volsum = h0/dz
!-------------------------------------------------------------------------------------------

ENDSUBROUTINE HEIGHT_SMOOTH

