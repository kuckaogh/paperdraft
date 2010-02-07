SUBROUTINE HEIGHT_SMOOTH(volsum,NXm,NYm)  

USE COMDAT
INTEGER(4)   ,INTENT(IN)                  :: NXm,NYm
REAL(8)    ,DIMENSION( 0:NXm+1, 0:NYm+1)        ,INTENT(INOUT)    :: volsum
REAL(8)    ,DIMENSION(-2:NXm+3,-2:NYm+3)              :: h_sm,h_sm_star 

h_sm          = 0
h_sm_star     = 0

! smoothing height
!-------------------------------------------------------------------------------------------
DO i=1,NXm
DO m=1,NYm
h_sm(i,m)   = volsum(i,m)
ENDDO
ENDDO

!h_sm(-1)  = volsumOUT(1)
!h_sm( 0)  = volsumOUT(1)
!h_sm(NX+1)= volsumOUT(NX)
!h_sm(NX+2)= volsumOUT(NX)

DO m=1,NYm
h_sm(-1,m)  = volsum(1,m)
h_sm(-0,m)  = volsum(1,m)
h_sm(NXm+1,m)= volsum(NXm,m)
h_sm(NXm+2,m)= volsum(NXm,m)
ENDDO

DO i=1,NXm
h_sm(i,-1)  = volsum(i,1)
h_sm(i,-0)  = volsum(i,1)
h_sm(i,NYm+1)= volsum(i,NYm)
h_sm(i,NYm+2)= volsum(i,NYm)
ENDDO

h_sm(0,0)      = 0.5* (h_sm(0,1)     + h_sm(1,0) )
h_sm(NXm+1,0)   = 0.5* (h_sm(NXm,   0) + h_sm(NXm+1, 1) )
h_sm(NXm+1,NYm+1)= 0.5* (h_sm(NXm,NYm+1) + h_sm(NXm+1,NYm) )
h_sm(0,NYm+1)   = 0.5* (h_sm( 1,NYm+1) + h_sm(   0,NYm) )

h_sm(-1,-1)     =   h_sm(0,0) 
h_sm(NXm+2,-1)   = h_sm(NXm+1,0) 
h_sm(NXm+2,NYm+2) = h_sm(NXm+1,NYm+1)
h_sm(-1,NYm+2)   = h_sm(0,NYm+1) 

!DO i=1,NX
!h_sm_star(i)= -0.0625*h_sm(i-2)+0.25*h_sm(i-1)+0.625*h_sm(i)+0.25*h_sm(i+1)-0.0625*h_sm(i+2)
!ENDDO

DO i=1,NXm
DO m=NYm,1,-1
!h_sm_star(i,m)= 0.375*h_sm(i,m)+0.125*(h_sm(i-1,m)+h_sm(i+1,m)+h_sm(i,m-1)+h_sm(i,m+1)) &
!               +1./64.*(h_sm(i-2,m)+h_sm(i+2,m)+h_sm(i,m-2)+h_sm(i,m+2))      &
!               +1./64.*(h_sm(i-1,m-1)+h_sm(i+1,m+1)+h_sm(i+1,m-1)+h_sm(i-1,m+1))
IF (NYm >= 3) THEN

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

ELSEIF (NYm>=2) THEN

h_sm_star(i,m)= 0.5*h_sm(i,m)+0.125*(h_sm(i-1,m)+h_sm(i+1,m)+h_sm(i,m-1)+h_sm(i,m+1))

ELSE

h_sm_star(i,m)= -0.0625*h_sm(i-2,m)+0.25*h_sm(i-1,m)+0.625*h_sm(i,m)+0.25*h_sm(i+1,m)-0.0625*h_sm(i+2,m)

ENDIF
!h_sm_star(i,m)= 0.2*h_sm(i,m)+0.2*(h_sm(i-1,m)+h_sm(i+1,m)+h_sm(i,m-1)+h_sm(i,m+1))
ENDDO
ENDDO

DO i=1,NXm
DO m=1,NYm
 volsum(i,m) = h_sm_star(i,m)
ENDDO
ENDDO

DO m=1,NYm
volsum(NXm+1,m)=volsum(NXm,m)
volsum(0,m)=volsum(1,m)
ENDDO

DO i=1,NXm
volsum(i,NYm+1)=volsum(i,NYm)
volsum(i,0)=volsum(i,1)
ENDDO

!  volsum = h0/dz
!-------------------------------------------------------------------------------------------

ENDSUBROUTINE HEIGHT_SMOOTH

