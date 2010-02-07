SUBROUTINE OUTPUT_INIT_colocated(voll,volsum,uu,vv,ww,pp,ff,den,NXm,NYm,NZm)

USE COMDAT
INTEGER(4),INTENT(IN)           :: NXm,NYm,NZm
REAL(8)   ,DIMENSION(0:NXm+1,0:NYm+1) ,INTENT(IN) :: volsum
REAL(8)   ,DIMENSION(0:NXm+1,0:NYm+1,0:NZm+1)    :: den
REAL(8)   ,DIMENSION(0:NXm+1,0:NYm+1,0:NZm+1)    :: pp,ff
REAL(8)   ,DIMENSION(0:NXm+1,0:NYm+1,0:NZm+1)    :: x, y, z
REAL(8)   ,DIMENSION(0:NXm+1,0:NYm+1,0:NZm+1)    :: voll
REAL(8)   ,DIMENSION(0:NXm+1,0:NYm+1)            :: volsum_new
REAL(8)   ,DIMENSION(0:NXm+1,0:NYm+1,0:NZm+1)    :: uu,vv,ww

x=0
y=0
z=0

!Plot Average
!-----------------------------------------------------------------------
DO i=1,NXm
DO m=1,NYm
DO j=1,NZm
          x(i,m,j) = dx*(i-1)+dx/2                    ! Assigning x, y for plotting
          y(i,m,j) = dy*(m-1)+dy/2
          z(i,m,j) = dz*(j-1)+dz/2
ENDDO
ENDDO
ENDDO

z=z-h0
!volsum_new=volsum*dz-h0
volsum_new=volsum-h0

!NZ=20
WRITE(15,100)
WRITE(15,101) step*dt, NXm, NYm, NZm
WRITE(15,102) ((( x(i,m,j), i=1,nxm),m=1,NYm),j=1,NZm)
WRITE(15,102) ((( y(i,m,j), i=1,nxm),m=1,NYm),j=1,NZm)
WRITE(15,102) ((( z(i,m,j), i=1,nxm),m=1,NYm),j=1,NZm)
WRITE(15,102) (((den(i,m,j), i=1,nxm),m=1,NYm),j=1,NZm)
WRITE(15,102) ((( pp(i,m,j), i=1,nxm),m=1,NYm),j=1,NZm)
WRITE(15,102) ((( ff(i,m,j), i=1,nxm),m=1,NYm),j=1,NZm)
WRITE(15,102) ((( voll(i,m,j), i=1,nxm),m=1,NYm),j=1,NZm)

WRITE(15,102) ((( uu(i,m,j), i=1,nxm),m=1,NYm),j=1,NZm)
WRITE(15,102) ((( vv(i,m,j), i=1,nxm),m=1,NYm),j=1,NZm)
WRITE(15,102) ((( ww(i,m,j), i=1,nxm),m=1,NYm),j=1,NZm)

DO j=1,NZm ;DO m=1,NYm ;DO i=1,NXm
WRITE(15,102) 0.5*( (uu(i,m,j+1)-uu(i,m,j-1))/dz - (ww(i+1,m,j)-ww(i-1,m,j))/dx ) 
ENDDO;ENDDO;ENDDO
100   format('VARIABLES = x, y, z, den, p, f, lvol, u, v, w, curl')
101   format(' zone t="',f9.4,'", i=',i5, ', j=',i5,', k=',i5,', f=block')
102   format(E15.6) !1pe15.6 could cause bug
!NZ=25
ENDSUBROUTINE OUTPUT_INIT_colocated