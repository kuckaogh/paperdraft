SUBROUTINE OUTPUT_INIT(voll,volsum,uu,vv,ww,pp,ff,den)
!SUBROUTINE OUTPUT_INIT(ff,voll,volsum,den,uu,ww,p)

USE COMDAT

REAL(8)      ,DIMENSION(0:NX+1,0:NY+1) ,INTENT(IN) :: volsum
REAL(8)      ,DIMENSION(0:NX+1,0:NY+1,0:NZ+1)      :: den
REAL(8)      ,DIMENSION(0:NX+1,0:NY+1,0:NZ+1)      :: pp,ff
REAL(8)      ,DIMENSION(0:NX+1,0:NY+1,0:NZ+1)      :: x, y, z
REAL(8)      ,DIMENSION(0:NX+1,0:NY+1,0:NZ+1)      :: voll
REAL(8)      ,DIMENSION(0:NX+1,0:NY+1)             :: volsum_new
REAL(8)       uu(0:NX,0:NY+1,0:NZ+1), ww(0:NX+1,0:NY+1,0:NZ), vv(0:NX+1,0:NY,0:NZ+1)

x=0
y=0
z=0

!Plot Average
!-----------------------------------------------------------------------
DO i=1,NX
DO m=1,NY
DO j=1,NZ
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
WRITE(15,101) step*dt, NX, NY, NZ
WRITE(15,102) ((( x(i,m,j), i=1,nx),m=1,NY),j=1,NZ)
WRITE(15,102) ((( y(i,m,j), i=1,nx),m=1,NY),j=1,NZ)
WRITE(15,102) ((( z(i,m,j), i=1,nx),m=1,NY),j=1,NZ)
WRITE(15,102) (((den(i,m,j), i=1,nx),m=1,NY),j=1,NZ)
WRITE(15,102) ((( pp(i,m,j), i=1,nx),m=1,NY),j=1,NZ)
WRITE(15,102) ((( ff(i,m,j), i=1,nx),m=1,NY),j=1,NZ)
WRITE(15,102) ((( voll(i,m,j), i=1,nx),m=1,NY),j=1,NZ)

WRITE(15,102) ((( (uu(i-1,m,j)+uu(i,m,j))/2, i=1,nx),m=1,NY),j=1,NZ)
WRITE(15,102) ((( (vv(i,m-1,j)+vv(i,m,j))/2, i=1,nx),m=1,NY),j=1,NZ)
WRITE(15,102) ((( (ww(i,m,j-1)+ww(i,m,j))/2, i=1,nx),m=1,NY),j=1,NZ)
!WRITE(15,102) ((( volsum_new(i,m), i=1,nx),m=1,NY),j=1,NZ)
CALL VelocityBC(uu,vv,ww)
WRITE(15,102) ((( 0.25*((uu(i-1,m,j+1)+uu(i,m,j+1)-uu(i-1,m,j-1)-uu(i,m,j-1))/dz-(ww(i+1,m,j)+ww(i+1,m,j-1)-ww(i-1,m,j)-ww(i-1,m,j-1))/dx ), i=1,NX),m=1,NY),j=1,NZ)

100   format('VARIABLES = x, y, z, den, p, f, lvol, u, v, w, curl')
101   format(' zone t="',f9.4,'", i=',i5, ', j=',i5,', k=',i5,', f=block')
102   format(E15.6) !1pe15.6 could cause bug
!NZ=25
ENDSUBROUTINE OUTPUT_INIT