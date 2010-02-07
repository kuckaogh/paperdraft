SUBROUTINE STABILITY(uu,vv,ww)

USE COMDAT
REAL(8)    ,DIMENSION(0:NX  ,0:NY+1,0:NZ+1) ,INTENT(IN)    :: uu
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1,0:NZ  ) ,INTENT(IN)    :: ww
REAL(8)    ,DIMENSION(0:NX+1,0:NY  ,0:NZ+1) ,INTENT(IN)    :: vv


!INTEGER ,DIMENSION(3)  :: uu_max_loc, vv_max_loc, ww_max_loc
!REAL(8) :: uu_Max, vv_Max, ww_Max, tempcoeff

!MaxVelocity = Max( ABS(MAXVAL(u)),ABS(MAXVAL(v)),ABS(MINVAL(u)),ABS(MINVAL(v)) )
!MaxVelocity = MAX( MAXVAL(ABS(uu)),MAXVAL(ABS(MAXVAL(v)),ABS(MINVAL(u)),ABS(MINVAL(v)) )

!uu_max_loc = MAXLOC((uu))
!vv_max_loc = MAXLOC((vv))
!ww_max_loc = MAXLOC((ww))
!
!uu_Max = uu(uu_max_loc(1),uu_max_loc(2),uu_max_loc(3))
!vv_Max = vv(vv_max_loc(1),vv_max_loc(2),vv_max_loc(3))
!ww_Max = ww(ww_max_loc(1),ww_max_loc(2),ww_max_loc(3))

uu_Max_CL = MAXVAL(ABS(uu))*dt_o/dx
vv_Max_CL = MAXVAL(ABS(vv))*dt_o/dy
ww_Max_CL = MAXVAL(ABS(ww))*dt_o/dz
!uu2vv2ww2_Max = MAXVAL(uu*uu+vv*vv+ww*ww) 

!Do i= 1, nx
! Do j= 1, ny

!   UV(i,j) = u(i,j)*u(i,j) + v(i,j)*v(i,j) 

! Enddo
!Enddo

!print *, uu_Max,vv_Max,ww_Max

!IF(RBM) THEN
!    WRITE(18,*) step*subiteration,INT(REAL(step*subiteration)/REAL(NT) * 100),'%'
!ELSE
!	WRITE(18,*) step,INT(REAL(Step)/REAL(NT) * 100),'%' 
!ENDIF

Write(18,*) '----------------------------------------------'
Write(18,'(" Max UU Corant# =         ",E11.4)') uu_Max_CL
Write(18,'(" Max VV Corant# =         ",E11.4)') vv_Max_CL
Write(18,'(" Max WW Corant# =         ",E11.4)') ww_Max_CL
Write(18,'(" Max (6D*dt/Min(dx,dy,dz)^2) =       ",E11.4)') 6*Visc*dt_o/MIN(dx,dy,dz)**2
Write(18,*) '----------------------------------------------'  
Write(*,*) '----------------------------------------------'
Write(*,'(" Max UU Corant# =         ",E11.4)') uu_Max_CL
Write(*,'(" Max VV Corant# =         ",E11.4)') vv_Max_CL
Write(*,'(" Max WW Corant# =         ",E11.4)') ww_Max_CL
Write(*,'(" Max (6D*dt/Min(dx,dy,dz)^2) =       ",E11.4)') 6*Visc*dt_o/MIN(dx,dy,dz)**2
Write(*,*) '----------------------------------------------' 

!!Write(*,'(" Max 4D*dt/dx^2 =      ",E11.4)') 4*Visc*dt/(dx*dx)
!!Write(*,'(" Max (U^2+V^2)*dt/4D = ",E11.4)') MAXVAL(UV)*dt/(Visc*4.)

!tempcoeff = MAX(uu_Max*dt_o/dx,vv_Max*dt_o/dy,ww_Max*dt_o/dz,6*Visc*dt_o/MIN(dx,dy,dz)**2)
!
!IF (tempcoeff > 0.95) THEN
!dt_o = dt_o*0.5
!ELSE
!dt_o = ((0.75-tempcoeff) + 1 )*dt_o
!ENDIF


ENDSUBROUTINE