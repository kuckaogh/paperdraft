MODULE COMDAT

IMPLICIT NONE
REAL(8)     ,PARAMETER                      :: PI=3.141592653589793238462643383279502884197169399375

LOGICAL                                     :: pp_reset
INTEGER(4)                                  :: buffer, SubIteration,Receding, CASENZ ,SideSlip
INTEGER(4)                                  :: NX, NY, NZ, NXS, NYS, NZS, NZS2  !NX=60,NY=6,NZ=40
INTEGER(4)                                  :: NXS_start, NXS_end



LOGICAL                                     :: restart
REAL(4)                                     :: time_begin,time_end


REAL(8)                                     :: err,errsum
!------------- From input.dat -----------------------------------------------(
INTEGER(4)                                  :: Bossi,i_dfdxN,i_dfdx_vN
INTEGER(4)                                  :: ipltn,NT,height_smoothing
INTEGER(8)                                  :: ii_SOR,ii_SOR_Pos, iSOR, iSORPOS, itersum, itersum1, itersum2
INTEGER(8)                                  :: ii_SOR_Min,ISOR_ErrMax
INTEGER(4)                                  :: InterX
INTEGER(4)                                  :: CIP2D
INTEGER(4)                                  :: ANH 
!( 0: Hydrostatic 1: Nonhydrostatic without H Correction 2: Adjustable-Nonhydrostatic)
REAL(8)                                     :: beta, errmax, dt_o,dt,visc, uuwall
REAL(8)                                     :: dx,dy,dz,dz_high,dz_low
REAL(8)                                     :: den_w
REAL(8)                                     :: gravity,Am0,h0
REAL(8)                                     :: WBaC
REAL(8)                                     :: KdifX,KdifZ,KdifY
REAL(8)                                     :: ff_initial,ff_strati
INTEGER(4)                                  :: conservation ! One-Domain only, ff conservation
REAL(8)                                     :: conservation_power
LOGICAL                                     :: OPENZ
INTEGER(4)                                  :: SurfaceCorrection !(0:no correction Otherwise: common correct )
!------------- From input.dat -----------------------------------------------)
!REAL(8)    ,DIMENSION(0:NX+1,0:NY+1,0:NZ+1) :: kp
REAL(8)                                     :: dx2,dz2,dy2
INTEGER(4)                                  :: step,step_1,step_2,i,m,j
 CHARACTER                                   :: dummy
REAL(8)                                     :: time_sec
REAL(8)                                     :: transp_x, transp_y, transp_z
LOGICAL                                     :: Height_field
INTEGER(4)                                  :: irestartn
INTEGER(4)                                  :: TIME_ORDER
INTEGER(4)                                  :: pp_solve_1,pp_solve_2,pp_solve_3,pp_solve_4
!------------- For Output ------------
REAL(8)                                     :: uu_Max_CL,   vv_Max_CL,   ww_Max_CL
REAL(8)                                     :: uu_Max_CL_T, vv_Max_CL_T, ww_Max_CL_T
!------------------------------------simple----
INTEGER(4)             :: NXm,NYm,NZm
INTEGER(4)             :: NXm2,NYm2,NZm2
INTEGER(4)             :: NXm3,NYm3,NZm3
INTEGER(4)             :: NXm4,NYm4,NZm4
INTEGER(4)             :: NXm5,NYm5,NZm5
INTEGER(4)             :: MultiLevel,i_Multi1,i_Multi2,i_Multi3,i_Multi4,i_Multi5
REAL(8)                :: errLevel1,errLevel2,errLevel3,errLevel4
ENDMODULE COMDAT