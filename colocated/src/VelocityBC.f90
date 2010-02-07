SUBROUTINE VelocityBC(uuB,vvB,wwB,NXm,NYm,NZm) !INOUT  
!_______________________________________________________________________________________

USE COMDAT
INTEGER(4)   ,INTENT(IN)                       :: NXm,NYm,NZm
REAL(8)    ,DIMENSION(0:NXm  ,0:NYm+1,0:NZm+1) ,INTENT(INOUT)   :: uuB
REAL(8)    ,DIMENSION(0:NXm+1,0:NYm  ,0:NZm+1) ,INTENT(INOUT)   :: vvB
REAL(8)    ,DIMENSION(0:NXm+1,0:NYm+1,0:NZm  ) ,INTENT(INOUT)   :: wwB

!________________________________________________________________

!VBC  1:side 2:top 3:left 4:bottom 5:right    >0 slip =<0 nonslip

IF(VBC(1)>0) THEN ! Side Slip 
        !---------------------------------- uuB --------
        DO i=1,NXm-1 ;DO j=0,NZm+1
        uuB(i,NYm+1,j)  =     uuB(i,NYm,j) 
        uuB(i,0,j)      =     uuB(i,1,j)
        ENDDO;ENDDO
        !---------------------------------- wwB --------
        DO j=1,NZm-1 ;DO i=0,NXm+1
        wwB(i,NYm+1,j)  =     wwB(i,NYm,j)
        wwB(i,0,j)      =     wwB(i,1,j)
        ENDDO;ENDDO
ELSE  !side nonslip
        !---------------------------------- uuB --------
        DO i=1,NXm-1 ;DO j=0,NZm+1
        uuB(i,NYm+1,j)  = 0 - uuB(i,NYm,j) 
        uuB(i,0,j)      = 0 - uuB(i,1,j)
        ENDDO;ENDDO
        !---------------------------------- wwB --------
        DO j=1,NZm-1 ;DO i=0,NXm+1
        wwB(i,NYm+1,j)  = 0 - wwB(i,NYm,j)
        wwB(i,0,j)      = 0 - wwB(i,1,j)
        ENDDO;ENDDO
ENDIF

IF(VBC(2)>0) THEN ! TopSlip 
        !---------------------------------- uuB --------
        DO i=1,NXm-1 ;DO m=0,NYm+1
        uuB(i,m,NZm+1)  =     uuB(i,m,NZm) 
        ENDDO;ENDDO
        !---------------------------------- vvB --------
        DO m=1,NYm-1 ;DO i=1-1,NXm+1 
        vvB(i,m,NZm+1)  =     vvB(i,m,NZm) 
        ENDDO;ENDDO
ELSE !Top nonlip
        !---------------------------------- uuB --------
        DO i=1,NXm-1 ;DO m=0,NYm+1
        uuB(i,m,NZm+1)  = 0 - uuB(i,m,NZm) 
        ENDDO;ENDDO
        !---------------------------------- vvB --------
        DO m=1,NYm-1 ;DO i=1-1,NXm+1 
        vvB(i,m,NZm+1)  = 0 - vvB(i,m,NZm) 
        ENDDO;ENDDO
ENDIF        

IF(VBC(3)>0) THEN ! Left Slip 
        !---------------------------------- vvB --------
        DO m=1,NYm-1 ;DO j=0,NZm+1
        vvB(0,m,j)      =     vvB(1,m,j)
        ENDDO;ENDDO
        !---------------------------------- wwB --------
        DO j=1,NZm-1 ;DO m=0,NYm+1
        wwB(0,m,j)      =     wwB(1,m,j)
        ENDDO;ENDDO
ELSE ! Left Nonslip
        !---------------------------------- vvB --------
        DO m=1,NYm-1 ;DO j=0,NZm+1
        vvB(0,m,j)      = 0 - vvB(1,m,j)
        ENDDO;ENDDO
        !---------------------------------- wwB --------
        DO j=1,NZm-1 ;DO m=0,NYm+1
        wwB(0,m,j)      = 0 - wwB(1,m,j)
        ENDDO;ENDDO
ENDIF        

IF(VBC(4)>0) THEN ! Bottom Slip 
        !---------------------------------- uuB --------
        DO i=1,NXm-1 ;DO m=0,NYm+1 
        uuB(i,m,0)      =     uuB(i,m,1)
        ENDDO;ENDDO
        !---------------------------------- vvB --------
        DO m=1,NYm-1 ;DO i=1-1,NXm+1  
        vvB(i,m,0)      =     vvB(i,m,1)
        ENDDO;ENDDO
ELSE !Bottom nonlip
        !---------------------------------- uuB --------
        DO i=1,NXm-1 ;DO m=0,NYm+1 
        uuB(i,m,0)      = 0 - uuB(i,m,1)
        ENDDO;ENDDO
        !---------------------------------- vvB --------
        DO m=1,NYm-1 ;DO i=1-1,NXm+1  
        vvB(i,m,0)      = 0 - vvB(i,m,1)
        ENDDO;ENDDO
ENDIF

IF (VBC(5)>0) THEN ! Right Slip
        !---------------------------------- vvB --------
        DO m=1,NYm-1 ;DO j=0,NZm+1
        vvB(NXm+1,m,j)  =     vvB(NXm,m,j) 
        ENDDO;ENDDO
        !---------------------------------- wwB --------
        DO j=1,NZm-1 ;DO m=0,NYm+1
        wwB(NXm+1,m,j)  =     wwB(NXm,m,j)
        ENDDO;ENDDO
ELSE ! Right Nonslip
        !---------------------------------- vvB --------
        DO m=1,NYm-1 ;DO j=0,NZm+1
        vvB(NXm+1,m,j)  = 0 - vvB(NXm,m,j) 
        ENDDO;ENDDO
        !---------------------------------- wwB --------
        DO j=1,NZm-1 ;DO m=0,NYm+1
        wwB(NXm+1,m,j)  = 0 - wwB(NXm,m,j)
        ENDDO;ENDDO
ENDIF

ENDSUBROUTINE VelocityBC 