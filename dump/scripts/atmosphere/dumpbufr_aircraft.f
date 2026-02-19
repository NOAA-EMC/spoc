!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      subroutine rddump(file,ret)

      character(255) file
      character(8)   subset
      integer        lunin/20/
      real(8)        tab
      integer(4),intent(out) :: ret  
     
      if(file/='readns') then
         open(lunin,file=file,form='unformatted') 
         call openbf(lunin,'IN',lunin)
         call ufbtab(-lunin,tab,1,1,iret,'count')
         ret = iret
      elseif(ireadns(lunin,subset,idate)==0) then
         call rdsubs(lunin,subset,idate)
         ret = 0
      else
         ret = -1
      endif
        
      end subroutine
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      subroutine rdsubs(lunin,subset,idate)

      implicit none

      REAL(8),PARAMETER :: RSTAR = 1.98        
      REAL(8),PARAMETER :: TZERO = 273.16
      REAL(8),PARAMETER :: EVLAT = 597.3
      REAL(8),PARAMETER :: VMOLW = 18.016
      REAL(8),PARAMETER :: DMOLW = 28.966
      REAL(8),PARAMETER :: EZERO = 6.11
      REAL(8),PARAMETER :: EPSLN = VMOLW/DMOLW 

      common/data/  stid,acfn,actn
      common/data/  otyp,styp,poaf,acns
      common/data/  year,mnth,days,hour,minu,seco
      common/data/  flat,flon,pres,elev,ialr
      common/data/  ovat,ovsh,ovew,ovnw
      common/data/  qmat,qmsh,qmew,qmnw
      common/data/  oeat,oesh,oeew,oenw

      character(8)  stid,acfn,actn
      real(8)       otyp,styp,poaf,acns
      real(8)       year,mnth,days,hour,minu,seco
      real(8)       flat,flon,pres,elev,ialr
      real(8)       ovat,ovsh,ovew,ovnw
      real(8)       qmat,qmsh,qmew,qmnw
      real(8)       oeat,oesh,oeew,oenw

      integer lunin,idate,iret,ioret,ibfms,ityp
      real(8) prlc,flvl,ialt,psal,heit,hmsl,flvlst
      real(8) mixr,rehu,tmdp,rahu,qmdd,mstq
      real(8) tmdb,tmdx,WDIR,WSPD,QMWN
      real(8) rpid,acid,acrn,dpof
      real(8) pi180/0.017453293/

      character(8) sid,subset
      equivalence (sid,rid)
      real(8) rid
      real(8) tiny/.00000001/,fill/10d10/

      real(8) HGTF_HI,HGTF_LO,AS,TFRMQP,ES,QFRMTP,P8,PR,PRS,Q,Z
      real(8) ESPH,RELH,ERLH,QSPH,EMIX,EDEW,WMIX,DPAL,RMIX,VIRT,SENT
      real(8) RH,DP,TS,TV,E,P,T,R,D,W

!-----------------------------------------------------------------------
!  FCNS HGTF_HI, HGTF_LO CALC. Z FROM P < 226.3MB AND P > 226.3MB; RESP
!-----------------------------------------------------------------------
      HGTF_HI(P8) = 11000 - LOG(P8/226.3_8)/1.576106E-4
      HGTF_LO(P8) = (1.-(P8/1013.25)**(1./5.256))*(288.15/.0065)
!-----------------------------------------------------------------------
!  Fcns below estimate pressure (mb) using indicated altitude (m) via
!-----------------------------------------------------------------------
      PR(Z) = 1013.25 * (((288.15 - (.0065 * Z))/288.15)**5.256)
      PRS(Z) = 226.3 * EXP(1.576106E-4 * (11000. - Z))
!-----------------------------------------------------------------------
!  Fcns below calculate various moisture variabels
!-----------------------------------------------------------------------
      ES(T)        = EZERO*EXP(17.269*T/(T+237.3))
      RH(R)        = .01*R
      ERLH(P,R,T)  = P/(1.+P/(RH(R)*ES(T))-1./RH(R))
      EMIX(P,W)    = W*P/(EPSLN+W)
      EDEW(DP)     = ES(DP)
      QSPH(P,E)    = EPSLN*E/(P-E*(1.-EPSLN))
      !RELH(P,E,T)  = (E/(P-E))/(ES(T)/(P-ES(T)))
      !WMIX(P,E)    = EPSLN*E/(P-E)
      !ESPH(P,Q)    = Q*P/(EPSLN+Q*(1.-EPSLN))
      !DPAL(E)      = 237.3* LOG(E/EZERO)/(17.269-LOG(E/EZERO))
      !RMIX(P,E)    = EPSLN*E/(P-E)
      !VIRT(P,E,TS) = (TS+TZERO)*(1.+RMIX/EPSLN)/(1.+RMIX)-TZERO
      !SENT(P,E,TV) = (TV+TZERO)*(1.+RMIX)/(1.+RMIX/EPSLN)-TZERO
!-----------------------------------------------------------------------

!  read basic elements from the aircraft dumpfile
!  ----------------------------------------------

      call ufbint(lunin,RPID,1,1,iret,'RPID')
      call ufbint(lunin,ACID,1,1,iret,'ACID')
      call ufbint(lunin,ACRN,1,1,iret,'ACRN')
      call ufbint(lunin,YEAR,1,1,iret,'YEAR')
      call ufbint(lunin,MNTH,1,1,iret,'MNTH')
      call ufbint(lunin,DAYS,1,1,iret,'DAYS')
      call ufbint(lunin,HOUR,1,1,iret,'HOUR')
      call ufbint(lunin,MINU,1,1,iret,'MINU')
      call ufbint(lunin,SECO,1,1,iret,'SECO')
      call ufbint(lunin,TMDB,1,1,iret,'TMDB')
      call ufbint(lunin,TMDX,1,1,iret,'TMDBST')
      call ufbint(lunin,WDIR,1,1,iret,'WDIR')
      call ufbint(lunin,WSPD,1,1,iret,'WSPD')
      call ufbint(lunin,QMAT,1,1,iret,'QMAT')
      call ufbint(lunin,QMWN,1,1,iret,'QMWN')
      call ufbint(lunin,MIXR,1,1,iret,'MIXR')
      call ufbint(lunin,REHU,1,1,iret,'REHU')
      call ufbint(lunin,RAHU,1,1,iret,'RAWHU')
      call ufbint(lunin,QMDD,1,1,iret,'QMDD')
      call ufbint(lunin,TMDP,1,1,iret,'TMDP')
      call ufbint(lunin,MSTQ,1,1,iret,'MSTQ')
      call ufbint(lunin,IALR,1,1,iret,'IALR')
      call ufbint(lunin,ACNS,1,1,iret,'ACNS')

      if(ibfms(acns)==1) acns = 7

!  aircraft call ids
!  -----------------

      stid=' ';acfn=' ';actn=' '
      if(ibfms(rpid)==0) write(stid,'(a8)') rpid
      if(ibfms(acid)==0) write(acfn,'(a8)') acid
      if(ibfms(acrn)==0) write(actn,'(a8)') acrn
      if(ibfms(rpid)==1) write(stid,'(a8)') 'MISSING '
      if(ibfms(acid)==1) write(acfn,'(a8)') 'MISSING '
      if(ibfms(acrn)==1) write(actn,'(a8)') 'MISSING '

!  CHECK THE DATE TIME ELEMENTS
!  ----------------------------

      if(ibfms(year)==1) year = 1 
      if(ibfms(mnth)==1) mnth = 1 
      if(ibfms(days)==1) days = 1 
      if(ibfms(hour)==1) hour = 0 
      if(ibfms(minu)==1) minu = 0 
      if(ibfms(seco)==1) seco = 0 

!  TEMPERATURE AND QUALITY MARKS
!  -----------------------------
      
      IF(IBFMS(TMDB)==1) TMDB=TMDX
      IF(IBFMS(TMDB)==0) then
         OVAT=TMDB
         IF(IBFMS(QMAT)==1) QMAT=2
      else
         OVAT=FILL
         QMAT=FILL 
      ENDIF

!  WIND AND QUALITY MARKS
!  ----------------------
      
      IF(IBFMS(MAX(WDIR,WSPD))==0) THEN
         IF(WSPD<=0.0)  THEN
            OVEW = 0.0
            OVNW = 0.0
         ELSE
            OVEW = -WSPD*SIN(WDIR*PI180)
            OVNW = -WSPD*COS(WDIR*PI180)
         END IF
         IF(IBFMS(QMWN)==1) THEN
            QMEW=2
            QMNW=2
         ELSE
            QMEW=QMWN
            QMNW=QMWN
         ENDIF
      ELSE
         OVEW=FILL
         OVNW=FILL
         QMEW=FILL
         QMNW=FILL
      ENDIF

!  ASSIGN A REPORT TYPE EACH DUMPFILE SUBSET CATEGORY                          |
!  --------------------------------------------------

      READ(SUBSET,'(5x,I3)',iostat=ioret) ITYP; styp=ityp

      IF(IORET==0) then  
         IF(STYP==001) OTYP=30    ! MTYP 004-001  Manual AIREP & ADS (AIREP)
         IF(STYP==002) OTYP=30    ! MTYP 004-002  Manual PIREP (PIREP)       
         IF(STYP==003) OTYP=31    ! MTYP 004-003  Automated AMDAR (FM-42 AMDAR)
         IF(STYP==004) OTYP=33    ! MTYP 004-004  Automated MDCRS (ARINC to NCEP) (BUFR)
         IF(STYP==005) OTYP=32    ! MTYP 004-005  Flight level reconnaissance (RECCO)
         IF(STYP==006) OTYP=31    ! MTYP 004-006  Automated European AMDAR (BUFR) 
         IF(STYP==007) OTYP=33    ! MTYP 004-007  Auto MDCRS (ARINC to AFWA to NCEP)
         IF(STYP==008) OTYP=34    ! MTYP 004-008  TAMDAR from MADIS (Mesaba) (NetCDF) 
         IF(STYP==009) OTYP=35    ! MTYP 004-009  Automated Canadian AMDAR (BUFR)
         IF(STYP==010) OTYP=34    ! MTYP 004-010  TAMDAR from Panasonic (BUFR)
         IF(STYP==011) OTYP=31    ! MTYP 004-011  Automated Korean AMDAR (BUFR) 
         IF(STYP==012) OTYP=34    ! MTYP 004-012  TAMDAR from MADIS (PenAir) (NetCDF)
         IF(STYP==013) OTYP=34    ! MTYP 004-013  TAMDAR from MADIS (Chautauqua)(NetCDF)
         IF(STYP==014) OTYP=31    ! MTYP 004-014  Automated French AMDAR (BUFR)
         IF(STYP==015) OTYP=32    ! MTYP 004-015  High density recconnaissance obs (HDOB)
         IF(STYP==103) OTYP=31    ! MTYP 004-103  All other automated AMDAR (BUFR)
      ELSE
         STYP=FILL; OTYP=FILL 
      ENDIF

!  IF LOW RES LAT/LON MISSING, REPORT LIKELY CONTAINS HI RES LAT/LON
!  -----------------------------------------------------------------
 
      CALL UFBINT(LUNIN,FLON,1,1,IRET,'CLON')
      CALL UFBINT(LUNIN,FLAT,1,1,IRET,'CLAT')

      IF(IBFMS(FLON)/=0) CALL UFBINT(LUNIN,FLON,1,1,IRET,'CLONH')
      IF(IBFMS(FLAT)/=0) CALL UFBINT(LUNIN,FLAT,1,1,IRET,'CLATH')

!  TRY TO FIND THE FLIGHT LEVEL HEIGHT
!  -----------------------------------

      call ufbint(LUNIN,PSAL,1,1,iret,'PSAL')
      call ufbint(LUNIN,FLVL,1,1,iret,'FLVL')
      call ufbint(LUNIN,IALT,1,1,iret,'IALT')
      call ufbint(LUNIN,PRLC,1,1,iret,'PRLC')
      call ufbint(LUNIN,HEIT,1,1,iret,'HEIT')
      call ufbint(LUNIN,HMSL,1,1,iret,'HMSL')
      call ufbint(LUNIN,FLVLST,1,1,iret,'FLVLST')

      IF(IBFMS(PRLC)==0)  THEN
         IF(PRLC.LT.22630) ELEV = HGTF_HI(PRLC*.01)
         IF(PRLC.GE.22630) ELEV = HGTF_LO(PRLC*.01)
      ELSEIF(IBFMS(IALT)==0)  THEN
         ELEV = IALT + SIGN(tiny,IALT)
      ELSEIF(IBFMS(PSAL)==0) THEN
         ELEV = PSAL + SIGN(tiny,PSAL)
      ELSEIF(IBFMS(FLVL)==0)  THEN
         ELEV = FLVL + SIGN(tiny,FLVL)
      ELSEIF(IBFMS(HEIT)==0)  THEN
         ELEV = HEIT + SIGN(tiny,HEIT)
      ELSEIF(IBFMS(HMSL)==0)  THEN
         ELEV = HMSL + SIGN(tiny,HMSL)
      ELSEIF(IBFMS(FLVLST)==0)  THEN
         ELEV = FLVLST + SIGN(tiny,FLVLST)
      ENDIF

!  CALCULATE PRESSURE IF PRLC IS MISSING
!  -------------------------------------

      IF(IBFMS(PRLC)==1) THEN
         IF(NINT(ELEV).LE.11000) PRES = PR(ELEV)
         IF(NINT(ELEV).GT.11000) PRES = PRS(ELEV)
      ELSE
         PRES = PRLC*0.1
      END If

!  CALCULATE Q SPECIFIC HUMIDITY
!  -----------------------------

      IF(IBFMS(MAX(MIXR,pres))==0) THEN
         P=PRES; W=MIXR
         OVSH=QSPH(P,EMIX(P,W))*1.e3
      ELSEIF(IBFMS(MAX(REHU,pres,TMDB))==0) THEN
         P=PRES; R=REHU; T=TMDB
         OVSH=QSPH(P,ERLH(P,R,T))*1.e3 
      ELSEIF(IBFMS(MAX(RAHU,pres,TMDB))==0) THEN
         P=PRES; R=RAHU; T=TMDB
         OVSH=QSPH(P,ERLH(P,R,T))*1.e3
      ELSEIF(IBFMS(MAX(TMDP,pres))==0) THEN
         P=PRES; D=TMDP
         OVSH=QSPH(P,EDEW(D))*1.e3
      ELSE
         OVSH=FILL  
      ENDIF

      IF(IBFMS(MSTQ)==0) QMSH=MSTQ 
      IF(IBFMS(QMDD)==0) QMSH=QMDD
      IF(OVSH==FILL)     QMSH=FILL
      IF(OVSH/=FILL)     QMSH=2   

C  GET PHASE OF FLIGHT from DPOF or POAF elements
C  ----------------------------------------------

         CALL UFBINT(LUNIN,DPOF,1,1,IRET,'DPOF')
         CALL UFBINT(LUNIN,POAF,1,1,IRET,'POAF')

         IF(IBFMS(DPOF)==0) THEN
           POAF = DPOF
           IF(INT(DPOF)>=7 .AND.INT(DPOF)<=10) POAF=5
           IF(INT(DPOF)>=11.AND.INT(DPOF)<=14) POAF=6
         ENDIF

      END SUBROUTINE
!----------------------------------------------------------------------
!----------------------------------------------------------------------
