      program MERGER

! version 091
! BJW, Riverside 2006, Praha 2007, 2008, 2009, Brno 2015
! simulates a minor merger (giant elliptical + dwarf elliptical)
! with Plummer and Hernquist potential, simulates gas with 
! sticky-particles scheme


      IMPLICIT REAL*8 (A-H,O-Z)     
      INCLUDE 'merge.inp'

      CHARACTER ProfileP*1, ProfileS*1, PhiP*1, Display*3
      REAL*4 TIMEI, TIMEF 
      COMMON /UNITS/ TUNIT,VUNIT,YVUNIT,GMUNIT,GMUNIT6
      COMMON /GXIES/  Xp(3),Xs(3),Vp(3),Vs(3)
      COMMON /GXIES2/ GMp, Bp, Bp2, YBp2, CTrBp, Sigma0p, Vesc0p, Sr0p,
     &                GMs, Bs, Bs2, YBs2, CTrBs, Sigma0s, Vesc0s, Sr0s,
     &                GMeff, GMStoP, GMpart1,GMpart2, Eps2,Eps2s,Eps2p, 
     &                Omega23, Dini, Vini, Trmin, DTrBp, DTrBs, Eps,
     &                RtrBp, RtrBs, BgasP, BgasS, BgasP2, BgasS2,
     &                YBgasP2, YBgasS2, CTrBgasP, CtrBgasS, Vesc0gasP,
     &                Vesc0gasS, Sr0gasP, Sr0gasS, DTrBgasP, DTrBgasS,
     &                RTrBgasP, RTrBgasS, Eps2gasP, Eps2gasS
      COMMON /TestP/  XT (NX), YT(NX), ZT(NX), VXT(NX), VYT(NX), VZT(NX)
      COMMON /TestP1/ N1, N2, NALL
      COMMON /Gas/    Xg(NX), Yg(NX), Zg(NX), VXg(NX), VYg(NX), VZg(NX),
     &                Xhalf(NX), Yhalf(NX), Zhalf(NX), Mesh(NX)
      COMMON /Gas1/   N1g, N2g, NALLg
      
      CALL SetUnits  ! set code units
c      print*, 'Units set'
      CALL INI1    (ProfileP, ProfileS, PhiP)! input par. and init. conditions
c      print*, 'INI1 done'
      CALL InitTPs (ProfileP, ProfileS)      ! initialize test particles 
c      print*, 'IniTPs done'
!      CALL InitPGmap ! initialize PGPLOT

      CALL CPU_TIME(TIMEI)

!      print*,'Display = xwin / gif / ps ?'
!      read*, Display(1:1)
      Display(1:1) = 'x'
      IF (Display(1:1) .EQ. 'g') THEN
        Display = 'gif'
      ELSEIF (Display(1:1) .EQ. 'p') THEN
        Display = 'cps'
      ELSE
        Display = 'win'
      ENDIF 

      ITend = 200000

      IF (PhiP .EQ. 'P') THEN     ! primary's potential = Plummer

        CALL EqMovG0P   ! initialize leapfrog for the primary+secondary 
        !print*, 'eqmovG0P'
        CALL EqMovT0P   ! initialize leapfrog for the test particles
        !print*, 'eqmovT0P'
        CALL EqMovTgas0P
        !print*, 'eqmovTgas0P'
c        CALL StickyParticles
        !print*, 'StickyParticles'

        DO 100 IT = 0, ITend

c          print*, IT*TUNIT
 
          CALL Disrupt (Drel, Vrel)    ! test for catastrophic discruption 
          !IF (MOD(IT,40000) .EQ. 0 .AND. IT .NE. 0) THEN
          IF (MOD(IT,10000) .EQ. 0) THEN
            print*, IT*TUNIT
            CALL EqMovG1P
            CALL EqMovT1Pprint
            CALL EqMovTgas1Pprint
            CALL StickyParticles
          ELSE
            CALL EqMovG1P   ! integrate eq. of motion for galaxies
            CALL EqMovT1P   ! integrate eq. of motion for the test particles
            CALL EqMovTgas1P
            CALL StickyParticles
          ENDIF

!          IF (MOD(IT,10) .EQ. 0) THEN
!            CALL CompSurfAP               ! compute surf. density analytically 
!            CALL CompSurfN (0.d0, 0.d0)   ! compute surf. density numerically
                                        !   (test particles)
!            CALL MAPSURF (Display, IT*TUNIT, Drel, Vrel*VUNIT)
!          ENDIF

 100    CONTINUE


      ELSEIF (PhiP .EQ. 'V') THEN    ! primary's potential = de Vaucouleurs 

        CALL Young76c0               ! initialize de Vaucouleurs table 

        CALL EqMovG0V   ! initialize leapfrog for the primary+secondary 
        CALL EqMovT0V   ! initialize leapfrog for the test particles
        
        DO 200 IT = 0, ITend

          print*, IT*TUNIT

          CALL Disrupt  (Drel, Vrel)    ! test for catastrophic discruption 
          CALL EqMovG1V (Drel)  !  eq. of motion for galaxies
          CALL EqMovT1V   ! integrate eq. of motion for the test particles 

!          IF (MOD(IT,10) .EQ. 0) THEN
!            CALL CompSurfAV               ! compute surf. density analytically 
!            CALL CompSurfN (0.d0, 0.d0)   ! compute surf. density numerically
!                                        !   (test particles)
!            CALL MAPSURF (Display, IT*TUNIT, Drel, Vrel*VUNIT)
!          ENDIF

 200   CONTINUE
 
      ELSE                        !primary potential = Hernquist
        
        CALL EqMovG0H           !initialize leapfrog for the primary+secondary
        CALL EqMovT0H           !initialize leapfrog for the test particles
        CALL EqMovTgas0H
        
        
        DO 300 IT = 0, ITend
        
c          print*, IT*TUNIT
          
          CALL Disrupt(Drel, Vrel)
          IF (MOD(IT,10000) .EQ. 0) THEN
            print*, IT*TUNIT
            CALL EqMovG1H
            CALL EqMovT1Hprint
            CALL EqMovTgas1Hprint
            CALL StickyParticles
          ELSE
            CALL EqMovG1H   ! integrate eq. of motion for galaxies
            CALL EqMovT1H   ! integrate eq. of motion for the test particles
            CALL EqMovTgas1H
            CALL StickyParticles
          ENDIF


!          IF (MOD(IT,10) .EQ. 0) THEN
!            CALL CompSurfAH               ! compute surf. density analytically 
!            CALL CompSurfN (0.d0, 0.d0)   ! compute surf. density numerically
                                        !   (test particles)
!            CALL MAPSURF (Display, IT*TUNIT, Drel, Vrel*VUNIT)
!          ENDIF

 300   CONTINUE
      ENDIF 


      CALL CPU_TIME(TIMEF)
      print*, 'CPU time [min] = ', (TIMEF - TIMEI) / 60.  

      

!      CALL PGEND
  
      END


*************************************************************************
      SUBROUTINE MAPSURF (Display, T, D, V)
*************************************************************************
* calls mapping subroutines:                                                      
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER Display*3

      IF (Display .EQ. 'win') THEN
        CALL MapSurf1 ('win', T, D, V)
!       CALL MapSurf2 ('win', T, D, V)
!       CALL MapSurf3 ('win', T, D, V)

      ELSEIF (Display .EQ. 'gif') THEN

          CALL MapSurf1 ('gif', T, D, V)  ! surface densities
          CALL MapSurf3 ('gif', T, D, V)  ! particle plots + pV diagram
          CALL MapSurf4 ('gif', T, D, V)  ! projections
          CALL MapSurf5 ('gif', T, D, V)  ! projections (zoom)
          CALL MapSurf6 ('gif', T, D, V)  ! unsharp masking 
          CALL MapSurf7 ('gif', T, D, V)  ! unsharp masking (zoom)

      ELSE

          CALL MapSurf1 ('cps', T, D, V)  ! surface densities
          CALL MapSurf3 ('cps', T, D, V)  ! particle plots + pV diagram
          CALL MapSurf4 ('cps', T, D, V)  ! projections
          CALL MapSurf5 ('cps', T, D, V)  ! projections (zoom)
          CALL MapSurf6 ('cps', T, D, V)  ! unsharp masking 
          CALL MapSurf7 ('cps', T, D, V)  ! unsharp masking (zoom)


      ENDIF

      END

*************************************************************************
      SUBROUTINE InitPGmap 
*************************************************************************
*     initialize PGPLOT

      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'merge.inp'

      CHARACTER PLOTDEV*20
      REAL*4 TR(6), TRzoom(6), XMAXext, XMAXint, VMaxExt, VMaxInt
      REAL*4 DUMMY2D(ND,ND)

      COMMON /TR6/    TR, TRzoom, XMAXext, XMAXint, VMaxExt, VMaxInt 
      COMMON /PANELS/ NCOLS, NROWS, NCRS
      COMMON /DUMMY/  DUMMY2D
      SAVE   /TR6/, /PANELS/, /DUMMY/

*-------------------------------
      DO 10 IY = 1,ND               ! creat a dummy 2D array for clearing 
        DO 10 IX = 1, ND            !   the screen between successive particle
 10       DUMMY2D (IX,IY) = 0.e0    !   plots:
*-------------------------------

      XMAXext = SNGL(XR)
      XMAXint = SNGL(XRz)

      VMAXext = XMAXext                ! max V range for pV plots 
      VMAXint = XMaxInt

      TR(1) = SNGL(-CG2 / CG1);   TR(4) = TR(1)
      TR(2) = SNGL(1./CG1);       TR(6) = TR(2)

      TRzoom(1) = SNGL(-CG2z / CG1z);   TRzoom(4) = TRzoom(1)
      TRzoom(2) = SNGL(1./CG1z);        TRzoom(6) = TRzoom(2)

      PLOTDEV = '/xwin'

!     NCOLS = 2;   NROWS = 2;   NCRS = NCols * NRows  
      NCOLS = 3;   NROWS = 2;   NCRS = NCols * NRows  

      CALL PGBEGIN(0,PLOTDEV, NCOLS, NROWS)

!     CALL ENV3x2
!     CALL ENV3x2zoom  ! when both lower & upper pannels used for zomm
!     CALL ENV3x2Lzoom ! when lower pannels used for zomm
      CALL ENV3x2pV    ! for particle + pV plots 

      CALL PALETT(3, 1.0, 0.5)  ! 3 = heat

      END

*************************************************************************
      SUBROUTINE ENV3x2 
*************************************************************************
*     initialize PGPLOT

      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'merge.inp'

      REAL*4 TR(6), TRzoom(6), XMAXext, XMAXint, VMaxExt, VMaxInt

      COMMON /UNITS/  TUNIT,VUNIT,YVUNIT,GMUNIT,GMUNIT6
      COMMON /GXIES2/ GMp, Bp, Bp2, YBp2, CTrBp, Sigma0p, Vesc0p, Sr0p,
     &                GMs, Bs, Bs2, YBs2, CTrBs, Sigma0s, Vesc0s, Sr0s,
     &                GMeff, GMStoP, GMpart1,GMpart2, Eps2,Eps2s,Eps2p, 
     &                Omega23, Dini, Vini, Trmin, DTrBp, DTrBs, Eps,
     &                RtrBp, RtrBs, BgasP, BgasS, BgasP2, BgasS2,
     &                YBgasP2, YBgasS2, CTrBgasP, CtrBgasS, Vesc0gasP,
     &                Vesc0gasS, Sr0gasP, Sr0gasS, DTrBgasP, DTrBgasS,
     &                RTrBgasP, RTrBgasS, Eps2gasP, Eps2gasS 
      COMMON /TR6/    TR, TRzoom, XMAXext, XMAXint, VMaxExt, VMaxInt 
      SAVE   /TR6/
 
      CALL PGASK (.FALSE.)

* upper pannels:
      CALL PGENV (-XMAXext, XMAXext, -XMAXext, XMAXext, 1, 0)     
      CALL LabelMap0 (GMp*GMUNIT, GMs/GMp, Bp, Bs/Bp, 
     &                Sr0p*VUNIT, Sr0s*VUNIT,
     &                Sigma0p,    Sigma0s,  Vesc0p*VUNIT,  
     &                Omega23, Dini,  Vini*VUNIT, Trmin)

      CALL PGENV (-XMAXext, XMAXext, -XMAXext, XMAXext, 1, 0)       
      CALL PGENV (-XMAXext, XMAXext, -XMAXext, XMAXext, 1, 0)     

* lower pannels:
      CALL PGENV (-XMAXext, XMAXext, -XMAXext, XMAXext, 1, 0)
      CALL PGENV (-XMAXext, XMAXext, -XMAXext, XMAXext, 1, 0) 
      CALL PGENV (-XMAXext, XMAXext, -XMAXext, XMAXext, 1, 0)       

      END


*************************************************************************
      SUBROUTINE ENV3x2Lzoom 
*************************************************************************
*     initialize PGPLOT

      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'merge.inp'

      REAL*4 TR(6), TRzoom(6), XMAXext, XMAXint, VMaxExt, VMaxInt

      COMMON /UNITS/  TUNIT,VUNIT,YVUNIT,GMUNIT,GMUNIT6
      COMMON /GXIES2/ GMp, Bp, Bp2, YBp2, CTrBp, Sigma0p, Vesc0p, Sr0p,
     &                GMs, Bs, Bs2, YBs2, CTrBs, Sigma0s, Vesc0s, Sr0s,
     &                GMeff, GMStoP, GMpart1,GMpart2, Eps2,Eps2s,Eps2p, 
     &                Omega23, Dini, Vini, Trmin, DTrBp, DTrBs, Eps,
     &                RtrBp, RtrBs, BgasP, BgasS, BgasP2, BgasS2,
     &                YBgasP2, YBgasS2, CTrBgasP, CtrBgasS, Vesc0gasP,
     &                Vesc0gasS, Sr0gasP, Sr0gasS, DTrBgasP, DTrBgasS,
     &                RTrBgasP, RTrBgasS, Eps2gasP, Eps2gasS    
      COMMON /TR6/    TR, TRzoom, XMAXext, XMAXint, VMaxExt, VMaxInt 
      SAVE   /TR6/
 
      CALL PGASK (.FALSE.)

* upper pannels:
      CALL PGENV (-XMAXext, XMAXext, -XMAXext, XMAXext, 1, 0)     
      CALL LabelMap0 (GMp*GMUNIT, GMs/GMp, Bp, Bs/Bp, 
     &                Sr0p*VUNIT, Sr0s*VUNIT,
     &                Sigma0p,    Sigma0s,  Vesc0p*VUNIT,  
     &                Omega23, Dini,  Vini*VUNIT, Trmin)

      CALL PGENV (-XMAXext, XMAXext, -XMAXext, XMAXext, 1, 0)       
      CALL PGENV (-XMAXext, XMAXext, -XMAXext, XMAXext, 1, 0)     

* lower pannels:
      CALL PGENV (-XMAXint, XMAXint, -XMAXint, XMAXint, 1, 0)
      CALL PGENV (-XMAXint, XMAXint, -XMAXint, XMAXint, 1, 0)
      CALL PGENV (-XMAXint, XMAXint, -XMAXint, XMAXint, 1, 0)

      END

*************************************************************************
      SUBROUTINE ENV3x2pV
*************************************************************************
*     initialize PGPLOT

      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'merge.inp'

      REAL*4 TR(6), TRzoom(6), XMAXext, XMAXint, VMaxExt, VMaxInt

      COMMON /UNITS/  TUNIT,VUNIT,YVUNIT,GMUNIT,GMUNIT6
      COMMON /GXIES2/ GMp, Bp, Bp2, YBp2, CTrBp, Sigma0p, Vesc0p, Sr0p,
     &                GMs, Bs, Bs2, YBs2, CTrBs, Sigma0s, Vesc0s, Sr0s,
     &                GMeff, GMStoP, GMpart1,GMpart2, Eps2,Eps2s,Eps2p, 
     &                Omega23, Dini, Vini, Trmin, DTrBp, DTrBs, Eps,
     &                RtrBp, RtrBs, BgasP, BgasS, BgasP2, BgasS2,
     &                YBgasP2, YBgasS2, CTrBgasP, CtrBgasS, Vesc0gasP,
     &                Vesc0gasS, Sr0gasP, Sr0gasS, DTrBgasP, DTrBgasS,
     &                RTrBgasP, RTrBgasS, Eps2gasP, Eps2gasS   
      COMMON /TR6/    TR, TRzoom, XMAXext, XMAXint, VMaxExt, VMaxInt 
      SAVE   /TR6/
 
      CALL PGASK (.FALSE.)

* upper pannels:
      CALL PGENV (-XMAXext, XMAXext, -XMAXext, XMAXext, 1, 0)     
      CALL LabelMap0 (GMp*GMUNIT, GMs/GMp, Bp, Bs/Bp, 
     &                Sr0p*VUNIT, Sr0s*VUNIT,
     &                Sigma0p,    Sigma0s,  Vesc0p*VUNIT,  
     &                Omega23, Dini,  Vini*VUNIT, Trmin)

      CALL PGENV (-XMAXext, XMAXext, -XMAXext, XMAXext, 1, 0)       
      CALL PGENV (-XMAXext, XMAXext, -VMAXext, VMAXext, 1, 0)       

* lower pannels:
      CALL PGENV (-XMAXint, XMAXint, -XMAXint, XMAXint, 1, 0)
      CALL PGENV (-XMAXint, XMAXint, -XMAXint, XMAXint, 1, 0)
      CALL PGENV (-XMAXint, XMAXint, -VMAXint, VMAXint, 1, 0)

      END

*************************************************************************
      SUBROUTINE ENV3x2zoom 
*************************************************************************
*     initialize PGPLOT

      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'merge.inp'

      REAL*4 TR(6), TRzoom(6), XMAXext, XMAXint, VMaxExt, VMaxInt

      COMMON /UNITS/  TUNIT,VUNIT,YVUNIT,GMUNIT,GMUNIT6
      COMMON /GXIES2/ GMp, Bp, Bp2, YBp2, CTrBp, Sigma0p, Vesc0p, Sr0p,
     &                GMs, Bs, Bs2, YBs2, CTrBs, Sigma0s, Vesc0s, Sr0s,
     &                GMeff, GMStoP, GMpart1,GMpart2, Eps2,Eps2s,Eps2p, 
     &                Omega23, Dini, Vini, Trmin, DTrBp, DTrBs, Eps,
     &                RtrBp, RtrBs, BgasP, BgasS, BgasP2, BgasS2,
     &                YBgasP2, YBgasS2, CTrBgasP, CtrBgasS, Vesc0gasP,
     &                Vesc0gasS, Sr0gasP, Sr0gasS, DTrBgasP, DTrBgasS,
     &                RTrBgasP, RTrBgasS, Eps2gasP, Eps2gasS   
      COMMON /TR6/    TR, TRzoom, XMAXext, XMAXint, VMaxExt, VMaxInt 
      SAVE   /TR6/
 
      CALL PGASK (.FALSE.)

* upper pannels:
      CALL PGENV (-XMAXint, XMAXint, -XMAXint, XMAXint, 1, 0)     
      CALL LabelMap0 (GMp*GMUNIT, GMs/GMp, Bp, Bs/Bp, 
     &                Sr0p*VUNIT, Sr0s*VUNIT,
     &                Sigma0p,    Sigma0s,  Vesc0p*VUNIT,  
     &                Omega23, Dini,  Vini*VUNIT, Trmin)

      CALL PGENV (-XMAXint, XMAXint, -XMAXint, XMAXint, 1, 0)       
      CALL PGENV (-XMAXint, XMAXint, -XMAXint, XMAXint, 1, 0)     

* lower pannels:
      CALL PGENV (-XMAXint, XMAXint, -XMAXint, XMAXint, 1, 0)
      CALL PGENV (-XMAXint, XMAXint, -XMAXint, XMAXint, 1, 0)
      CALL PGENV (-XMAXint, XMAXint, -XMAXint, XMAXint, 1, 0)

      END


*************************************************************************
      SUBROUTINE CompSurfAP  
*************************************************************************
*     computes surface density of the primary and secondary (analytically)
*     valid only for the motion in the (x,y) plane - 
*      --> remake for the general case

      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'merge.inp'

      DIMENSION XG(2)
      COMMON /GXIES/  Xp(3),Xs(3),Vp(3),Vs(3)
      COMMON /GXIES2/ GMp, Bp, Bp2, YBp2, CTrBp, Sigma0p, Vesc0p, Sr0p,
     &                GMs, Bs, Bs2, YBs2, CTrBs, Sigma0s, Vesc0s, Sr0s,
     &                GMeff, GMStoP, GMpart1,GMpart2, Eps2,Eps2s,Eps2p, 
     &                Omega23, Dini, Vini, Trmin, DTrBp, DTrBs, Eps,
     &                RtrBp, RtrBs, BgasP, BgasS, BgasP2, BgasS2,
     &                YBgasP2, YBgasS2, CTrBgasP, CtrBgasS, Vesc0gasP,
     &                Vesc0gasS, Sr0gasP, Sr0gasS, DTrBgasP, DTrBgasS,
     &                RTrBgasP, RTrBgasS, Eps2gasP, Eps2gasS   
      COMMON /DENS2/  RHOpA (ND,ND), RHOpN (ND,ND), 
     &                RHOsA (ND,ND), RHOsN (ND,ND), 
     &                RHOtA (ND,ND), RHOtN (ND,ND), Ratio(ND,ND)
      COMMON /DENS2z/ RHOpAz(NDz,NDz), RHOpNz (NDz,NDz), ! zoom grid
     &                RHOsAz(NDz,NDz), RHOsNz (NDz,NDz), 
     &                RHOtAz(NDz,NDz), RHOtNz (NDz,NDz), Ratioz(NDz,NDz)
      SAVE /DENS2/, /DENS2z/ 
                      ! RHO[ps]A - (anal.) surf. dens. of the primary/secondary
                      ! RHO[ps]N - surface density of the test particles
                      ! RHOtA = RHOpA + RHOsA 
                      ! RHOtN = RHOpN + RHOsN 


! large grid
      DO 10 IY = 1,ND
        DO 10 IX = 1,ND
          XG(1) = DX * ( IX - 0.5d0*ND1 ) ! projected grid (observer) coor.
          XG(2) = DX * ( IY - 0.5d0*ND1 )
        
          R2 = (XG(1) - Xp(1))**2 + (XG(2) - Xp(2))**2 
          RHOpA (IX,IY) = PlummerSD (Sigma0p, YBp2, R2)  

          R2 = (XG(1) - Xs(1))**2 + (XG(2) - Xs(2))**2 
          RHOsA (IX,IY) = PlummerSD (Sigma0s, YBs2, R2)  

          RHOtA (IX,IY) = RHOpA (IX,IY) + RHOsA (IX,IY)
 10   CONTINUE    


! zoom grid
      DO 20 IY = 1,NDz
        DO 20 IX = 1,NDz
          XG(1) = DXz * ( IX - 0.5d0*ND1z ) ! projected grid (observer) coor.
          XG(2) = DXz * ( IY - 0.5d0*ND1z )
        
          R2 = (XG(1) - Xp(1))**2 + (XG(2) - Xp(2))**2 
          RHOpAz (IX,IY) = PlummerSD (Sigma0p, YBp2, R2)  

          R2 = (XG(1) - Xs(1))**2 + (XG(2) - Xs(2))**2 
          RHOsAz (IX,IY) = PlummerSD (Sigma0s, YBs2, R2)  

          RHOtAz (IX,IY) = RHOpAz (IX,IY) + RHOsAz (IX,IY)
 20   CONTINUE    

      END

*************************************************************************
      SUBROUTINE CompSurfAH  
*************************************************************************
*     computes surface density of the primary and secondary (analytically)
*     valid only for the motion in the (x,y) plane - 
*      --> remake for the general case

      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'merge.inp'

      DIMENSION XG(2)
      COMMON /GXIES/  Xp(3),Xs(3),Vp(3),Vs(3)
      COMMON /GXIES2/ GMp, Bp, Bp2, YBp2, CTrBp, Sigma0p, Vesc0p, Sr0p,
     &                GMs, Bs, Bs2, YBs2, CTrBs, Sigma0s, Vesc0s, Sr0s,
     &                GMeff, GMStoP, GMpart1,GMpart2, Eps2,Eps2s,Eps2p, 
     &                Omega23, Dini, Vini, Trmin, DTrBp, DTrBs, Eps,
     &                RtrBp, RtrBs, BgasP, BgasS, BgasP2, BgasS2,
     &                YBgasP2, YBgasS2, CTrBgasP, CtrBgasS, Vesc0gasP,
     &                Vesc0gasS, Sr0gasP, Sr0gasS, DTrBgasP, DTrBgasS,
     &                RTrBgasP, RTrBgasS, Eps2gasP, Eps2gasS   
      COMMON /DENS2/  RHOpA (ND,ND), RHOpN (ND,ND), 
     &                RHOsA (ND,ND), RHOsN (ND,ND), 
     &                RHOtA (ND,ND), RHOtN (ND,ND), Ratio(ND,ND)
      COMMON /DENS2z/ RHOpAz(NDz,NDz), RHOpNz (NDz,NDz), ! zoom grid
     &                RHOsAz(NDz,NDz), RHOsNz (NDz,NDz), 
     &                RHOtAz(NDz,NDz), RHOtNz (NDz,NDz), Ratioz(NDz,NDz)
      SAVE /DENS2/, /DENS2z/ 
                      ! RHO[ps]A - (anal.) surf. dens. of the primary/secondary
                      ! RHO[ps]N - surface density of the test particles
                      ! RHOtA = RHOpA + RHOsA 
                      ! RHOtN = RHOpN + RHOsN 


! large grid
      DO 10 IY = 1,ND
        DO 10 IX = 1,ND
          XG(1) = DX * ( IX - 0.5d0*ND1 ) ! projected grid (observer) coor.
          XG(2) = DX * ( IY - 0.5d0*ND1 )
        
          R = DSQRT((XG(1) - Xp(1))**2 + (XG(2) - Xp(2))**2)
          RHOpA (IX,IY) = HernquistSD (Sigma0p, Bp, R)  

          R = DSQRT((XG(1) - Xs(1))**2 + (XG(2) - Xs(2))**2)
          RHOsA (IX,IY) = HernquistSD (Sigma0s, Bs, R)  

          RHOtA (IX,IY) = RHOpA (IX,IY) + RHOsA (IX,IY)
 10   CONTINUE    


! zoom grid
      DO 20 IY = 1,NDz
        DO 20 IX = 1,NDz
          XG(1) = DXz * ( IX - 0.5d0*ND1z ) ! projected grid (observer) coor.
          XG(2) = DXz * ( IY - 0.5d0*ND1z )
        
          R = DSQRT((XG(1) - Xp(1))**2 + (XG(2) - Xp(2))**2) 
          RHOpAz (IX,IY) = HernquistSD (Sigma0p, Bp, R)  

          R = DSQRT((XG(1) - Xs(1))**2 + (XG(2) - Xs(2))**2)
          RHOsAz (IX,IY) = HernquistSD (Sigma0s, Bs, R)  

          RHOtAz (IX,IY) = RHOpAz (IX,IY) + RHOsAz (IX,IY)
 20   CONTINUE    

      END
      
*************************************************************************
      SUBROUTINE CompSurfAV  
*************************************************************************
*     computes surface density of the primary and secondary (analytically)
*     valid only for the motion in the (x,y) plane - 
*      --> remake for the general case
* like ComSurfAP but Primary = de Vaucouleurs

      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'merge.inp'

      DIMENSION XG(2)
      COMMON /GXIES/  Xp(3),Xs(3),Vp(3),Vs(3)
      COMMON /GXIES2/ GMp, Bp, Bp2, YBp2, CTrBp, Sigma0p, Vesc0p, Sr0p,
     &                GMs, Bs, Bs2, YBs2, CTrBs, Sigma0s, Vesc0s, Sr0s,
     &                GMeff, GMStoP, GMpart1,GMpart2, Eps2,Eps2s,Eps2p, 
     &                Omega23, Dini, Vini, Trmin, DTrBp, DTrBs, Eps,
     &                RtrBp, RtrBs, BgasP, BgasS, BgasP2, BgasS2,
     &                YBgasP2, YBgasS2, CTrBgasP, CtrBgasS, Vesc0gasP,
     &                Vesc0gasS, Sr0gasP, Sr0gasS, DTrBgasP, DTrBgasS,
     &                RTrBgasP, RTrBgasS, Eps2gasP, Eps2gasS   
      COMMON /DENS2/  RHOpA (ND,ND), RHOpN (ND,ND), 
     &                RHOsA (ND,ND), RHOsN (ND,ND), 
     &                RHOtA (ND,ND), RHOtN (ND,ND), Ratio(ND,ND)
      COMMON /DENS2z/ RHOpAz(NDz,NDz), RHOpNz (NDz,NDz), ! zoom grid
     &                RHOsAz(NDz,NDz), RHOsNz (NDz,NDz), 
     &                RHOtAz(NDz,NDz), RHOtNz (NDz,NDz), Ratioz(NDz,NDz)
      SAVE /DENS2/, /DENS2z/ 
                      ! RHO[ps]A - (anal.) surf. dens. of the primary/secondary
                      ! RHO[ps]N - surface density of the test particles
                      ! RHOtA = RHOpA + RHOsA 
                      ! RHOtN = RHOpN + RHOsN 


! large grid
      DO 10 IY = 1,ND
        DO 10 IX = 1,ND
          XG(1) = DX * ( IX - 0.5d0*ND1 ) ! projected grid (observer) coor.
          XG(2) = DX * ( IY - 0.5d0*ND1 )
        
          R2 = (XG(1) - Xp(1))**2 + (XG(2) - Xp(2))**2 
          RHOpA (IX,IY) = DeVacSD (Sigma0p, YBp2, R2)  

          R2 = (XG(1) - Xs(1))**2 + (XG(2) - Xs(2))**2 
          RHOsA (IX,IY) = PlummerSD (Sigma0s, YBs2, R2)  

          RHOtA (IX,IY) = RHOpA (IX,IY) + RHOsA (IX,IY)
 10   CONTINUE    


! zoom grid
      DO 20 IY = 1,NDz
        DO 20 IX = 1,NDz
          XG(1) = DXz * ( IX - 0.5d0*ND1z ) ! projected grid (observer) coor.
          XG(2) = DXz * ( IY - 0.5d0*ND1z )
        
          R2 = (XG(1) - Xp(1))**2 + (XG(2) - Xp(2))**2 
          RHOpAz (IX,IY) = DeVacSD (Sigma0p, YBp2, R2)  

          R2 = (XG(1) - Xs(1))**2 + (XG(2) - Xs(2))**2 
          RHOsAz (IX,IY) = PlummerSD (Sigma0s, YBs2, R2)  

          RHOtAz (IX,IY) = RHOpAz (IX,IY) + RHOsAz (IX,IY)
 20   CONTINUE    

      END


*************************************************************************
      SUBROUTINE CompSurfN (PAdeg, AIdeg)
*************************************************************************
*     computes surface density of test particles
*     valid only for the motion in the (x,y) plane - 
*      --> remake for the general case

      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'merge.inp'

      COMMON /GXIES2/ GMp, Bp, Bp2, YBp2, CTrBp, Sigma0p, Vesc0p, Sr0p,
     &                GMs, Bs, Bs2, YBs2, CTrBs, Sigma0s, Vesc0s, Sr0s,
     &                GMeff, GMStoP, GMpart1,GMpart2, Eps2,Eps2s,Eps2p, 
     &                Omega23, Dini, Vini, Trmin, DTrBp, DTrBs, Eps,
     &                RtrBp, RtrBs, BgasP, BgasS, BgasP2, BgasS2,
     &                YBgasP2, YBgasS2, CTrBgasP, CtrBgasS, Vesc0gasP,
     &                Vesc0gasS, Sr0gasP, Sr0gasS, DTrBgasP, DTrBgasS,
     &                RTrBgasP, RTrBgasS, Eps2gasP, Eps2gasS   
      COMMON /TestP1/ N1, N2, NALL       
      COMMON /DENS2/  RHOpA (ND,ND), RHOpN (ND,ND), 
     &                RHOsA (ND,ND), RHOsN (ND,ND), 
     &                RHOtA (ND,ND), RHOtN (ND,ND), Ratio(ND,ND)
      COMMON /DENS2z/ RHOpAz(NDz,NDz), RHOpNz (NDz,NDz), ! zoom grid
     &                RHOsAz(NDz,NDz), RHOsNz (NDz,NDz), 
     &                RHOtAz(NDz,NDz), RHOtNz (NDz,NDz), Ratioz(NDz,NDz)
      SAVE /DENS2/, /DENS2z/ 
                      ! RHO[ps]A - (anal.) surf. dens. of the primary/secondary
                      ! RHO[ps]N - surface density of the test particles
                      ! RHOtA = RHOpA + RHOsA 
                      ! RHOtN = RHOpN + RHOsN 

      CALL ProjectXYZ (0, NALL, PAdeg, AIdeg)

! large grid: comp. surf.dens. of test part.(primary) 
      CALL RHO1CIC (0,  N1, RHOpN, ND, CG1, CG2, XR) 
!             comp. surf.dens. of test part.(secondary)
      CALL RHO1CIC (N1, N2, RHOsN, ND, CG1, CG2, XR) 

! zoom grid: comp. surf.dens. of test part.(primary) 
      CALL RHO1CIC (0,  N1, RHOpNz, NDz, CG1z, CG2z, XRz)
!            comp. surf.dens. of test part.(secondary)
      CALL RHO1CIC (N1, N2, RHOsNz, NDz, CG1z, CG2z, XRz)


! large grid
      DO 10 IY = 1,ND
        DO 10 IX = 1,ND
          RHOpN (IX,IY) = GMpart1 * RHOpN (IX,IY) 
          RHOsN (IX,IY) = GMpart2 * RHOsN (IX,IY) 
!         RHOtN (IX,IY) = RHOpN (IX,IY) + RHOsN (IX,IY)
          RHOtN (IX,IY) = RHOpA (IX,IY) + RHOsN (IX,IY)
          Ratio (IX,IY) = RhosN (IX,IY) / RHOpA (IX,IY)
 10   CONTINUE    


! zoom grid
      DO 20 IY = 1,NDz
        DO 20 IX = 1,NDz
          RHOpNz (IX,IY) = GMpart1 * ZoomRatio * RHOpNz (IX,IY) 
          RHOsNz (IX,IY) = GMpart2 * ZoomRatio * RHOsNz (IX,IY) 
!         RHOtNz (IX,IY) = RHOpNz (IX,IY) + RHOsNz (IX,IY)
          RHOtNz (IX,IY) = RHOpAz (IX,IY) + RHOsNz (IX,IY)
          Ratioz (IX,IY) = RhosNz (IX,IY) / RHOpAz (IX,IY)
 20   CONTINUE    
 

      END


*************************************************************************
      SUBROUTINE ProjectXYZ (Nother, N, PA0, AI0)
*************************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'merge.inp'

      PARAMETER (PI=3.141592653589, DegToRad=PI/180., M = 3)
      DIMENSION R1(M,M), R2(M,M), R3(M,M), R(M,M)
      COMMON /TestP/  X(NX), Y(NX), Z(NX), VX(NX), VY(NX), VZ(NX)
      COMMON /TestP3/ Xpro (NX), Ypro (NX)
      SAVE /TestP3/

      PA = PA0 * DegToRad;   AI = AI0 * DegToRad
      COSPA = DCOS(PA);      SINPA = DSIN(PA)
      COSAI = DCOS(AI);      SINAI = DSIN(AI) 

      R1(1,1) =  COSPA;   R1(1,2) =  SINPA;   R1(1,3) =   0.d0
      R1(2,1) = -SINPA;   R1(2,2) =  COSPA;   R1(2,3) =   0.d0
      R1(3,1) =   0.d0;   R1(3,2) =   0.d0;   R1(3,3) =   1.d0

      R2(1,1) =  COSAI;   R2(1,2) =   0.d0;   R2(1,3) = -SINAI
      R2(2,1) =   0.d0;   R2(2,2) =   1.d0;   R2(2,3) =   0.d0
      R2(3,1) =  SINAI;   R2(3,2) =   0.d0;   R2(3,3) =  COSAI
      
      R3(1,1) =  COSPA;   R3(1,2) = -SINPA;   R3(1,3) =   0.d0
      R3(2,1) =  SINPA;   R3(2,2) =  COSPA;   R3(2,3) =   0.d0
      R3(3,1) =   0.d0;   R3(3,2) =   0.d0;   R3(3,3) =   1.d0

      CALL MPRODUCT(R2, R1, R , M)
      CALL MPRODUCT(R3, R , R1, M) 

      DO 10 I = Nother + 1, Nother + N
        Xpro(I)   = R1(1,1) * X(I) + R1(1,2) * Y(I) + R1(1,3) * Z(I) 
        Ypro(I)   = R1(2,1) * X(I) + R1(2,2) * Y(I) + R1(2,3) * Z(I) 
 10   CONTINUE

      END  

*************************************************************************
      SUBROUTINE RHO1CIC (Nother, N, RHO, ND0, CG10, CG20, XR0)
*************************************************************************
*     computes surface density on a grid via CIC scheme

      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'merge.inp'

      DIMENSION RHO (ND0,ND0)
      COMMON /TestP3/ X(NX), Y(NX)
      SAVE /TestP3/

      DO 10, IY = 1, ND0
       DO 10, IX = 1, ND0
10       RHO (IX,IY) = 0.d0

      DO 100, I = Nother + 1, Nother + N

        IF (DABS(X(I)) .LT. XR0 .AND. DABS(Y(I)) .LT. XR0) THEN 
!       IF ((X(I)**2 + Y(I)**2) .LT. XR20) THEN        ! particle inside grid
          XX=X(I)*CG10+CG20;    YY=Y(I)*CG10+CG20

          JX=XX;              JY=YY
          IX=JX+1;            IY=JY+1
          XI=XX-JX;           YI=YY-JY
          XJ=IX-XX;           YJ=IY-YY

          RHO (IX,IY) = RHO (IX,IY) + XI*YI
          RHO (IX,JY) = RHO (IX,JY) + XI*YJ
          RHO (JX,IY) = RHO (JX,IY) + XJ*YI
          RHO (JX,JY) = RHO (JX,JY) + XJ*YJ
        ENDIF

100   CONTINUE

      END

*************************************************************************
      REAL*8 FUNCTION PlummerSD (Sigma0, YB2, R2) 
*************************************************************************
! analytically computes surface density for a Plummer sphere

      IMPLICIT REAL*8 (A-H,O-Z)

      PlummerSD = Sigma0 / (1.d0 + R2 * YB2)**2

      END 

*************************************************************************
      REAL*8 FUNCTION HernquistSD (Sigma0, B, R) 
*************************************************************************
! analytically computes surface density for a Hernquist sphere

      IMPLICIT REAL*8 (A-H,O-Z)

      S = R/B
      
      IF (S.GE.1) THEN
        Xfn = ACOS(1.d0/S) / DSQRT(S**2.d0 - 1.d0)
        HernquistSD=Sigma0*((2.d0+S**2.) *Xfn-3.d0) /(1.d0-S**2.)**2.
      ELSE
        Xfn = LOG((1.d0+ DSQRT(1.d0-S**2.d0))/S) / DSQRT(1.d0-S**2.d0)
        HernquistSD=Sigma0*((2.d0+S**2.) *Xfn-3.d0) /(1.d0-S**2.)**2.
      ENDIF

      END 
      
*************************************************************************
      REAL*8 FUNCTION DeVacSD (Sigma0, YB2, R2) 
*************************************************************************
! analytically computes surface density for de Vaucouleurs law

      IMPLICIT REAL*8 (A-H,O-Z)

      S = (R2 * YB2)**0.125        ! (R/R_eff)^1/4
      DeVacSD = Sigma0 * DEXP(-7.67d0 * (S-1.d0))

      END 

*************************************************************************
      REAL*8 FUNCTION GenPosPL0 (Xrandom, BPlummer, CPlummer) 
*************************************************************************
! generates particle position within Plummer sphere 
!    of scale-length BPlummer; 
! the sphere can be truncated, CPlummer is the ratio of the total mass
! (of the untruncated) and of the truncated mass (if CPlummer = 1, the
! the sphere is untruncated  

      IMPLICIT REAL*8 (A-H,O-Z)

      GenPosPL0 = BPlummer / DSQRT ((CPlummer / Xrandom)**(2./3.) - 1.)  

      END

*************************************************************************
      REAL*8 FUNCTION GenPosHS0 (Xrandom, Rhs) 
*************************************************************************
! generates particle position within a homoheneous sphere of radius Rhs 

      IMPLICIT REAL*8 (A-H,O-Z)

      GenPosHS0 = Rhs * Xrandom**(1./3.)                            

      END 

*************************************************************************
      REAL*8 FUNCTION GenPosHQ0 (Xrandom, BHernquist, CHernquist) 
*************************************************************************
! generates particle position within Hernquist sphere 
!    of scale-length BHernquist; 
! the sphere can be truncated, CHernquist is the ratio of the total mass
! (of the untruncated) and of the truncated mass (if CHernquist = 1, the
! the sphere is untruncated  

      IMPLICIT REAL*8 (A-H,O-Z)

      GenPosHQ0 = BHernquist / ((CHernquist / Xrandom)**(1./2.) - 1.)  

      END

*************************************************************************
      SUBROUTINE MapSurf1 (Display, TMyr, Drel, Vrel)
*************************************************************************
*     visualizes surface density
*     display = 3x2 pannels:
*                   P+S   S     S/P
*                   zoom  zoom  zoom      


      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'merge.inp'

      REAL*4 TR(6), TRzoom(6), XMAXext, XMAXint, VMaxExt, VMaxInt
      REAL*4 ArrMin, ArrMax1, ArrMax2, ArrMax3
      CHARACTER GifFile*9, Display*3

      COMMON /GXIES2/ GMp, Bp, Bp2, YBp2, CTrBp, Sigma0p, Vesc0p, Sr0p,
     &                GMs, Bs, Bs2, YBs2, CTrBs, Sigma0s, Vesc0s, Sr0s,
     &                GMeff, GMStoP, GMpart1,GMpart2, Eps2,Eps2s,Eps2p, 
     &                Omega23, Dini, Vini, Trmin, DTrBp, DTrBs, Eps,
     &                RtrBp, RtrBs, BgasP, BgasS, BgasP2, BgasS2,
     &                YBgasP2, YBgasS2, CTrBgasP, CtrBgasS, Vesc0gasP,
     &                Vesc0gasS, Sr0gasP, Sr0gasS, DTrBgasP, DTrBgasS,
     &                RTrBgasP, RTrBgasS, Eps2gasP, Eps2gasS   
      COMMON /DENS2/  RHOpA (ND,ND), RHOpN (ND,ND), 
     &                RHOsA (ND,ND), RHOsN (ND,ND), 
     &                RHOtA (ND,ND), RHOtN (ND,ND), Ratio(ND,ND)
      COMMON /DENS2z/ RHOpAz(NDz,NDz), RHOpNz (NDz,NDz), ! zoom grid
     &                RHOsAz(NDz,NDz), RHOsNz (NDz,NDz), 
     &                RHOtAz(NDz,NDz), RHOtNz (NDz,NDz), Ratioz(NDz,NDz)
      SAVE   /DENS2/, /DENS2z/
                      ! RHOp  - surface density of the primary
                      ! RHOsA - (analytical) surface density of the secondary
                      ! RHOsN - surface density of the test particles
                      ! RHOtA = RHOp + RHOsA 
                      ! RHOtN = RHOp + RHOsN 
      COMMON /TR6/    TR, TRzoom, XMAXext, XMAXint, VMaxExt, VMaxInt 
      COMMON /PANELS/ NCOLS, NROWS, NCRS
      SAVE   /TR6/,/PANELS/

      DATA NoFile1 /0/
      SAVE NoFile1    

      ArrMin = 0.d0;   ArrMax1 = SNGL(1.*Sigma0p);  
                     ! ArrMax2 = 5.e-3*ArrMax1
                       ArrMax2 = SNGL(0.002*Sigma0s)
                       ArrMax3 = 0.1

      IF (Display .EQ. 'gif') THEN
        CALL FILENAME3 (NoFile1, 'sh1-', 'gif', GifFile)
        CALL PGBEGIN(0,GifFile//'.gif/gif', NCOLS, NROWS)
        CALL ENV3x2Lzoom
        CALL PALETT(3, 1.0, 0.5)  ! 3 = heat
      ELSEIF (Display .EQ. 'cps') THEN
        CALL FILENAME3 (NoFile1, 'sh1-', 'gif', GifFile)
        CALL PGBEGIN(0,GifFile//'.cps/cps', NCOLS, NROWS)
        CALL ENV3x2Lzoom
        CALL PALETT(3, 1.0, 0.5)  ! 3 = heat
      ENDIF

* upper pannels:

      CALL FindPanel;  CALL PGSWIN(-XMaxExt, XMaxExt, -XMaxExt, XMaxExt)  
      CALL MAPRHO (RHOtN, ND, TR, ArrMin, ArrMax1, 1) ! map primary
      CALL LabelMap (TMyr, Drel, Vrel)

      CALL FindPanel;  CALL PGSWIN(-XMaxExt, XMaxExt, -XMaxExt, XMaxExt)  
      CALL MAPRHO (RHOsN, ND, TR, ArrMin, ArrMax2, 1) 

      CALL FindPanel;  CALL PGSWIN(-XMaxExt, XMaxExt, -XMaxExt, XMaxExt)  
      CALL MAPRHO (Ratio, ND, TR, ArrMin, ArrMax3, 1) ! Flux Ratio

* lower pannels:
      CALL FindPanel;  CALL PGSWIN(-XMaxInt, XMaxInt, -XMaxInt, XMaxInt)  
      CALL MAPRHO (RHOtNz, NDz, TRzoom, ArrMin, ArrMax1, 1) ! zoom primary

      CALL FindPanel;  CALL PGSWIN(-XMaxInt, XMaxInt, -XMaxInt, XMaxInt)  
      CALL MAPRHO (RHOsNz, NDz, TRzoom, ArrMin, ArrMax2, 1) 
 
      CALL FindPanel;  CALL PGSWIN(-XMaxInt, XMaxInt, -XMaxInt, XMaxInt)  
      CALL MAPRHO (Ratioz, NDz, TRzoom, ArrMin, ArrMax3, 1) ! Flux Ratio


* save density to ascii files
c      CALL FILENAME3 (NoFile1, 'sh1-', 'rho', GifFile)
c      WRITE (20,1) ( (RHOtN(I,J), I=1,ND), J=1,ND)
c      WRITE (21,1) ( (RHOsN(I,J), I=1,ND), J=1,ND)
c      CLOSE(20);  CLOSE(21) 

      IF (Display .EQ. 'gif' .OR. Display .EQ. 'cps') THEN
        CALL PGEND
      ENDIF 

      NoFile1 = NoFile1 + 1       

  1   FORMAT (E11.4)

      END


*************************************************************************
      SUBROUTINE MapSurf2 (Display, TMyr, Drel, Vrel)
*************************************************************************
*     modification of MapSurf1gif
*     visualizes surface density and particle plots

*     display = 3x2 pannels:
*                   P+S (surf. dens.)   S (surf. dens.)  S (particles)
*                   zoom                zoom             zoom      

*     falling particles are blue, expanding red


      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'merge.inp'

      REAL*4 TR(6), TRzoom(6), XMAXext, XMAXint, VMaxExt, VMaxInt
      REAL*4 ArrMin, ArrMax1, ArrMax2, ArrMax3
      REAL*4 DUMMY2D(ND,ND)

      CHARACTER GifFile*9, Display*3

      COMMON /GXIES2/ GMp, Bp, Bp2, YBp2, CTrBp, Sigma0p, Vesc0p, Sr0p,
     &                GMs, Bs, Bs2, YBs2, CTrBs, Sigma0s, Vesc0s, Sr0s,
     &                GMeff, GMStoP, GMpart1,GMpart2, Eps2,Eps2s,Eps2p, 
     &                Omega23, Dini, Vini, Trmin, DTrBp, DTrBs, Eps,
     &                RtrBp, RtrBs, BgasP, BgasS, BgasP2, BgasS2,
     &                YBgasP2, YBgasS2, CTrBgasP, CtrBgasS, Vesc0gasP,
     &                Vesc0gasS, Sr0gasP, Sr0gasS, DTrBgasP, DTrBgasS,
     &                RTrBgasP, RTrBgasS, Eps2gasP, Eps2gasS   
      COMMON /DENS2/  RHOpA (ND,ND), RHOpN (ND,ND), 
     &                RHOsA (ND,ND), RHOsN (ND,ND), 
     &                RHOtA (ND,ND), RHOtN (ND,ND), Ratio(ND,ND)
      COMMON /DENS2z/ RHOpAz(NDz,NDz), RHOpNz (NDz,NDz), ! zoom grid
     &                RHOsAz(NDz,NDz), RHOsNz (NDz,NDz), 
     &                RHOtAz(NDz,NDz), RHOtNz (NDz,NDz), Ratioz(NDz,NDz)
      SAVE   /DENS2/, /DENS2z/
                      ! RHOp  - surface density of the primary
                      ! RHOsA - (analytical) surface density of the secondary
                      ! RHOsN - surface density of the test particles
                      ! RHOtA = RHOp + RHOsA 
                      ! RHOtN = RHOp + RHOsN
      COMMON /TR6/    TR, TRzoom, XMAXext, XMAXint, VMaxExt, VMaxInt 
      COMMON /PANELS/ NCOLS, NROWS, NCRS
      COMMON /TestP1/ N1, N2, NALL    
      COMMON /DUMMY/  DUMMY2D

      SAVE   /TR6/, /PANELS/, /DUMMY/
      DATA NoFile2 /0/
      SAVE NoFile2    

      ArrMin = 0.d0;   ArrMax1 = SNGL(1.*Sigma0p);  
                     ! ArrMax2 = 5.e-3*ArrMax1
                       ArrMax2 = SNGL(0.002*Sigma0s)
                       ArrMax3 = 0.1

      IF (Display .EQ. 'gif') THEN
        CALL FILENAME3 (NoFile2, 'sh2-', 'gif', GifFile)
        CALL PGBEGIN(0,GifFile//'.gif/gif', NCOLS, NROWS)
        CALL ENV3x2Lzoom
        CALL PALETT(3, 1.0, 0.5)  ! 3 = heat
      ELSEIF (Display .EQ. 'cps') THEN
        CALL FILENAME3 (NoFile2, 'sh2-', 'gif', GifFile)
        CALL PGBEGIN(0,GifFile//'.cps/cps', NCOLS, NROWS)
        CALL ENV3x2
        CALL PALETT(3, 1.0, 0.5)  ! 3 = heat
      ENDIF

* upper pannels:

      CALL FindPanel;  CALL PGSWIN(-XMaxExt, XMaxExt, -XMaxExt, XMaxExt)  
      CALL MAPRHO (RHOtN, ND, TR, ArrMin, ArrMax1, 1) ! map primary
      CALL LabelMap (TMyr, Drel, Vrel)

      CALL FindPanel;  CALL PGSWIN(-XMaxExt, XMaxExt, -XMaxExt, XMaxExt)  
      CALL MAPRHO (RHOsN, ND, TR, ArrMin, ArrMax2, 1) ! map primary

      CALL FindPanel;  CALL PGSWIN(-XMaxExt, XMaxExt, -XMaxExt, XMaxExt)  
      CALL PGIMAG(DUMMY2D, ND, ND, 1, ND, 1, ND, ArrMin, ArrMax1, TR)  
      CALL PARTPLOT2 (N1, N2, 4, 2) ! falling particles blue, expanding red

* lower pannels:
      CALL FindPanel;  CALL PGSWIN(-XMaxInt, XMaxInt, -XMaxInt, XMaxInt)  
      CALL MAPRHO (RHOtNz, NDz, TRzoom, ArrMin, ArrMax1, 1) ! zoom primary

      CALL FindPanel;  CALL PGSWIN(-XMaxInt, XMaxInt, -XMaxInt, XMaxInt)  
      CALL MAPRHO (RHOsNz, NDz, TRzoom, ArrMin, ArrMax2, 1) ! zoom primary

      CALL FindPanel;  CALL PGSWIN(-XMaxInt, XMaxInt, -XMaxInt, XMaxInt)  
      CALL PGIMAG(DUMMY2D, ND, ND, 1, ND, 1, ND, ArrMin, ArrMax1, TR)  
      CALL PARTPLOT2 (N1, N2, 4, 2) ! falling particles blue, expanding red

      IF (Display .EQ. 'gif' .OR. Display .EQ. 'cps') THEN
        NoFile2 = NoFile2 + 1       
        CALL PGEND
      ENDIF 

      END


*************************************************************************
      SUBROUTINE MapSurf3 (Display, TMyr, Drel, Vrel)
*************************************************************************
*     modification of MapSurf2gif
*     visualizes surface density and particle plots

*     display = 3x2 pannels:
*                   P+S (surf. dens.)   S (surf. dens.)  S (particles)
*                   zoom                zoom             zoom      

*     falling particles are blue, expanding red


      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'merge.inp'

      REAL*4 TR(6), TRzoom(6), XMAXext, XMAXint, VMaxExt, VMaxInt
      REAL*4 ArrMin, ArrMax1, ArrMax2, ArrMax3
      REAL*4 DUMMY2D(ND,ND)

      CHARACTER GifFile*9, Display*3

      COMMON /GXIES2/ GMp, Bp, Bp2, YBp2, CTrBp, Sigma0p, Vesc0p, Sr0p,
     &                GMs, Bs, Bs2, YBs2, CTrBs, Sigma0s, Vesc0s, Sr0s,
     &                GMeff, GMStoP, GMpart1,GMpart2, Eps2,Eps2s,Eps2p, 
     &                Omega23, Dini, Vini, Trmin, DTrBp, DTrBs, Eps,
     &                RtrBp, RtrBs, BgasP, BgasS, BgasP2, BgasS2,
     &                YBgasP2, YBgasS2, CTrBgasP, CtrBgasS, Vesc0gasP,
     &                Vesc0gasS, Sr0gasP, Sr0gasS, DTrBgasP, DTrBgasS,
     &                RTrBgasP, RTrBgasS, Eps2gasP, Eps2gasS   
      COMMON /DENS2/  RHOpA (ND,ND), RHOpN (ND,ND), 
     &                RHOsA (ND,ND), RHOsN (ND,ND), 
     &                RHOtA (ND,ND), RHOtN (ND,ND), Ratio(ND,ND)
      COMMON /DENS2z/ RHOpAz(NDz,NDz), RHOpNz (NDz,NDz), ! zoom grid
     &                RHOsAz(NDz,NDz), RHOsNz (NDz,NDz), 
     &                RHOtAz(NDz,NDz), RHOtNz (NDz,NDz), Ratioz(NDz,NDz)
      SAVE   /DENS2/, /DENS2z/
                      ! RHOp  - surface density of the primary
                      ! RHOsA - (analytical) surface density of the secondary
                      ! RHOsN - surface density of the test particles
                      ! RHOtA = RHOp + RHOsA 
                      ! RHOtN = RHOp + RHOsN
      COMMON /TR6/    TR, TRzoom, XMAXext, XMAXint, VMaxExt, VMaxInt 
      COMMON /PANELS/ NCOLS, NROWS, NCRS
      COMMON /TestP1/ N1, N2, NALL    
      COMMON /DUMMY/  DUMMY2D

      SAVE   /TR6/, /PANELS/, /DUMMY/
      DATA NoFile3 /0/
      SAVE NoFile3    

      ArrMin = 0.d0;   ArrMax1 = SNGL(1.*Sigma0p);  
                     ! ArrMax2 = 5.e-3*ArrMax1
                       ArrMax2 = SNGL(0.002*Sigma0s)
                       ArrMax3 = 0.1

      IF (Display .EQ. 'gif') THEN
        CALL FILENAME3 (NoFile3, 'sh3-', 'gif', GifFile)
        CALL PGBEGIN(0,GifFile//'.gif/gif', NCOLS, NROWS)
        CALL ENV3x2pV
        CALL PALETT(3, 1.0, 0.5)  ! 3 = heat
      ELSEIF (Display .EQ. 'cps') THEN
        CALL FILENAME3 (NoFile3, 'sh3-', 'gif', GifFile)
        CALL PGBEGIN(0,GifFile//'.cps/cps', NCOLS, NROWS)
        CALL ENV3x2
        CALL PALETT(3, 1.0, 0.5)  ! 3 = heat
      ENDIF

      CALL PartPlot0 (N1, N2, Nfall, Nexp) ! separates particles to falling
                                           !    and outflowing
* upper pannels:

      CALL FindPanel;  CALL PGSWIN(-XMaxExt, XMaxExt, -XMaxExt, XMaxExt)  
      CALL PGIMAG(DUMMY2D, ND, ND, 1, ND, 1, ND, ArrMin, ArrMax1, TR)  
      CALL PartPlotXY (Nfall, Nexp, 4, 2) ! falling blue, expanding red
      CALL LabelMap (TMyr, Drel, Vrel)

      CALL FindPanel;  CALL PGSWIN(-XMaxExt, XMaxExt, -XMaxExt, XMaxExt)  
      CALL PGIMAG(DUMMY2D, ND, ND, 1, ND, 1, ND, ArrMin, ArrMax1, TR)  
      CALL PartPlotZY (Nfall, Nexp, 4, 2) ! falling blue, expanding red

      CALL FindPanel; CALL PGSWIN(-XMaxExt, XMaxExt, -VMaxExt, VMaxExt)  
      CALL PGIMAG(DUMMY2D, ND, ND, 1, ND, 1, ND, ArrMin, ArrMax1, TR)
      CALL PVPLOT2(N1, N2, 4,2)! compute & image position-velocity (r-v_r) plot

* lower pannels:
      CALL FindPanel;  CALL PGSWIN(-XMaxInt, XMaxInt, -XMaxInt, XMaxInt)  
      CALL PGIMAG(DUMMY2D, ND, ND, 1, ND, 1, ND, ArrMin, ArrMax1, TR)  
      CALL PartPlotXY (Nfall, Nexp, 4, 2) ! falling blue, expanding red

      CALL FindPanel;  CALL PGSWIN(-XMaxInt, XMaxInt, -XMaxInt, XMaxInt)  
      CALL PGIMAG(DUMMY2D, ND, ND, 1, ND, 1, ND, ArrMin, ArrMax1, TR)  
      CALL PartPlotZY (Nfall, Nexp, 4, 2) ! falling blue, expanding red

      CALL FindPanel; CALL PGSWIN(-XMaxInt, XMaxInt, -VMaxInt, VMaxInt)  
      CALL PGIMAG(DUMMY2D, ND, ND, 1, ND, 1, ND, ArrMin, ArrMax1, TR)
      CALL PVPLOT2(N1, N2, 4,2)! compute & image position-velocity (r-v_r) plot

      IF (Display .EQ. 'gif' .OR. Display .EQ. 'cps') THEN
        NoFile3 = NoFile3 + 1       
        CALL PGEND
      ENDIF 

      END 

*************************************************************************
      SUBROUTINE MapSurf4 (Display, TMyr, Drel, Vrel)
*************************************************************************
*     visualizes surface density for different inclinations
*     display = 3x2 pannels:

      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'merge.inp'

      REAL*4 TR(6), TRzoom(6), XMAXext, XMAXint, VMaxExt, VMaxInt
      REAL*4 ArrMin, ArrMax1, ArrMax2, ArrMax3
      CHARACTER GifFile*9, Display*3

      COMMON /GXIES2/ GMp, Bp, Bp2, YBp2, CTrBp, Sigma0p, Vesc0p, Sr0p,
     &                GMs, Bs, Bs2, YBs2, CTrBs, Sigma0s, Vesc0s, Sr0s,
     &                GMeff, GMStoP, GMpart1,GMpart2, Eps2,Eps2s,Eps2p, 
     &                Omega23, Dini, Vini, Trmin, DTrBp, DTrBs, Eps,
     &                RtrBp, RtrBs, BgasP, BgasS, BgasP2, BgasS2,
     &                YBgasP2, YBgasS2, CTrBgasP, CtrBgasS, Vesc0gasP,
     &                Vesc0gasS, Sr0gasP, Sr0gasS, DTrBgasP, DTrBgasS,
     &                RTrBgasP, RTrBgasS, Eps2gasP, Eps2gasS   
      COMMON /DENS2/  RHOpA (ND,ND), RHOpN (ND,ND), 
     &                RHOsA (ND,ND), RHOsN (ND,ND), 
     &                RHOtA (ND,ND), RHOtN (ND,ND), Ratio(ND,ND)
      COMMON /DENS2z/ RHOpAz(NDz,NDz), RHOpNz (NDz,NDz), ! zoom grid
     &                RHOsAz(NDz,NDz), RHOsNz (NDz,NDz), 
     &                RHOtAz(NDz,NDz), RHOtNz (NDz,NDz), Ratioz(NDz,NDz)
      SAVE   /DENS2/, /DENS2z/
                      ! RHOp  - surface density of the primary
                      ! RHOsA - (analytical) surface density of the secondary
                      ! RHOsN - surface density of the test particles
                      ! RHOtA = RHOp + RHOsA 
                      ! RHOtN = RHOp + RHOsN 
      COMMON /TR6/    TR, TRzoom, XMAXext, XMAXint, VMaxExt, VMaxInt 
      COMMON /PANELS/ NCOLS, NROWS, NCRS
      SAVE   /TR6/,/PANELS/

      DATA NoFile4 /0/
      SAVE NoFile4    

      ArrMin = 0.d0;   ArrMax1 = SNGL(1.*Sigma0p);  
                     ! ArrMax2 = 5.e-3*ArrMax1
                       ArrMax2 = SNGL(0.002*Sigma0s)
                       ArrMax3 = 0.1

      IF (Display .EQ. 'gif') THEN
        CALL FILENAME3 (NoFile4, 'sh4-', 'gif', GifFile)
        CALL PGBEGIN(0,GifFile//'.gif/gif', NCOLS, NROWS)
        CALL ENV3x2
        CALL PALETT(3, 1.0, 0.5)  ! 3 = heat
      ELSEIF (Display .EQ. 'cps') THEN
        CALL FILENAME3 (NoFile4, 'sh4-', 'gif', GifFile)
        CALL PGBEGIN(0,GifFile//'.cps/cps', NCOLS, NROWS)
        CALL ENV3x2
        CALL PALETT(3, 1.0, 0.5)  ! 3 = heat
      ENDIF

      DO 300 I = 1, 6       ! loop over inclinations

        AIdeg = (I-1) * 18.   ! inclination in degrees
        CALL CompSurfN (0.d0, AIdeg)   ! compute projected surf. density 

        CALL FindPanel;  CALL PGSWIN(-XMaxExt,XMaxExt, -XMaxExt,XMaxExt)
        CALL MAPRHO (RHOsN,  ND,  TR, ArrMin, ArrMax2, 1) 
        CALL LabelMapI (TMyr, AIdeg)

 300  CONTINUE  


      IF (Display .EQ. 'gif' .OR. Display .EQ. 'cps') THEN
        NoFile4 = NoFile4 + 1       
        CALL PGEND
      ENDIF 

      END

*************************************************************************
      SUBROUTINE MapSurf5 (Display, TMyr, Drel, Vrel)
*************************************************************************
*     visualizes surface density for different inclinations
*     display = 3x2 pannels:

      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'merge.inp'

      REAL*4 TR(6), TRzoom(6), XMAXext, XMAXint, VMaxExt, VMaxInt
      REAL*4 ArrMin, ArrMax1, ArrMax2, ArrMax3
      CHARACTER GifFile*9, Display*3

      COMMON /GXIES2/ GMp, Bp, Bp2, YBp2, CTrBp, Sigma0p, Vesc0p, Sr0p,
     &                GMs, Bs, Bs2, YBs2, CTrBs, Sigma0s, Vesc0s, Sr0s,
     &                GMeff, GMStoP, GMpart1,GMpart2, Eps2,Eps2s,Eps2p, 
     &                Omega23, Dini, Vini, Trmin, DTrBp, DTrBs, Eps,
     &                RtrBp, RtrBs, BgasP, BgasS, BgasP2, BgasS2,
     &                YBgasP2, YBgasS2, CTrBgasP, CtrBgasS, Vesc0gasP,
     &                Vesc0gasS, Sr0gasP, Sr0gasS, DTrBgasP, DTrBgasS,
     &                RTrBgasP, RTrBgasS, Eps2gasP, Eps2gasS   
      COMMON /DENS2/  RHOpA (ND,ND), RHOpN (ND,ND), 
     &                RHOsA (ND,ND), RHOsN (ND,ND), 
     &                RHOtA (ND,ND), RHOtN (ND,ND), Ratio(ND,ND)
      COMMON /DENS2z/ RHOpAz(NDz,NDz), RHOpNz (NDz,NDz), ! zoom grid
     &                RHOsAz(NDz,NDz), RHOsNz (NDz,NDz), 
     &                RHOtAz(NDz,NDz), RHOtNz (NDz,NDz), Ratioz(NDz,NDz)
      SAVE   /DENS2/, /DENS2z/
                      ! RHOp  - surface density of the primary
                      ! RHOsA - (analytical) surface density of the secondary
                      ! RHOsN - surface density of the test particles
                      ! RHOtA = RHOp + RHOsA 
                      ! RHOtN = RHOp + RHOsN 
      COMMON /TR6/    TR, TRzoom, XMAXext, XMAXint, VMaxExt, VMaxInt 
      COMMON /PANELS/ NCOLS, NROWS, NCRS
      SAVE   /TR6/,/PANELS/

      DATA NoFile5 /0/
      SAVE NoFile5    

      ArrMin = 0.d0;   ArrMax1 = SNGL(1.*Sigma0p);  
                     ! ArrMax2 = 5.e-3*ArrMax1
                       ArrMax2 = SNGL(0.002*Sigma0s)
                       ArrMax3 = 0.1

      IF (Display .EQ. 'gif') THEN
        CALL FILENAME3 (NoFile5, 'sh5-', 'gif', GifFile)
        CALL PGBEGIN(0,GifFile//'.gif/gif', NCOLS, NROWS)
        CALL ENV3x2zoom
        CALL PALETT(3, 1.0, 0.5)  ! 3 = heat
      ELSEIF (Display .EQ. 'cps') THEN
        CALL FILENAME3 (NoFile5, 'sh5-', 'gif', GifFile)
        CALL PGBEGIN(0,GifFile//'.cps/cps', NCOLS, NROWS)
        CALL ENV3x2zoom
        CALL PALETT(3, 1.0, 0.5)  ! 3 = heat
      ENDIF

      DO 300 I = 1, 6       ! loop over inclinations I>0

        AIdeg = (I-1) * 18.   ! inclination in degrees
        CALL CompSurfN (0.d0, AIdeg)   ! compute projected surf. density 


        CALL FindPanel;  CALL PGSWIN(-XMaxInt,XMaxInt, -XMaxInt,XMaxInt)  
        CALL MAPRHO (RHOsNz, NDz, TRzoom, ArrMin, ArrMax2, 1) 
        CALL LabelMapI (TMyr, AIdeg)

 300  CONTINUE  


      IF (Display .EQ. 'gif' .OR. Display .EQ. 'cps') THEN
        NoFile5 = NoFile5 + 1       
        CALL PGEND
      ENDIF 

      END


*************************************************************************
      SUBROUTINE MapSurf6 (Display, TMyr, Drel, Vrel)
*************************************************************************
*     visualizes surface density for different inclinations
*     like mapsurf4 but shows Ratio (unsharp masked image) 
*        instead of RhosN             

      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'merge.inp'

      REAL*4 TR(6), TRzoom(6), XMAXext, XMAXint, VMaxExt, VMaxInt
      REAL*4 ArrMin, ArrMax1, ArrMax2, ArrMax3
      CHARACTER GifFile*9, Display*3

      COMMON /GXIES2/ GMp, Bp, Bp2, YBp2, CTrBp, Sigma0p, Vesc0p, Sr0p,
     &                GMs, Bs, Bs2, YBs2, CTrBs, Sigma0s, Vesc0s, Sr0s,
     &                GMeff, GMStoP, GMpart1,GMpart2, Eps2,Eps2s,Eps2p, 
     &                Omega23, Dini, Vini, Trmin, DTrBp, DTrBs, Eps,
     &                RtrBp, RtrBs, BgasP, BgasS, BgasP2, BgasS2,
     &                YBgasP2, YBgasS2, CTrBgasP, CtrBgasS, Vesc0gasP,
     &                Vesc0gasS, Sr0gasP, Sr0gasS, DTrBgasP, DTrBgasS,
     &                RTrBgasP, RTrBgasS, Eps2gasP, Eps2gasS   
      COMMON /DENS2/  RHOpA (ND,ND), RHOpN (ND,ND), 
     &                RHOsA (ND,ND), RHOsN (ND,ND), 
     &                RHOtA (ND,ND), RHOtN (ND,ND), Ratio(ND,ND)
      COMMON /DENS2z/ RHOpAz(NDz,NDz), RHOpNz (NDz,NDz), ! zoom grid
     &                RHOsAz(NDz,NDz), RHOsNz (NDz,NDz), 
     &                RHOtAz(NDz,NDz), RHOtNz (NDz,NDz), Ratioz(NDz,NDz)
      SAVE   /DENS2/, /DENS2z/
                      ! RHOp  - surface density of the primary
                      ! RHOsA - (analytical) surface density of the secondary
                      ! RHOsN - surface density of the test particles
                      ! RHOtA = RHOp + RHOsA 
                      ! RHOtN = RHOp + RHOsN 
      COMMON /TR6/    TR, TRzoom, XMAXext, XMAXint, VMaxExt, VMaxInt 
      COMMON /PANELS/ NCOLS, NROWS, NCRS
      SAVE   /TR6/,/PANELS/

      DATA NoFile6 /0/
      SAVE NoFile6    

      ArrMin = 0.d0;   ArrMax1 = SNGL(1.*Sigma0p);  
                     ! ArrMax2 = 5.e-3*ArrMax1
                       ArrMax2 = SNGL(0.002*Sigma0s)
                       ArrMax3 = 0.1

      IF (Display .EQ. 'gif') THEN
        CALL FILENAME3 (NoFile6, 'sh6-', 'gif', GifFile)
        CALL PGBEGIN(0,GifFile//'.gif/gif', NCOLS, NROWS)
        CALL ENV3x2
        CALL PALETT(3, 1.0, 0.5)  ! 3 = heat
      ELSEIF (Display .EQ. 'cps') THEN
        CALL FILENAME3 (NoFile6, 'sh6-', 'gif', GifFile)
        CALL PGBEGIN(0,GifFile//'.cps/cps', NCOLS, NROWS)
        CALL ENV3x2
        CALL PALETT(3, 1.0, 0.5)  ! 3 = heat
      ENDIF

      DO 300 I = 1, 6       ! loop over inclinations

        AIdeg = (I-1) * 18.   ! inclination in degrees
        CALL CompSurfN (0.d0, AIdeg)   ! compute projected surf. density 

        CALL FindPanel;  CALL PGSWIN(-XMaxExt,XMaxExt, -XMaxExt,XMaxExt)
        CALL MAPRHO (Ratio,  ND,  TR, ArrMin, ArrMax3, 1) 
        CALL LabelMapI (TMyr, AIdeg)

 300  CONTINUE  


      IF (Display .EQ. 'gif' .OR. Display .EQ. 'cps') THEN
        NoFile6 = NoFile6 + 1       
        CALL PGEND
      ENDIF 

      END

*************************************************************************
      SUBROUTINE MapSurf7 (Display, TMyr, Drel, Vrel)
*************************************************************************
*     visualizes surface density for different inclinations
*     display = 3x2 pannels:
*     like mapsurf5 but shows Ratioz (zoomed unsharp masked image) 
*        instead of RhosNz             


      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'merge.inp'

      REAL*4 TR(6), TRzoom(6), XMAXext, XMAXint, VMaxExt, VMaxInt
      REAL*4 ArrMin, ArrMax1, ArrMax2, ArrMax3
      CHARACTER GifFile*9, Display*3

      COMMON /GXIES2/ GMp, Bp, Bp2, YBp2, CTrBp, Sigma0p, Vesc0p, Sr0p,
     &                GMs, Bs, Bs2, YBs2, CTrBs, Sigma0s, Vesc0s, Sr0s,
     &                GMeff, GMStoP, GMpart1,GMpart2, Eps2,Eps2s,Eps2p, 
     &                Omega23, Dini, Vini, Trmin, DTrBp, DTrBs, Eps,
     &                RtrBp, RtrBs, BgasP, BgasS, BgasP2, BgasS2,
     &                YBgasP2, YBgasS2, CTrBgasP, CtrBgasS, Vesc0gasP,
     &                Vesc0gasS, Sr0gasP, Sr0gasS, DTrBgasP, DTrBgasS,
     &                RTrBgasP, RTrBgasS, Eps2gasP, Eps2gasS   
      COMMON /DENS2/  RHOpA (ND,ND), RHOpN (ND,ND), 
     &                RHOsA (ND,ND), RHOsN (ND,ND), 
     &                RHOtA (ND,ND), RHOtN (ND,ND), Ratio(ND,ND)
      COMMON /DENS2z/ RHOpAz(NDz,NDz), RHOpNz (NDz,NDz), ! zoom grid
     &                RHOsAz(NDz,NDz), RHOsNz (NDz,NDz), 
     &                RHOtAz(NDz,NDz), RHOtNz (NDz,NDz), Ratioz(NDz,NDz)
      SAVE   /DENS2/, /DENS2z/
                      ! RHOp  - surface density of the primary
                      ! RHOsA - (analytical) surface density of the secondary
                      ! RHOsN - surface density of the test particles
                      ! RHOtA = RHOp + RHOsA 
                      ! RHOtN = RHOp + RHOsN 
      COMMON /TR6/    TR, TRzoom, XMAXext, XMAXint, VMaxExt, VMaxInt 
      COMMON /PANELS/ NCOLS, NROWS, NCRS
      SAVE   /TR6/,/PANELS/

      DATA NoFile7 /0/
      SAVE NoFile7    

      ArrMin = 0.d0;   ArrMax1 = SNGL(1.*Sigma0p);  
                     ! ArrMax2 = 5.e-3*ArrMax1
                       ArrMax2 = SNGL(0.002*Sigma0s)
                       ArrMax3 = 0.1

      IF (Display .EQ. 'gif') THEN
        CALL FILENAME3 (NoFile7, 'sh7-', 'gif', GifFile)
        CALL PGBEGIN(0,GifFile//'.gif/gif', NCOLS, NROWS)
        CALL ENV3x2zoom
        CALL PALETT(3, 1.0, 0.5)  ! 3 = heat
      ELSEIF (Display .EQ. 'cps') THEN
        CALL FILENAME3 (NoFile7, 'sh7-', 'gif', GifFile)
        CALL PGBEGIN(0,GifFile//'.cps/cps', NCOLS, NROWS)
        CALL ENV3x2zoom
        CALL PALETT(3, 1.0, 0.5)  ! 3 = heat
      ENDIF

      DO 300 I = 1, 6       ! loop over inclinations I>0

        AIdeg = (I-1) * 18.   ! inclination in degrees
        CALL CompSurfN (0.d0, AIdeg)   ! compute projected surf. density 


        CALL FindPanel;  CALL PGSWIN(-XMaxInt,XMaxInt, -XMaxInt,XMaxInt)  
        CALL MAPRHO (Ratioz, NDz, TRzoom, ArrMin, ArrMax3, 1) 
        CALL LabelMapI (TMyr, AIdeg)

 300  CONTINUE  


      IF (Display .EQ. 'gif' .OR. Display .EQ. 'cps') THEN
        NoFile7 = NoFile7 + 1       
        CALL PGEND
      ENDIF 

      END


*************************************************************************
      SUBROUTINE LabelMap0(P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11,P12,P13)
*************************************************************************
*     label the display at t=0 with input parameters

      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER C5*5, C6*20
      REAL*4 Vshift1, Vshift2
      PARAMETER (Vshift1 = +2.0, Vshift2 = +0.5 )

      CALL PGSLW (2); CALL PGSCH (2.)

      CALL PGNUMB(NINT(1.e-10*P1), 10,0,C6,NC);  CALL PGSCI(4)
      CALL PGMTXT ('T', Vshift1, 0.2, 0.5, 'M\dp\u='//C6)         ! Mp 

      CALL PGNUMB(NINT(1.e+3*P2), -3,0,C6,NC);  CALL PGSCI(4)
      CALL PGMTXT ('T', Vshift2, 0.2, 0.5, 'M\ds\u/M\dp\u='//C6)  ! Ms/Mp

      CALL PGNUMB(NINT(1.e+1*P3), -1,0,C5,NC);  CALL PGSCI(2)
      CALL PGMTXT ('T', Vshift1, 0.7, 0.5, 'a\dp\u='//C5)         ! a_p

      CALL PGNUMB(NINT(1.e+2*P4), -2,0,C5,NC);  CALL PGSCI(2)     ! a_s / a_p
      CALL PGMTXT ('T', Vshift2, 0.7, 0.5, 'a\ds\u/a\dp\u='//C5)  

      CALL PGNUMB(NINT(P7), 0,0,C5,NC);  CALL PGSCI(4)        
      CALL PGMTXT ('T', Vshift1, 1.1, 0.5, '\gS\dp\u='//C5)       ! Sigma_p

      CALL PGNUMB(NINT(P8), 0,0,C5,NC);  CALL PGSCI(4)
      CALL PGMTXT ('T', Vshift2, 1.1, 0.5, '\gS\ds\u='//C5)       ! Sigma_s

       CALL PGNUMB(NINT(P5), 0,0,C5,NC);  CALL PGSCI(2)            
      CALL PGMTXT ('T', Vshift1, 1.4, 0.5, '\gs\dp\u='//C5)       ! sigma_r,p

      CALL PGNUMB(NINT(P6), 0,0,C5,NC);  CALL PGSCI(2)
      CALL PGMTXT ('T', Vshift2, 1.4, 0.5, '\gs\ds\u='//C5)       ! sigma_r,s

      CALL PGNUMB(NINT(P9), 0,0,C5,NC);  CALL PGSCI(4)            
      CALL PGMTXT ('T', Vshift1, 1.8, 0.5, 'V\desc,p\u='//C5)     ! v_esc,p

      CALL PGNUMB(NINT(1.e+3*P10), -3,0,C6,NC);  CALL PGSCI(1)
      CALL PGMTXT ('T', Vshift1, 2.2, 0.5, '\gW\u2/3\d='//C6)

      CALL PGNUMB(NINT(1.e+3*P13), -3,0,C6,NC);  CALL PGSCI(2)            
      CALL PGMTXT ('T', Vshift2, 2.2, 0.5, 'T\dmin\u='//C6)   ! Trad,min [g.u.]

      CALL PGNUMB(NINT(P11), 0,0,C5,NC);  CALL PGSCI(4)
      CALL PGMTXT ('T', Vshift1, 2.6, 0.5, 'd\dini\u='//C5)

      CALL PGNUMB(NINT(P12), 0,0,C5,NC);  CALL PGSCI(2)
      CALL PGMTXT ('T', Vshift1, 3.0, 0.5, 'v\dini\u='//C5)

      CALL PGSLW (1); CALL PGSCH (1.);  CALL PGSCI(1)
 
      END

*************************************************************************
      SUBROUTINE LabelMap (Txt1, Txt2, Txt3)
*************************************************************************
*     labels 2D maps

      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER CTIME*5
       
      CALL PGSLW (2); CALL PGSCH (3.)

      CALL PGSCI(7) 
      CALL PGNUMB(NINT(Txt1), 0,0,CTIME,NSCAR)
      CALL PGMTXT ('T', -1.2, 0.2, 0.5, CTIME)

      CALL PGNUMB(NINT(Txt2), 0,0,CTIME,NSCAR)
      CALL PGMTXT ('T', -1.2, 0.6, 0.5, CTIME)

      CALL PGNUMB(NINT(Txt3), 0,0,CTIME,NSCAR)
      CALL PGMTXT ('T', -1.2, 0.8, 0.5, CTIME)

      CALL PGSLW (1); CALL PGSCH (1.); CALL PGSCI(1)
 
      END

*************************************************************************
      SUBROUTINE LabelMapI (Txt1, Txt2)
*************************************************************************
*     labels 2D maps

      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER CTIME*5
       
      CALL PGSLW (2); CALL PGSCH (3.)

      CALL PGSCI(7) 
      CALL PGNUMB(NINT(Txt1), 0,0,CTIME,NSCAR)
      CALL PGMTXT ('T', -1.2, 0.2, 0.5, CTIME)

      CALL PGNUMB(NINT(Txt2), 0,0,CTIME,NSCAR)
      CALL PGMTXT ('T', -1.2, 0.8, 0.5, CTIME)

      CALL PGSLW (1); CALL PGSCH (1.); CALL PGSCI(1) 
 
      END


*************************************************************************
      SUBROUTINE PARTPLOT2 (Nother, N, IColor1, IColor2)
*************************************************************************
*     plots particles

      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'merge.inp'

      PARAMETER (MaxPlot = 10000) ! max number of plotted particles
      REAL*4 Xfall(NX), Yfall(NX), Xexp(NX), Yexp(NX) 
      COMMON /TestP/  X(NX), Y(NX), Z(NX), VX(NX), VY(NX), VZ(NX)

      Dilution = FLOAT(MaxPlot)/N     ! dilution factor for particle plotting
                              
      Nfall = 0;  Nexp = 0    ! No. of particles at Vr<0 and Vr>0

      DO 10 I = Nother + 1, Nother + N

        IF (RAND(0) .GT. Dilution) GOTO 10   ! particle not plotted

        R = DSQRT (X(I)**2 + Y(I)**2 + Z(I)**2)  ! radial dist. from CM
        VR= ( X(I) * VX(I) + Y(I) * VY(I) + Z(I) * VZ(I) ) / R

        IF (VR .LE. 0.d0) THEN   ! falling particle, blue 
          Nfall = Nfall + 1
          Xfall (Nfall) = SNGL(X(I));    Yfall (Nfall) = SNGL(Y(I))
        ELSE
          Nexp = Nexp + 1
          Xexp  (Nexp) = SNGL(X(I));     Yexp  (Nexp) =  SNGL(Y(I))
        ENDIF


 10   CONTINUE 

!      print*, Nfall, Nexp

      CALL PGSCI (IColor1) 
      CALL PGPOINT(Nfall, Xfall, Yfall, -1)
      CALL PGSCI (IColor2) 
      CALL PGPOINT(Nexp,  Xexp,  Yexp,  -1)
      CALL PGSCI (1) 

      END

*************************************************************************
      SUBROUTINE PartPlotXY (Nfall, Nexp, IColor1, IColor2)
*************************************************************************
*     plots particles in XY plane

      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'merge.inp'

      REAL*4 Xfall(NX),Yfall(NX),Zfall(NX), Xexp(NX),Yexp(NX),Zexp(NX) 
      COMMON /TestP2/ Xfall, Yfall, Zfall, Xexp, Yexp, Zexp    
      SAVE /TestP2/

      CALL PGSCI (IColor1) 
      CALL PGPOINT(Nfall, Xfall, Yfall, -1)
      CALL PGSCI (IColor2) 
      CALL PGPOINT(Nexp,  Xexp,  Yexp,  -1)
      CALL PGSCI (1) 

      END


*************************************************************************
      SUBROUTINE PartPlotZY (Nfall, Nexp, IColor1, IColor2)
*************************************************************************
*     plots particles in ZY plane

      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'merge.inp'

      REAL*4 Xfall(NX),Yfall(NX),Zfall(NX), Xexp(NX),Yexp(NX),Zexp(NX) 
      COMMON /TestP2/ Xfall, Yfall, Zfall, Xexp, Yexp, Zexp    
      SAVE /TestP2/

      CALL PGSCI (IColor1) 
      CALL PGPOINT(Nfall, Zfall, Yfall, -1)
      CALL PGSCI (IColor2) 
      CALL PGPOINT(Nexp,  Zexp,  Yexp,  -1)
      CALL PGSCI (1) 

      END

*************************************************************************
      SUBROUTINE PartPlot0 (Nother, N, Nfall, Nexp)
*************************************************************************
*     separates particles to falling (blue) and outflowing (red)

      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'merge.inp'

      PARAMETER (MaxPlot = 10000) ! max number of plotted particles
      REAL*4 Xfall(NX),Yfall(NX),Zfall(NX), Xexp(NX),Yexp(NX),Zexp(NX) 
      COMMON /TestP/  X(NX), Y(NX), Z(NX), VX(NX), VY(NX), VZ(NX)
      COMMON /TestP2/ Xfall, Yfall, Zfall, Xexp, Yexp, Zexp    
      SAVE /TestP2/

      Dilution = FLOAT(MaxPlot)/N     ! dilution factor for particle plotting

      Nfall = 0;  Nexp = 0    ! No. of particles at Vr<0 and Vr>0

      DO 10 I = Nother + 1, Nother + N

        IF (RAND(0) .GT. Dilution) GOTO 10   ! particle not plotted

        R = DSQRT (X(I)**2 + Y(I)**2 + Z(I)**2)  ! radial dist. from CM
        VR= ( X(I) * VX(I) + Y(I) * VY(I) + Z(I) * VZ(I) ) / R

        IF (VR .LE. 0.d0) THEN   ! falling particle, blue 
          Nfall = Nfall + 1
          Xfall (Nfall) = SNGL(X(I));    Yfall (Nfall) = SNGL(Y(I))
                                         Zfall (Nfall) = SNGL(Z(I))
        ELSE
          Nexp = Nexp + 1
          Xexp  (Nexp) = SNGL(X(I));     Yexp  (Nexp) =  SNGL(Y(I))
                                         Zexp  (Nexp) =  SNGL(Z(I))
        ENDIF

 10   CONTINUE 

      END


*************************************************************************
      SUBROUTINE PVPLOT (Nother, N, IColor)
*************************************************************************
*     plots particles in the P-V diagram (radius vs. radial velocity)

      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'merge.inp'

      REAL*4 XS(NX), YS(NX) 
      COMMON /TestP/  X(NX), Y(NX), Z(NX), VX(NX), VY(NX), VZ(NX)
      COMMON /UNITS/ TUNIT,VUNIT,YVUNIT,GMUNIT,GMUNIT6


      DO 10 I = Nother + 1, Nother + N
        R = DSQRT (X(I)**2 + Y(I)**2 + Z(I)**2)  ! radial dist. from CM
        VR= ( X(I) * VX(I) + Y(I) * VY(I) + Z(I) * VZ(I) ) / R 
        XS (I-Nother) = SNGL(R);    
!       YS (I-Nother) = SNGL(VR*VUNIT/10) ! VR/10 plotted
        YS (I-Nother) = SNGL(VR*VUNIT   ) ! VR  plotted
                                    
 10   CONTINUE 

      CALL PGSCI (IColor) 
      CALL PGPOINT(N,XS,YS,-1)
      CALL PGSCI (1) 

      END

*************************************************************************
      SUBROUTINE PVPLOT1 (Nother, N, IColor1, IColor2)
*************************************************************************
*     plots particles in the P-V diagram (radius vs. radial velocity)
*     like PVPLOT but particles at X<0 and X>0 are plotted with diff. color 

      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'merge.inp'

      REAL*4 Rleft(NX), Rright(NX), Vleft(NX), Vright(NX)  
      COMMON /TestP/  X(NX), Y(NX), Z(NX), VX(NX), VY(NX), VZ(NX)
      COMMON /UNITS/ TUNIT,VUNIT,YVUNIT,GMUNIT,GMUNIT6

      Nleft = 0;  Nright = 0    ! No. of particles at X<0 and X>0

      DO 10 I = Nother + 1, Nother + N
        R = DSQRT (X(I)**2 + Y(I)**2 + Z(I)**2)  ! radial dist. from CM
        VR= ( X(I) * VX(I) + Y(I) * VY(I) + Z(I) * VZ(I) ) / R
        IF (X(I) .LE. 0.d0) THEN
          Nleft = Nleft + 1;
          Rleft (Nleft) = -SNGL(R)    
          Vleft (Nleft) = SNGL(VR*VUNIT/10)
        ELSE
          Nright = Nright + 1;
          Rright (Nright) = SNGL(R)    
          Vright (Nright) = SNGL(VR*VUNIT/10)
        ENDIF
                                    ! VR/10 plotted
 10   CONTINUE 

      CALL PGSCI (IColor1) 
      CALL PGPOINT(Nleft, Rleft, Vleft, -1)
      CALL PGSCI (IColor2) 
      CALL PGPOINT(Nright,Rright,Vright,-1)
      CALL PGSCI (1) 

      END

*************************************************************************
      SUBROUTINE PVPLOT2 (Nother, N, IColor1, IColor2)
*************************************************************************
*     plots particles in the P-V diagram (radius vs. radial velocity)
*     like PVPLOT1 but particles with Vr < 0 and Vr >0 are plotted 
*        with diff. color: 
*        Vr < 0: falling particles, blue;
*        Vr > 0: expanding particles, red
 

      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'merge.inp'

      REAL*4 Rfall(NX), Rexp(NX), Vfall(NX), Vexp(NX)  
      COMMON /TestP/  X(NX), Y(NX), Z(NX), VX(NX), VY(NX), VZ(NX)
      COMMON /UNITS/ TUNIT,VUNIT,YVUNIT,GMUNIT,GMUNIT6

      Nfall = 0;  Nexp = 0    ! No. of particles at Vr<0 and Vr>0

      DO 10 I = Nother + 1, Nother + N
        R = DSQRT (X(I)**2 + Y(I)**2 + Z(I)**2)  ! radial dist. from CM
        VR= ( X(I) * VX(I) + Y(I) * VY(I) + Z(I) * VZ(I) ) / R
        IF (VR .LE. 0.d0) THEN   ! falling particle, blue 
          Nfall = Nfall + 1
          IF (X(I) .LT. 0.d0) THEN
            Rfall (Nfall) = -SNGL(R)
          ELSE    
            Rfall (Nfall) =  SNGL(R)
          ENDIF
          Vfall (Nfall) = SNGL(VR*VUNIT/10) ! VR/10 plotted
!         Vfall (Nfall) = SNGL(VR*VUNIT   ) ! VR    plotted
        ELSE
          Nexp = Nexp + 1;
          IF (X(I) .LT. 0.d0) THEN
            Rexp (Nexp) = -SNGL(R)
          ELSE    
            Rexp (Nexp) =  SNGL(R)
          ENDIF
          Vexp (Nexp) = SNGL(VR*VUNIT/10)   ! VR/10 plotted
!         Vexp (Nexp) = SNGL(VR*VUNIT   )   ! VR    plotted
        ENDIF
                                    
 10   CONTINUE 

      CALL PGSCI (IColor1) 
      CALL PGPOINT(Nfall, Rfall, Vfall, -1)
      CALL PGSCI (IColor2) 
      CALL PGPOINT(Nexp,  Rexp,  Vexp,  -1)
      CALL PGSCI (1) 

      END

*************************************************************************
      SUBROUTINE FindPanel
*************************************************************************
*     find pgplot subpanel

      COMMON /PANELS/ NCOLS, NROWS, NCRS    ! NCRS = NCols * NRows
      SAVE /PANELS/, ICALL
      DATA ICALL /0/

      ICALL = ICALL + 1
      IPAN  = MOD(ICALL-1, NCRS) + 1 
      IF (IPAN .GT. NCOLS) THEN 
        IROW = 2; ICOL = IPAN-NCOLS
      ELSE
        IROW = 1; ICOL = IPAN
      ENDIF 

      CALL PGPANL   (ICOL, IROW)

      END

*************************************************************************
      SUBROUTINE MAPRHO (ARRAY, ND0, TR0, ARRMIN, ARRMAX, LSCALE) 
*************************************************************************
*     map 2D density arrays

      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'merge.inp'

      DIMENSION ARRAY(ND0,ND0)

      REAL*4 RHO4 (ND0,ND0), ARRMIN, ARRMAX
      REAL*4 TR0(6)

      IF (LSCALE .EQ. 1) THEN    ! linear scale 

        DO 10 IY=1, ND0
          DO 10 IX=1, ND0
            RHO4(IX,IY) = SNGL(ARRAY(IX,IY)) 
 10     CONTINUE
        
      ELSE                       ! logarithmic scale 
 
        DO 20 IY=1, ND0
          DO 20 IX=1, ND0
            IF (ARRAY(IX,IY) .GT. 1.d-16) THEN
              RHO4(IX,IY) = ALOG10(SNGL(ARRAY(IX,IY))) + 8.e0
            ELSE
              RHO4(IX,IY) = -8.e0
            ENDIF 
 20     CONTINUE

      ENDIF

      CALL PGIMAG(RHO4, ND0, ND0, 1, ND0, 1, ND0, ArrMin, ArrMax, TR0)

      END

************************************************************************
      SUBROUTINE INI1 (ProfileP, ProfileS, PhiP)  
************************************************************************

      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'merge.inp'

      PARAMETER (PI = 3.141592653589d0, PI2 = 2.d0*PI)
      CHARACTER ProfileP*1, ProfileS*1, PhiP*1
      COMMON /UNITS/ TUNIT,VUNIT,YVUNIT,GMUNIT,GMUNIT6
      COMMON /GXIES/  Xp(3),Xs(3),Vp(3),Vs(3)
      COMMON /GXIES2/ GMp, Bp, Bp2, YBp2, CTrBp, Sigma0p, Vesc0p, Sr0p,
     &                GMs, Bs, Bs2, YBs2, CTrBs, Sigma0s, Vesc0s, Sr0s,
     &                GMeff, GMStoP, GMpart1,GMpart2, Eps2,Eps2s,Eps2p, 
     &                Omega23, Dini, Vini, Trmin, DTrBp, DTrBs, Eps,
     &                RtrBp, RtrBs, BgasP, BgasS, BgasP2, BgasS2,
     &                YBgasP2, YBgasS2, CTrBgasP, CtrBgasS, Vesc0gasP,
     &                Vesc0gasS, Sr0gasP, Sr0gasS, DTrBgasP, DTrBgasS,
     &                RTrBgasP, RTrBgasS, Eps2gasP, Eps2gasS   
      COMMON /TestP1/ N1, N2, NALL
      COMMON /Gas1/   N1g, N2g, NALLg
      
           

! input parameters to be set in this routine:
!   PhiP- potential of the primary ('P'=Plummer, 'V'=de Vaucouleurs)
!   GMp - mass of the primary
!   Bp  - scale length of the primary
!   GMStoP - mass ratio secondary-to-primary
!   Bs  - scale length of the secondary
!   N1  - number of test particles for the primary
!   N2  - number of test particles for the secondary
!   ProfileP - density profile type for test particles of the primary
!   ProfileS - density profile type for test particles of the secondary
!   RtrBp    - truncation radius for test particles of the primary
!   RtrBs    - truncation radius for test particles of the primary
!   Dini     - initial separation of the galaxies
!   Vini       - initial relative velocity of the galaxies 

      

* primary galaxy
!     PhiP = 'V'   ! de Vaucouleurs potential
c      PhiP = 'P'   ! Plummer potential
      PhiP = 'H' ! Hernquist potential
      
      GMp  = 3.2d+11 / GMUNIT       ! mass
!     GMp  = 3.2d+12 / GMUNIT       ! mass

!     Bp   = 2.d0                  ! scale-length
      Bp   = 5.d0                  ! scale-length
!     Bp   = 28.d0                 ! Yoyo IR value

      Bp2  = Bp*Bp
      YBp2 = 1.d0 / Bp2 
      
      BgasP = 5.d0               ! scale-length for gas
      
      BgasP2 = BgasP*BgasP
      YBgasP2 = 1.d0 / BgasP2

      IF (PhiP .EQ. 'P') THEN   ! primary's potential = Plummer 
        Sigma0p = GMp * YBp2 / PI * GMUNIT6    ! ctrl surface density in Mo/pc2
        Vesc0p = DSQRT (GMp * 2.d0 / Bp)       ! ctrl escape velocity (Plummer = Hernquist)
        Vesc0gasP = DSQRT (GMp * 2.d0 / BgasP) ! same as above but for gas
        Sr0p   = DSQRT (GMp / 6.d0 / Bp)       ! ctrl velocity disp. 
        Sr0gasP = DSQRT (GMp / 6.d0 / BgasP)   ! same as above but for gas
        Trmin  = PI2 * DSQRT (Bp**3.d0 / GMp) * TUNIT * 1.e-3 ! Trad,min [Gyr]
      ELSEIF (PhiP .EQ. 'H') THEN ! primary's potential = Hernquist
        Sigma0p = GMp * YBp2 / (2.d0 * PI) * GMUNIT6 ! ctrl surface density in Mo/pc2
        Vesc0p = DSQRT (GMp * 2.d0 / Bp)       ! ctrl escape velocity (Plummer = Hernquist)
        Vesc0gasP = DSQRT (GMp * 2.d0 / BgasP) ! same as above but for gas
        Sr0p   = DSQRT (GMp / Bp)       ! ctrl velocity disp. 
        Sr0gasP = DSQRT (GMp / BgasP)   ! same as above but for gas
        Trmin  = PI2 * DSQRT (Bp**3.d0 / GMp) * TUNIT * 1.e-3 ! Trad,min [Gyr]
      ELSEIF (PhiP .EQ. 'V') THEN     ! primary's potential = de Vaucouleurs
        Q = 7.22    ! de Vauc. factor: 8!*exp(7.67)/7.67^8 
        Sigma0p = GMp * YBp2 / (Q*PI) * GMUNIT6  ! surf. dens. at Reff (deVauc)
        Vesc0p = DSQRT (GMp * 2.d0 / Bp)       ! ctrl escape velocity (Plummer = Hernquist)
        Sr0p   = DSQRT (GMp / 6.d0 / Bp)       ! ctrl velocity disp.  (Plummer)
        Trmin  = PI2 * DSQRT (Bp**3.d0 / GMp) * TUNIT * 1.e-3 ! Trad,min [Gyr] 
      ENDIF


* secondary-to-primary mass ratio
!     GMStoP  = 0.02d0
      GMStoP  = 0.010d0
      GMStoP1 = 1.d0 + GMStoP  
      YMStoP1 = 1.d0 / GMStoP1
      YMStoP2 = YMStoP1 * YMStoP1 

* secondary galaxy
      GMs  = GMstoP * GMp

      Bs   = 0.5d0
!     Bs   = 0.1d0 * Bp
!     Bs   = 0.027d0 * Bp    ! use if Bp=28 (Yoyo) to keep Bs fixed
!     Bs   = 0.76

      Bs2  = Bs*Bs
      YBs2 = 1.d0 / Bs2 
      
      BgasS = 0.5d0
      
      BgasS2 = BgasS * BgasS
      YBgasS2 = 1.d0 / BgasS2

      Sigma0s = GMs * YBs2 / PI * GMUNIT6
      Vesc0s = DSQRT (GMs * 2.d0 / Bs)  ! escape velocity at r=0
      Vesc0gasS = DSQRT (GMs * 2.d0 / BgasS)
      Sr0s   = DSQRT (GMs / 6.d0 / Bs)  ! rad. velocity disp. at r=0
      Sr0gasS = DSQRT (GMs / 6.d0 / BgasS)

      Omega23 = GMStoP * Bs/Bp             ! phase space volume^2/3 

      GMeff   = GMp * YMStoP2  

* softening for integrating motion of the secondary and of test particles
      Eps2s = Bs2     ! secondary's softening
      Eps2gasS = BgasS2
      IF (PhiP .EQ. 'P')  THEN   ! primary's potential = Plummer 
        Eps2p = Bp2   ! primary's softening = Plummer scale-length
        Eps2gasP = BgasP2
      ELSEIF (PhiP .EQ. 'V') THEN ! primary's potential = de Vaucouleurs
        Eps2p = DX2z  ! primary's softening = cell of the zoom mesh 
      ENDIF                                    
      Eps2    = (Eps2p  + Eps2s) * YMStoP2  ! total softening
      Eps = (Bp + Bs) * YMstoP1          ! total softening (for Herquist)                                               
*** test particles
      N1      = 0             ! number of test particles for the primary
      N2      = 1000000       ! number of test particles for the secondary
      NALL    = N1 + N2
      N1g     = 0
      N2g     = 0        ! number of gas particles for secondary
      NALLg   = N1g + N2g
      IF (NALL .GT. NX) THEN
!        print*,'Number of particles exceeds Nmax !!!'
        STOP 'Exit'
      ENDIF  

      GMpart1 = GMp * GMUNIT6 / N1 / DX2  ! particle mass per cell surface 
      GMpart2 = GMs * GMUNIT6 / N2 / DX2  ! particle mass per cell surface 
                       
                       ! type of the density profile:
                       !   P = Plummer, H = homog. sphere, S = spherical shell, Q = Hernquist 
      ProfileP = 'Q'   ! test particles of the primary 
      ProfileS = 'Q'   ! test particles of the secondary 

      IF     (ProfileP .EQ. 'P') THEN  
        RtrBp = 10.d0 * Bp            ! trunc. radius for Plummer
        CtrBp = (RtrBp**2.+Bp**2.)**1.5 / RtrBp**3.0  ! used for position gen
        DtrBp = 1.d0/(1.d0 + RtrBp**2.d0 * YBp2)**3.d0 ! used for velocity gen.
        RtrBgasP = 5.d0 * BgasP                                  ! values for gas
        CtrBgasP = (RtrBgasP**2.+BgasP**2.)**1.5 / RtrBgasP**3.0
        DtrBgasP = 1.d0/(1.d0 + RtrBgasP**2.d0 * YBgasP2)**3.d0
      ELSEIF (ProfileP .EQ. 'H') THEN
        RtrBp = 50.d0                 ! radius of the homogeneous sphere
        CtrBp = RtrBp
        DtrBp = 1.d0/(1.d0 + RtrBp**2.d0 * YBp2)**0.5d0 ! for vel. gen.
      ELSEIF (ProfileP .EQ. 'S') THEN 
        RtrBp = 50.d0                 ! radius of the shell
        CtrBp = RtrBp                 ! 
      ELSEIF (ProfileP .EQ. 'Q') THEN
        RtrBp = 10.d0 * Bp
        CtrBp = (RtrBp + Bp)**3.d0 / RtrBp**3.d0
        DtrBp = 1.d0 / (RtrBp + Bp)
        RtrBgasP = 5.d0 * BgasP
        CtrBgasP = (RtrBgasP + BgasP)**3.d0 / RtrBgasP**3.d0
        DtrBgasP = 1.d0 / (RtrBgasP + BgasP)
      ENDIF

      IF     (ProfileS .EQ. 'P') THEN  
        RtrBs = 5.d0 * Bs            ! trunc. radius for Plummer
        CtrBs = (RtrBs**2.+Bs**2.)**1.5 / RtrBs**3.0  ! used for Plummer
        DtrBs = 1.d0/(1.d0 + RtrBs**2.d0 * YBs2)**3.d0 ! used for velocity gen.
        RtrBgasS = 10.d0 * BgasS                                  ! values for gas
        CtrBgasS = (RtrBgasS**2.+BgasS**2.)**1.5 / RtrBgasS**3.0
        DtrBgasS = 1.d0/(1.d0 + RtrBgasS**2.d0 * YBgasS2)**3.d0
      ELSEIF (ProfileS .EQ. 'H') THEN
        RtrBs = 50.d0                 ! radius of the homogeneous sphere
        CtrBs = RtrBs                 ! 
        DtrBs = 1.d0/(1.d0 + RtrBs**2.d0 * YBs2)**0.5d0 ! for vel. gen.
      ELSEIF (ProfileS .EQ. 'S') THEN 
        RtrBs = 50.d0                 ! radius of the shell
        CtrBs = RtrBs                 ! 
      ELSEIF (ProfileS .EQ. 'Q') THEN
        RtrBs = 5.d0 * Bs
        CtrBs = (RtrBs + Bs)**3.d0 / RtrBs**3.d0
        DtrBs = 1.d0 / (RtrBs + Bs)
        RtrBgasS = 10.d0 * BgasS
        CtrBgasS = (RtrBgasS + BgasS)**3.d0 / RtrBgasS**3.d0
        DtrBgasS = 1.d0 / (RtrBgasS + BgasS)
      ENDIF
***


*  initial positions and velocities of galaxies
!     Dini = 12.*Bp   ! init. separ. as a multiple of primary's scale-length
      Dini = 91.2     ! 12 x 7.6 kpc
!      Dini = 20   

      Xp(2) = 0.d0;  Xp(3) = 0.d0;   ! both galaxies sitting on the x-axis
      Xs(2) = 0.d0;  Xs(3) = 0.d0;
                 
      Xs(1) = Dini * YMStoP1
      Xp(1) = - GMStoP * Xs(1)
                                     ! escape velocity of the secondary
      IF (ProfileP .EQ. 'Q') THEN
        Vesc = DSQRT (2.d0 * (GMp+GMs) / (Dini + Eps))
      ELSE
        Vesc  = DSQRT (2.d0 * (GMp+GMs) / DSQRT(Dini**2 + Eps2))
      ENDIF
      
      Torbit= 1.d0 ! orbit type (1=parabola, <1 bound orbit, >1 unbound orbit)
!     Torbit= 0.0d0 ! orbit type (1=parabola, <1 bound orbit, >1 unbound orbit)
      Vini    = Torbit * Vesc          ! initial relative velocity

      Vs2to1 = 0.0d0   ! Vy_secondary / Vx_secondary (if 0, then radial orbit) 
      YV     = 1.d0 / DSQRT(1.d0 + Vs2to1**2)

      Vs(1) =-Vini * YMStoP1 * YV;  Vs(2)= Vs(1) * Vs2to1;   Vs(3)= 0.d0
      Vp(1) =-GMStoP*Vs(1);  Vp(2) =-GMStoP*Vs(2);  Vp(3) =-GMStoP*Vs(3)     


c      WRITE(*,1) Sigma0p, Vesc0p*VUNIT, Sr0p*VUNIT, 
c     &           Sigma0s, Vesc0s*VUNIT, Sr0s*VUNIT
c      WRITE(*,2) Dini, Vesc*VUNIT, Vini*VUNIT
c      WRITE(*,3) Omega23
c
c 1    FORMAT (// 
c     &        '   Sigma (0)    Vesc (0)    sigma_r (0)'/
c     &        '   [Mo/pc^2]     [km/s]     [km/s]  '/
c     &        'P: ', F7.2, 6x,   F7.2, 5x,  F7.2      /        
c     &        'S: ', F7.2, 6x,   F7.2, 5x,  F7.2    / )  
c 2    FORMAT ('Initial separation, Esc. velocity, Rel. velocity   ',          
c     &         F7.2, 3x, F7.2, 3x, F7.2             / ) 
c 3    FORMAT ('Phase space volume of the secondary: ' /
c     &        'Omega^2/3 = Ms/Mp * Bs/Bp = ', F8.4 // )

      END


************************************************************************
      SUBROUTINE InitTPs (ProfileP, ProfileS)  ! initialize test particles
************************************************************************

      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'merge.inp'

      CHARACTER ProfileP*1, ProfileS*1
      COMMON /GXIES/  Xp(3),Xs(3),Vp(3),Vs(3)
      COMMON /GXIES2/ GMp, Bp, Bp2, YBp2, CTrBp, Sigma0p, Vesc0p, Sr0p,
     &                GMs, Bs, Bs2, YBs2, CTrBs, Sigma0s, Vesc0s, Sr0s,
     &                GMeff, GMStoP, GMpart1,GMpart2, Eps2,Eps2s,Eps2p, 
     &                Omega23, Dini, Vini, Trmin, DTrBp, DTrBs, Eps,
     &                RtrBp, RtrBs, BgasP, BgasS, BgasP2, BgasS2,
     &                YBgasP2, YBgasS2, CTrBgasP, CtrBgasS, Vesc0gasP,
     &                Vesc0gasS, Sr0gasP, Sr0gasS, DTrBgasP, DTrBgasS,
     &                RTrBgasP, RTrBgasS, Eps2gasP, Eps2gasS  
      COMMON /TestP1/ N1, N2, NALL 
      COMMON /Gas1/   N1g, N2g, NALLg    

* primary galaxy
      IF (ProfileP .EQ. 'P') THEN          ! Plummer  
        CALL GenPlummer    (0, N1, Bp, CTrBp, DTrBp, YBp2, Vesc0p, Sr0p)
        CALL GenPlummerGas (0, N1g, BgasP, CTrBgasP, DTrBgasP, YBgasP2,
     &                       Vesc0gasP, Sr0gasP)
      ELSEIF (ProfileP .EQ. 'Q') THEN !Hernquist
        CALL GenHernquist(0, N1, Bp, CTrBp, DTrBp, Vesc0p, RtrBp, GMp)
        CALL GenHernquistGas(0,N1g, BgasP, CTrBgasP, DTrBgasP,Vesc0gasP,
     &                        RtrBgasP, GMp)
      ELSEIF (ProfileP .EQ. 'H') THEN      ! Homog. sphere, radius = CTrBp   
        CALL GenHomSph  (0, N1,     CTrBp, DTrBp, YBp2, Vesc0p, Sr0p) 
      ELSE                                 ! Uniform shell, radius = CTrBp
        CALL GenShell   (0, N1,     CTrBp,        YBp2, Vesc0p      ) 
      ENDIF 

      IF (ProfileS .EQ. 'P') THEN          ! Plummer  
        CALL GenPlummer (N1, N2, Bs, CTrBs, DTrBs, YBs2, Vesc0s, Sr0s)
        CALL GenPlummerGas (N1g, N2g, BgasS, CTrBgasS, DTrBgasS,YBgasS2,
     &                       Vesc0gasS, Sr0gasS)
      ELSEIF (ProfileS .EQ. 'Q') THEN !Hernquist
        CALL GenHernquist(N1, N2, Bs, CTrBs, DTrBs, Vesc0s, RtrBs, GMs)
        CALL GenHernquistGas(N1g,N2g,BgasS,CTrBgasS, DTrBgasS,Vesc0gasS,
     &                        RtrBgasS, GMp)
      ELSEIF (ProfileS .EQ. 'H') THEN      ! Homog. sphere, radius = CTrBp   
        CALL GenHomSph  (N1, N2,     CTrBs, DTrBs, YBs2, Vesc0s, Sr0s) 
      ELSE                                 ! Uniform shell, radius = CTrBp
        CALL GenShell   (N1, N2,     CTrBs,        YBs2, Vesc0s      ) 
      ENDIF 


* transform positions and velocities to the CM frame:
      CALL CMFRAME  (0,  N1, Xp, Vp)
      CALL CMFRAME  (N1, N2, Xs, Vs)
      CALL CMFRAMEgas  (0,  N1g, Xp, Vp)
      CALL CMFRAMEgas  (N1g, N2g, Xs, Vs)

      END


************************************************************************
      SUBROUTINE CMFRAME (Nother, NPlummer, Xgal, Vgal)
************************************************************************
! transform particle positions and velocities from the galaxy frame 
!   to the CM frame

      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'merge.inp'

      DIMENSION Xgal(3), Vgal(3)
      COMMON /TestP/  XT(NX), YT(NX), ZT(NX), VXT(NX), VYT(NX), VZT(NX)

      DO 20 I = Nother + 1, Nother + NPlummer

        XT (I) = XT(I) + Xgal(1)       
        YT (I) = YT(I) + Xgal(2)       
        ZT (I) = ZT(I) + Xgal(3)       

        VXT (I) = VXT(I) + Vgal(1)       
        VYT (I) = VYT(I) + Vgal(2)       
        VZT (I) = VZT(I) + Vgal(3)       

 20   CONTINUE  

      END

************************************************************************
      SUBROUTINE CMFRAMEgas (Nother, NPlummer, Xgal, Vgal)
************************************************************************
! transform gas-particle positions and velocities from the galaxy frame 
!   to the CM frame

      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'merge.inp'

      DIMENSION Xgal(3), Vgal(3)
      COMMON /Gas/  Xg(NX), Yg(NX), Zg(NX), VXg(NX), VYg(NX), VZg(NX),
     &              Xhalf(NX), Yhalf(NX), Zhalf(NX), Mesh(NX) 

      DO 20 I = Nother + 1, Nother + NPlummer

        Xg (I) = Xg(I) + Xgal(1)       
        Yg (I) = Yg(I) + Xgal(2)       
        Zg (I) = Zg(I) + Xgal(3)       

        VXg (I) = VXg(I) + Vgal(1)       
        VYg (I) = VYg(I) + Vgal(2)       
        VZg (I) = VZg(I) + Vgal(3)       

 20   CONTINUE  

      END

************************************************************************
      SUBROUTINE GenPlummer (Nother, N, B, C, D, YB2, Vesc0, Sr0)
************************************************************************
* generates particle positions within Plummer sphere

      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'merge.inp'

      PARAMETER (PI=3.141592653589d0, PI2 = 2.d0*PI)
      REAL*4 RAND    
      COMMON /UNITS/ TUNIT,VUNIT,YVUNIT,GMUNIT,GMUNIT6      
      COMMON /TestP/  X(NX), Y(NX), Z(NX), VX(NX), VY(NX), VZ(NX)
      COMMON /GXIES2/ GMp, Bp, Bp2, YBp2, CTrBp, Sigma0p, Vesc0p, Sr0p,
     &                GMs, Bs, Bs2, YBs2, CTrBs, Sigma0s, Vesc0s, Sr0s,
     &                GMeff, GMStoP, GMpart1,GMpart2, Eps2,Eps2s,Eps2p, 
     &                Omega23, Dini, Vini, Trmin, DTrBp, DTrBs, Eps,
     &                RtrBp, RtrBs, BgasP, BgasS, BgasP2, BgasS2,
     &                YBgasP2, YBgasS2, CTrBgasP, CtrBgasS, Vesc0gasP,
     &                Vesc0gasS, Sr0gasP, Sr0gasS, DTrBgasP, DTrBgasS,
     &                RTrBgasP, RTrBgasS, Eps2gasP, Eps2gasS  

c     CALL SRAND(ISEED1)

      IF (N .EQ. 0) RETURN    

      NVhigh = 0  ! No. of repeatedly drawn particles (V > Vesc elimination)

      XCM  = 0.d0;   YCM  = 0.d0;   ZCM  = 0.d0    ! center-of-mass 
      VXCM = 0.d0;   VYCM = 0.d0;   VZCM = 0.d0    ! center-of-mass velocity 


* loop over particles:
      DO 2 I = Nother + 1, Nother + N, 2 

* generate angles
        PHI      = PI2   * DBLE(RAND(0))  
        COSTHETA = 2.d0  * DBLE(RAND(0)) - 1.d0
        SINTHETA = DSQRT(1.-COSTHETA*COSTHETA)

* generate radius       
 1913   R = GenPosPL0 (DBLE(Rand(0)), B, C)    
!       IF (R. LT. 0.001) GOTO 1913

* generate 3D position
        X(I) = R * DCOS(PHI) * SINTHETA
        Y(I) = R * DSIN(PHI) * SINTHETA
        Z(I) = R * COSTHETA

* velocity dispersion Sr
        Q = 1.d0 + R * R * YB2    
!       Sr    = Sr0   / Q**0.25d0      ! truncation not taken into account
        Sr    = Sr0   * DSQRT ( Q**2.5d0 * (Q**(-3.d0) - D) )

* escape velocity     
        Vesc  = Vesc0 / Q**0.25d0
        Vhigh = 0.99  * Vesc     

* generate velocities:
 95     VX(I) = Sr * GASDEV()          
        VY(I) = Sr * GASDEV()          
        VZ(I) = Sr * GASDEV()

        Vtot = DSQRT (VX(I)**2. + VY(I)**2. + VZ(I)**2.)

* check if particle not too fast
        IF (Vtot .GT. Vhigh) THEN 
          NVhigh = NVhigh + 1 
          GOTO 95            ! particle too fast --> draw velocities again
!          print*, 'Elimination: NVhigh = ', NVhigh 
        ENDIF

* pair particles by symmetry
        X (I+1) = -X (I);  Y (I+1) = -Y(I);   Z (I+1) = -Z(I)  
        VX(I+1) = -VX(I);  VY(I+1) = -VY(I);  VZ(I+1) = -VZ(I) 

* update center-of-mass
        XCM = XCM + X(I);    YCM = YCM + Y(I);    ZCM = ZCM + Z(I)
        XCM = XCM + X(I+1);  YCM = YCM + Y(I+1);  ZCM = ZCM + Z(I+1)

* update center-of-mass velocity 
        VXCM = VXCM +VX(I  ); VYCM = VYCM +VY(I  ); VZCM = VZCM +VZ(I  )
        VXCM = VXCM +VX(I+1); VYCM = VYCM +VY(I+1); VZCM = VZCM +VZ(I+1)
 2    CONTINUE
  
!      PRINT*,'No. of close-to-escape-velocity eliminations: ', NVhigh  



* center-of-mass:
      XCM  = XCM/N;    YCM  = YCM/N;    ZCM = ZCM/N
!      PRINT*,'CENTER-OF-MASS [kpc]: ', XCM, YCM, ZCM
      VXCM = VXCM/N;   VYCM = VYCM/N;   VZCM = VZCM/N
!      PRINT*,'CENTER-OF-MASS [km/s]: ',VXCM*VUNIT,VYCM*VUNIT,VZCM*VUNIT

* correct particle positions to the center-of-mass frame:
      DO 3, I = Nother + 1, Nother + N
        X (I) = X (I)-XCM;   Y(I)  = Y (I)-YCM;   Z(I)  = Z (I)-ZCM
        VX(I) = VX(I)-VXCM;  VY(I) = VY(I)-VYCM;  VZ(I) = VZ(I)-VZCM
 3    CONTINUE

      END
      
************************************************************************
      SUBROUTINE GenPlummerGas (Nother, N, B, C, D, YB2, Vesc0, Sr0)
************************************************************************
* generates gas-particle positions within Plummer sphere

      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'merge.inp'

      PARAMETER (PI=3.141592653589d0, PI2 = 2.d0*PI)
      REAL*4 RAND    
      COMMON /UNITS/ TUNIT,VUNIT,YVUNIT,GMUNIT,GMUNIT6      
      COMMON /Gas/   Xg(NX), Yg(NX), Zg(NX), VXg(NX), VYg(NX), VZg(NX),
     &               Xhalf(NX), Yhalf(NX), Zhalf(NX), Mesh(NX)
      COMMON /GXIES2/ GMp, Bp, Bp2, YBp2, CTrBp, Sigma0p, Vesc0p, Sr0p,
     &                GMs, Bs, Bs2, YBs2, CTrBs, Sigma0s, Vesc0s, Sr0s,
     &                GMeff, GMStoP, GMpart1,GMpart2, Eps2,Eps2s,Eps2p, 
     &                Omega23, Dini, Vini, Trmin, DTrBp, DTrBs, Eps,
     &                RtrBp, RtrBs, BgasP, BgasS, BgasP2, BgasS2,
     &                YBgasP2, YBgasS2, CTrBgasP, CtrBgasS, Vesc0gasP,
     &                Vesc0gasS, Sr0gasP, Sr0gasS, DTrBgasP, DTrBgasS,
     &                RTrBgasP, RTrBgasS, Eps2gasP, Eps2gasS  

c     CALL SRAND(ISEED1)

      IF (N .EQ. 0) RETURN    

      NVhigh = 0  ! No. of repeatedly drawn particles (V > Vesc elimination)

      XCM  = 0.d0;   YCM  = 0.d0;   ZCM  = 0.d0    ! center-of-mass 
      VXCM = 0.d0;   VYCM = 0.d0;   VZCM = 0.d0    ! center-of-mass velocity 


* loop over particles:
      DO 2 I = Nother + 1, Nother + N, 2 

* generate angles
        PHI      = PI2   * DBLE(RAND(0))  
        COSTHETA = 2.d0  * DBLE(RAND(0)) - 1.d0
        SINTHETA = DSQRT(1.-COSTHETA*COSTHETA)

* generate radius       
 1913   R = GenPosPL0 (DBLE(Rand(0)), B, C)    
c       IF (R. GT. B) GOTO 1913

* generate 3D position
        Xg(I) = R * DCOS(PHI) * SINTHETA
        Yg(I) = R * DSIN(PHI) * SINTHETA
        Zg(I) = R * COSTHETA

* velocity dispersion Sr
        Q = 1.d0 + R * R * YB2    
!       Sr    = Sr0   / Q**0.25d0      ! truncation not taken into account
        Sr    = Sr0   * DSQRT ( Q**2.5d0 * (Q**(-3.d0) - D) )

* escape velocity     
        Vesc  = Vesc0 / Q**0.25d0
        Vhigh = 0.99  * Vesc     

* generate velocities:
 95     VXg(I) = Sr * GASDEV()          
        VYg(I) = Sr * GASDEV()          
        VZg(I) = Sr * GASDEV()

        Vtot = DSQRT (VXg(I)**2. + VYg(I)**2. + VZg(I)**2.)

* check if particle not too fast
        IF (Vtot .GT. Vhigh) THEN 
          NVhigh = NVhigh + 1 
          GOTO 95            ! particle too fast --> draw velocities again
!          print*, 'Elimination: NVhigh = ', NVhigh 
        ENDIF

* pair particles by symmetry
        Xg (I+1) = -Xg (I);  Yg (I+1) = -Yg(I);   Zg (I+1) = -Zg(I)  
        VXg(I+1) = -VXg(I);  VYg(I+1) = -VYg(I);  VZg(I+1) = -VZg(I) 

* update center-of-mass
        XCM = XCM + Xg(I);    YCM = YCM + Yg(I);    ZCM = ZCM + Zg(I)
        XCM = XCM + Xg(I+1);  YCM = YCM + Yg(I+1);  ZCM = ZCM + Zg(I+1)

* update center-of-mass velocity 
        VXCM=VXCM+VXg(I  ); VYCM=VYCM+VYg(I  ); VZCM=VZCM +VZg(I  )
        VXCM=VXCM+VXg(I+1); VYCM=VYCM+VYg(I+1); VZCM=VZCM +VZg(I+1)
 2    CONTINUE
  
!      PRINT*,'No. of close-to-escape-velocity eliminations: ', NVhigh  



* center-of-mass:
      XCM  = XCM/N;    YCM  = YCM/N;    ZCM = ZCM/N
!      PRINT*,'CENTER-OF-MASS [kpc]: ', XCM, YCM, ZCM
      VXCM = VXCM/N;   VYCM = VYCM/N;   VZCM = VZCM/N
!      PRINT*,'CENTER-OF-MASS [km/s]: ',VXCM*VUNIT,VYCM*VUNIT,VZCM*VUNIT

* correct particle positions to the center-of-mass frame:
      DO 3, I = Nother + 1, Nother + N
        Xg (I) = Xg (I)-XCM;   Yg(I)  = Yg (I)-YCM;  Zg(I) = Zg (I)-ZCM
        VXg(I) = VXg(I)-VXCM;  VYg(I) = VYg(I)-VYCM; VZg(I)= VZg(I)-VZCM
 3    CONTINUE

      END
      
************************************************************************
      SUBROUTINE GenHernquist (Nother, N, B, C, D, Vesc0, Rtr, GM)
************************************************************************
* generates particle positions within Hernquist sphere

      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'merge.inp'

      PARAMETER (PI=3.141592653589d0, PI2 = 2.d0*PI)
      REAL*4 RAND    
      COMMON /UNITS/ TUNIT,VUNIT,YVUNIT,GMUNIT,GMUNIT6      
      COMMON /TestP/  X(NX), Y(NX), Z(NX), VX(NX), VY(NX), VZ(NX)
      COMMON /GXIES2/ GMp, Bp, Bp2, YBp2, CTrBp, Sigma0p, Vesc0p, Sr0p,
     &                GMs, Bs, Bs2, YBs2, CTrBs, Sigma0s, Vesc0s, Sr0s,
     &                GMeff, GMStoP, GMpart1,GMpart2, Eps2,Eps2s,Eps2p, 
     &                Omega23, Dini, Vini, Trmin, DTrBp, DTrBs, Eps,
     &                RtrBp, RtrBs, BgasP, BgasS, BgasP2, BgasS2,
     &                YBgasP2, YBgasS2, CTrBgasP, CtrBgasS, Vesc0gasP,
     &                Vesc0gasS, Sr0gasP, Sr0gasS, DTrBgasP, DTrBgasS,
     &                RTrBgasP, RTrBgasS, Eps2gasP, Eps2gasS
     
c     CALL SRAND(ISEED1)

      IF (N .EQ. 0) RETURN    

      NVhigh = 0  ! No. of repeatedly drawn particles (V > Vesc elimination)

      XCM  = 0.d0;   YCM  = 0.d0;   ZCM  = 0.d0    ! center-of-mass 
      VXCM = 0.d0;   VYCM = 0.d0;   VZCM = 0.d0    ! center-of-mass velocity 


* loop over particles:
      DO 2 I = Nother + 1, Nother + N, 2 

* generate angles
        PHI      = PI2   * DBLE(RAND(0))  
        COSTHETA = 2.d0  * DBLE(RAND(0)) - 1.d0
        SINTHETA = DSQRT(1.-COSTHETA*COSTHETA)

* generate radius       
 1913   R = GenPosHQ0(DBLE(Rand(0)), B, C)    
!       IF (R. LT. 0.001) GOTO 1913

* generate 3D position
        X(I) = R * DCOS(PHI) * SINTHETA
        Y(I) = R * DSIN(PHI) * SINTHETA
        Z(I) = R * COSTHETA

* velocity dispersion Sr
        Sr6 = B**(-4.) * DLOG(Rtr/R)
        Sr7 = B**(-4.) * DLOG((Rtr+B)/(R+B))
        Sr1 = DSQRT(GM/B * R*(R+B)**3.)
        Sr2 = B**(-3.) * (D - 1.d0/(R+B))
        Sr3 = 1.d0/(2*B**2.) * (D**2. - 1.d0/(R+B)**2.)
        Sr4 = 1.d0/(3*B) * (D**3. - 1.d0/(R+B)**3.)
        Sr5 = 1.d0/4.d0 * (D**4. - 1.d0/(R+B)**4.)
        Sr = Sr1 * DSQRT(Sr6 - Sr7 + Sr2 + Sr3 + Sr4 + Sr5)
        
* escape velocity
        Q = DSQRT(1.0 / (1. + R/B))
        Vesc  = Vesc0 / Q
        Vhigh = 0.99  * Vesc     

* generate velocities:
 95     VX(I) = Sr * GASDEV()          
        VY(I) = Sr * GASDEV()          
        VZ(I) = Sr * GASDEV()

        Vtot = DSQRT (VX(I)**2. + VY(I)**2. + VZ(I)**2.)

* check if particle not too fast
        IF (Vtot .GT. Vhigh) THEN 
          NVhigh = NVhigh + 1 
          GOTO 95            ! particle too fast --> draw velocities again
!          print*, 'Elimination: NVhigh = ', NVhigh 
        ENDIF

* pair particles by symmetry
        X (I+1) = -X (I);  Y (I+1) = -Y(I);   Z (I+1) = -Z(I)  
        VX(I+1) = -VX(I);  VY(I+1) = -VY(I);  VZ(I+1) = -VZ(I) 

* update center-of-mass
        XCM = XCM + X(I);    YCM = YCM + Y(I);    ZCM = ZCM + Z(I)
        XCM = XCM + X(I+1);  YCM = YCM + Y(I+1);  ZCM = ZCM + Z(I+1)

* update center-of-mass velocity 
        VXCM=VXCM+VX(I  ); VYCM=VYCM+VY(I  ); VZCM=VZCM +VZ(I )
        VXCM=VXCM+VX(I+1); VYCM=VYCM+VY(I+1); VZCM=VZCM +VZ(I+1)
 2    CONTINUE
  
!      PRINT*,'No. of close-to-escape-velocity eliminations: ', NVhigh  



* center-of-mass:
      XCM  = XCM/N;    YCM  = YCM/N;    ZCM = ZCM/N
!      PRINT*,'CENTER-OF-MASS [kpc]: ', XCM, YCM, ZCM
      VXCM = VXCM/N;   VYCM = VYCM/N;   VZCM = VZCM/N
!      PRINT*,'CENTER-OF-MASS [km/s]: ',VXCM*VUNIT,VYCM*VUNIT,VZCM*VUNIT

* correct particle positions to the center-of-mass frame:
      DO 3, I = Nother + 1, Nother + N
        X (I) = X (I)-XCM;   Y(I)  = Y (I)-YCM;  Z(I) = Z (I)-ZCM
        VX(I) = VX(I)-VXCM;  VY(I) = VY(I)-VYCM; VZ(I)= VZ(I)-VZCM
 3    CONTINUE

      END

************************************************************************
      SUBROUTINE GenHernquistGas (Nother, N, B, C, D, Vesc0, Rtr, GM)
************************************************************************
* generates particle positions within Hernquist sphere

      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'merge.inp'

      PARAMETER (PI=3.141592653589d0, PI2 = 2.d0*PI)
      REAL*4 RAND    
      COMMON /UNITS/ TUNIT,VUNIT,YVUNIT,GMUNIT,GMUNIT6      
      COMMON /Gas/   Xg(NX), Yg(NX), Zg(NX), VXg(NX), VYg(NX), VZg(NX),
     &               Xhalf(NX), Yhalf(NX), Zhalf(NX), Mesh(NX)
      COMMON /GXIES2/ GMp, Bp, Bp2, YBp2, CTrBp, Sigma0p, Vesc0p, Sr0p,
     &                GMs, Bs, Bs2, YBs2, CTrBs, Sigma0s, Vesc0s, Sr0s,
     &                GMeff, GMStoP, GMpart1,GMpart2, Eps2,Eps2s,Eps2p, 
     &                Omega23, Dini, Vini, Trmin, DTrBp, DTrBs, Eps,
     &                RtrBp, RtrBs, BgasP, BgasS, BgasP2, BgasS2,
     &                YBgasP2, YBgasS2, CTrBgasP, CtrBgasS, Vesc0gasP,
     &                Vesc0gasS, Sr0gasP, Sr0gasS, DTrBgasP, DTrBgasS,
     &                RTrBgasP, RTrBgasS, Eps2gasP, Eps2gasS

c     CALL SRAND(ISEED1)

      IF (N .EQ. 0) RETURN    

      NVhigh = 0  ! No. of repeatedly drawn particles (V > Vesc elimination)

      XCM  = 0.d0;   YCM  = 0.d0;   ZCM  = 0.d0    ! center-of-mass 
      VXCM = 0.d0;   VYCM = 0.d0;   VZCM = 0.d0    ! center-of-mass velocity 


* loop over particles:
      DO 2 I = Nother + 1, Nother + N, 2 

* generate angles
        PHI      = PI2   * DBLE(RAND(0))  
        COSTHETA = 2.d0  * DBLE(RAND(0)) - 1.d0
        SINTHETA = DSQRT(1.-COSTHETA*COSTHETA)

* generate radius       
 1913   R = GenPosHQ0(DBLE(Rand(0)), B, C)    
c       IF (R. GT. B) GOTO 1913

* generate 3D position
        Xg(I) = R * DCOS(PHI) * SINTHETA
        Yg(I) = R * DSIN(PHI) * SINTHETA
        Zg(I) = R * COSTHETA

* velocity dispersion Sr
        Sr6 = B**(-4.) * DLOG(Rtr/R)
        Sr7 = B**(-4.) * DLOG((Rtr+B)/(R+B))
        Sr1 = DSQRT(GM/B * R*(R+B)**3.)
        Sr2 = B**(-3.) * (D - 1.d0/(R+B))
        Sr3 = 1.d0/(2*B**2.) * (D**2. - 1.d0/(R+B)**2.)
        Sr4 = 1.d0/(3*B) * (D**3. - 1.d0/(R+B)**3.)
        Sr5 = 1.d0/4.d0 * (D**4. - 1.d0/(R+B)**4.)
        Sr = Sr1 * DSQRT(Sr6 - Sr7 + Sr2 + Sr3 + Sr4 + Sr5)
        
* escape velocity
        Q = DSQRT(1.0 / (1. + R/B))
        Vesc  = Vesc0 / Q
        Vhigh = 0.99  * Vesc     

* generate velocities:
 95     VXg(I) = Sr * GASDEV()          
        VYg(I) = Sr * GASDEV()          
        VZg(I) = Sr * GASDEV()

        Vtot = DSQRT (VXg(I)**2. + VYg(I)**2. + VZg(I)**2.)

* check if particle not too fast
        IF (Vtot .GT. Vhigh) THEN 
          NVhigh = NVhigh + 1 
          GOTO 95            ! particle too fast --> draw velocities again
!          print*, 'Elimination: NVhigh = ', NVhigh 
        ENDIF

* pair particles by symmetry
        Xg (I+1) = -Xg (I);  Yg (I+1) = -Yg(I);   Zg (I+1) = -Zg(I)  
        VXg(I+1) = -VXg(I);  VYg(I+1) = -VYg(I);  VZg(I+1) = -VZg(I) 

* update center-of-mass
        XCM = XCM + Xg(I);    YCM = YCM + Yg(I);    ZCM = ZCM + Zg(I)
        XCM = XCM + Xg(I+1);  YCM = YCM + Yg(I+1);  ZCM = ZCM + Zg(I+1)

* update center-of-mass velocity 
        VXCM=VXCM+VXg(I  ); VYCM=VYCM+VYg(I  ); VZCM=VZCM +VZg(I  )
        VXCM=VXCM+VXg(I+1); VYCM=VYCM+VYg(I+1); VZCM=VZCM +VZg(I+1)
 2    CONTINUE
  
!      PRINT*,'No. of close-to-escape-velocity eliminations: ', NVhigh  



* center-of-mass:
      XCM  = XCM/N;    YCM  = YCM/N;    ZCM = ZCM/N
!      PRINT*,'CENTER-OF-MASS [kpc]: ', XCM, YCM, ZCM
      VXCM = VXCM/N;   VYCM = VYCM/N;   VZCM = VZCM/N
!      PRINT*,'CENTER-OF-MASS [km/s]: ',VXCM*VUNIT,VYCM*VUNIT,VZCM*VUNIT

* correct particle positions to the center-of-mass frame:
      DO 3, I = Nother + 1, Nother + N
        Xg (I) = Xg (I)-XCM;   Yg(I)  = Yg (I)-YCM;  Zg(I) = Zg (I)-ZCM
        VXg(I) = VXg(I)-VXCM;  VYg(I) = VYg(I)-VYCM; VZg(I)= VZg(I)-VZCM
 3    CONTINUE

      END

************************************************************************
      SUBROUTINE GenHomSph (Nother, N, Rhs, D, YB2, Vesc0, Sr0)
************************************************************************
* generates particle positions within homogeneous sphere of radius Rhs
* generates particles velocities assuming these are test particles within
*   Plummer potential

      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'merge.inp'

      PARAMETER (PI=3.141592653589d0, PI2 = 2.d0*PI)
      REAL*4 RAND    
      COMMON /UNITS/ TUNIT,VUNIT,YVUNIT,GMUNIT,GMUNIT6      
      COMMON /TestP/  X(NX), Y(NX), Z(NX), VX(NX), VY(NX), VZ(NX)
      COMMON /GXIES2/ GMp, Bp, Bp2, YBp2, CTrBp, Sigma0p, Vesc0p, Sr0p,
     &                GMs, Bs, Bs2, YBs2, CTrBs, Sigma0s, Vesc0s, Sr0s,
     &                GMeff, GMStoP, GMpart1,GMpart2, Eps2,Eps2s,Eps2p, 
     &                Omega23, Dini, Vini, Trmin, DTrBp, DTrBs, Eps,
     &                RtrBp, RtrBs, BgasP, BgasS, BgasP2, BgasS2,
     &                YBgasP2, YBgasS2, CTrBgasP, CtrBgasS, Vesc0gasP,
     &                Vesc0gasS, Sr0gasP, Sr0gasS, DTrBgasP, DTrBgasS,
     &                RTrBgasP, RTrBgasS, Eps2gasP, Eps2gasS   

c     CALL SRAND(ISEED1)

      IF (N .EQ. 0) RETURN    

      NVhigh = 0  ! No. of repeatedly drawn particles (V > Vesc elimination)

      XCM  = 0.d0;   YCM  = 0.d0;   ZCM  = 0.d0    ! center-of-mass 
      VXCM = 0.d0;   VYCM = 0.d0;   VZCM = 0.d0    ! center-of-mass velocity 


      DO 2 I = Nother + 1, Nother + N, 2 

* generate angles
        PHI      = PI2   * DBLE (RAND(0))  
        COSTHETA = 2.d0  * DBLE (RAND(0)) - 1.d0
        SINTHETA = DSQRT(1.-COSTHETA*COSTHETA)

* generate radius       
 1913   R = GenPosHS0 (DBLE(Rand(0)), Rhs)    
!       IF (R. LT. 0.001) GOTO 1913

* generate 3D position
        X(I) = R * DCOS(PHI) * SINTHETA
        Y(I) = R * DSIN(PHI) * SINTHETA
        Z(I) = R * COSTHETA

* velocity dispersion Sr
        Q = 1.d0 + R * R * YB2    
        Sr    = Sr0   * DSQRT (6.d0 * (Q**(-0.5d0) - D) )
        Sr = 0.d0           ! 

* escape velocity     
        Vesc  = Vesc0 / Q**0.25d0
        Vhigh = 0.99   * Vesc     

        
* generate velocities:
 95     VX(I) = Sr * GASDEV()          
        VY(I) = Sr * GASDEV()          
        VZ(I) = Sr * GASDEV()

        Vtot = DSQRT (VX(I)**2. + VY(I)**2. + VZ(I)**2.)

* check if particle not too fast
        IF (Vtot .GT. Vhigh) THEN 
          NVhigh = NVhigh + 1 
          GOTO 95            ! particle too fast --> draw velocities again
!          print*, 'Elimination: NVhigh = ', NVhigh 
        ENDIF

* pair particles by symmetry
        X (I+1) = -X (I);  Y (I+1) = -Y (I);   Z(I+1) = - Z(I)  
        VX(I+1) = -VX(I);  VY(I+1) = -VY(I);  VZ(I+1) = -VZ(I) 

* update center-of-mass
        XCM = XCM + X(I);    YCM = YCM + Y(I);    ZCM = ZCM + Z(I)
        XCM = XCM + X(I+1);  YCM = YCM + Y(I+1);  ZCM = ZCM + Z(I+1)

* update center-of-mass velocity 
        VXCM = VXCM +VX(I  ); VYCM = VYCM +VY(I  ); VZCM = VZCM +VZ(I  )
        VXCM = VXCM +VX(I+1); VYCM = VYCM +VY(I+1); VZCM = VZCM +VZ(I+1)
 2    CONTINUE
  
!      PRINT*,'No. of close-to-escape-velocity eliminations: ', NVhigh  



* center-of-mass:
      XCM  = XCM/N;    YCM  = YCM/N;    ZCM = ZCM/N
!      PRINT*,'CENTER-OF-MASS [kpc]: ', XCM, YCM, ZCM
      VXCM = VXCM/N;   VYCM = VYCM/N;   VZCM = VZCM/N
!      PRINT*,'CENTER-OF-MASS [km/s]: ',VXCM*VUNIT,VYCM*VUNIT,VZCM*VUNIT

* correct particle positions to the center-of-mass frame:
      DO 3, I = Nother + 1, Nother + N
        X (I) = X (I)-XCM;   Y(I)  = Y (I)-YCM;   Z(I)  = Z (I)-ZCM
        VX(I) = VX(I)-VXCM;  VY(I) = VY(I)-VYCM;  VZ(I) = VZ(I)-VZCM
 3    CONTINUE

      END


************************************************************************
      SUBROUTINE GenShell (Nother, N, Rshell, YB2, Vesc0)
************************************************************************
* generates particle positions within a uniform shell of radius Rshell

      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'merge.inp'

      PARAMETER (PI=3.141592653589d0, PI2 = 2.d0*PI)
      REAL*4 RAND    
      COMMON /UNITS/ TUNIT,VUNIT,YVUNIT,GMUNIT,GMUNIT6      
      COMMON /TestP/  X(NX), Y(NX), Z(NX), VX(NX), VY(NX), VZ(NX)
      COMMON /GXIES2/ GMp, Bp, Bp2, YBp2, CTrBp, Sigma0p, Vesc0p, Sr0p,
     &                GMs, Bs, Bs2, YBs2, CTrBs, Sigma0s, Vesc0s, Sr0s,
     &                GMeff, GMStoP, GMpart1,GMpart2, Eps2,Eps2s,Eps2p, 
     &                Omega23, Dini, Vini, Trmin, DTrBp, DTrBs, Eps,
     &                RtrBp, RtrBs, BgasP, BgasS, BgasP2, BgasS2,
     &                YBgasP2, YBgasS2, CTrBgasP, CtrBgasS, Vesc0gasP,
     &                Vesc0gasS, Sr0gasP, Sr0gasS, DTrBgasP, DTrBgasS,
     &                RTrBgasP, RTrBgasS, Eps2gasP, Eps2gasS   

c     CALL SRAND(ISEED1)

      IF (N .EQ. 0) RETURN    

      NVhigh = 0  ! No. of repeatedly drawn particles (V > Vesc elimination)

      XCM  = 0.d0;   YCM  = 0.d0;   ZCM  = 0.d0    ! center-of-mass 
      VXCM = 0.d0;   VYCM = 0.d0;   VZCM = 0.d0    ! center-of-mass velocity 


      DO 2 I = Nother + 1, Nother + N, 2 

* generate angles
        PHI      = PI2   * DBLE (RAND(0))  
        COSTHETA = 2.d0  * DBLE (RAND(0)) - 1.d0
        SINTHETA = DSQRT(1.-COSTHETA*COSTHETA)

* generate radius       
        R = Rshell

* generate 3D position
        X(I) = R * DCOS(PHI) * SINTHETA
        Y(I) = R * DSIN(PHI) * SINTHETA
        Z(I) = R * COSTHETA

* velocity dispersion Sr
        Sr    = 0.d0

* escape velocity
        Q     = 1.d0 + R * R * YB2         
        Vesc  = Vesc0 / Q**0.25d0
        Vhigh = 0.9   * Vesc     

* generate velocities:
 95     VX(I) = Sr * GASDEV()          
        VY(I) = Sr * GASDEV()          
        VZ(I) = Sr * GASDEV()

        Vtot = DSQRT (VX(I)**2. + VY(I)**2. + VZ(I)**2.)

* check if particle not too fast
        IF (Vtot .GT. Vhigh) THEN 
          NVhigh = NVhigh + 1 
          GOTO 95            ! particle too fast --> draw velocities again
!          print*, 'Elimination: NVhigh = ', NVhigh 
        ENDIF

* pair particles by symmetry
        X (I+1) = -X (I);  Y (I+1) = -Y(I);   Z (I+1) = -Z(I)  
        VX(I+1) = -VX(I);  VY(I+1) = -VY(I);  VZ(I+1) = -VZ(I) 

* update center-of-mass
        XCM = XCM + X(I);    YCM = YCM + Y(I);    ZCM = ZCM + Z(I)
        XCM = XCM + X(I+1);  YCM = YCM + Y(I+1);  ZCM = ZCM + Z(I+1)

* update center-of-mass velocity 
        VXCM = VXCM +VX(I  ); VYCM = VYCM +VY(I  ); VZCM = VZCM +VZ(I  )
        VXCM = VXCM +VX(I+1); VYCM = VYCM +VY(I+1); VZCM = VZCM +VZ(I+1)
 2    CONTINUE
  
!      PRINT*,'No. of close-to-escape-velocity eliminations: ', NVhigh  



* center-of-mass:
      XCM  = XCM/N;    YCM  = YCM/N;    ZCM = ZCM/N
!      PRINT*,'CENTER-OF-MASS [kpc]: ', XCM, YCM, ZCM
      VXCM = VXCM/N;   VYCM = VYCM/N;   VZCM = VZCM/N
!      PRINT*,'CENTER-OF-MASS [km/s]: ',VXCM*VUNIT,VYCM*VUNIT,VZCM*VUNIT

* correct particle positions to the center-of-mass frame:
      DO 3, I = Nother + 1, Nother + N
        X (I) = X (I)-XCM;   Y(I)  = Y (I)-YCM;   Z(I)  = Z (I)-ZCM
        VX(I) = VX(I)-VXCM;  VY(I) = VY(I)-VYCM;  VZ(I) = VZ(I)-VZCM
 3    CONTINUE

      END

*************************************************************************
      SUBROUTINE EqMovG1P
********************************************************y*****************
! integrates eq. of motion of the two galaxies via leapfrog
! Primary's potential = Plummer

      IMPLICIT REAL*8 (A-H,O-Z)

      COMMON /GXIES/  Xp(3),Xs(3),Vp(3),Vs(3)
      COMMON /GXIES2/ GMp, Bp, Bp2, YBp2, CTrBp, Sigma0p, Vesc0p, Sr0p,
     &                GMs, Bs, Bs2, YBs2, CTrBs, Sigma0s, Vesc0s, Sr0s,
     &                GMeff, GMStoP, GMpart1,GMpart2, Eps2,Eps2s,Eps2p, 
     &                Omega23, Dini, Vini, Trmin, DTrBp, DTrBs, Eps,
     &                RtrBp, RtrBs, BgasP, BgasS, BgasP2, BgasS2,
     &                YBgasP2, YBgasS2, CTrBgasP, CtrBgasS, Vesc0gasP,
     &                Vesc0gasS, Sr0gasP, Sr0gasS, DTrBgasP, DTrBgasS,
     &                RTrBgasP, RTrBgasS, Eps2gasP, Eps2gasS   


      Rs2    = Xs(1)*Xs(1) + Xs(2)*Xs(2) + Xs(3)*Xs(3)

      DO 10 I = 1, 3

        As     = ACCPLUMM (GMeff, Xs(I), Rs2, Eps2)
        Vs (I) = Vs(I) + As                              ! update velocity
        Xs (I) = Xs(I) + Vs (I)                          ! update position

        Xp (I) = - GMStoP * Xs (I)                       ! adjust primary
        Vp (I) = - GMStoP * Vs (I)                       ! adjust primary

 10   CONTINUE

      END

*************************************************************************
      SUBROUTINE EqMovG1H
********************************************************y*****************
! integrates eq. of motion of the two galaxies via leapfrog
! Primary's potential = Hernquist

      IMPLICIT REAL*8 (A-H,O-Z)

      COMMON /GXIES/  Xp(3),Xs(3),Vp(3),Vs(3)
      COMMON /GXIES2/ GMp, Bp, Bp2, YBp2, CTrBp, Sigma0p, Vesc0p, Sr0p,
     &                GMs, Bs, Bs2, YBs2, CTrBs, Sigma0s, Vesc0s, Sr0s,
     &                GMeff, GMStoP, GMpart1,GMpart2, Eps2,Eps2s,Eps2p, 
     &                Omega23, Dini, Vini, Trmin, DTrBp, DTrBs, Eps,
     &                RtrBp, RtrBs, BgasP, BgasS, BgasP2, BgasS2,
     &                YBgasP2, YBgasS2, CTrBgasP, CtrBgasS, Vesc0gasP,
     &                Vesc0gasS, Sr0gasP, Sr0gasS, DTrBgasP, DTrBgasS,
     &                RTrBgasP, RTrBgasS, Eps2gasP, Eps2gasS  
     
      Rs = DSQRT(Xs(1)*Xs(1) + Xs(2)*Xs(2) + Xs(3)*Xs(3) ) 

      DO 10 I = 1, 3

        As     = ACCHERNQ (GMeff, Xs(I), Rs, Eps)
        Vs (I) = Vs(I) + As                              ! update velocity
        Xs (I) = Xs(I) + Vs (I)                          ! update position

        Xp (I) = - GMStoP * Xs (I)                       ! adjust primary
        Vp (I) = - GMStoP * Vs (I)                       ! adjust primary

 10   CONTINUE

      END

*************************************************************************
      SUBROUTINE EqMovG1V (Drel)
********************************************************y*****************
! integrates eq. of motion of the two galaxies via leapfrog
! like EqMovG1P but Primary's potential = de Vaucouleurs

      IMPLICIT REAL*8 (A-H,O-Z)

      COMMON /GXIES/  Xp(3),Xs(3),Vp(3),Vs(3)
      COMMON /GXIES2/ GMp, Bp, Bp2, YBp2, CTrBp, Sigma0p, Vesc0p, Sr0p,
     &                GMs, Bs, Bs2, YBs2, CTrBs, Sigma0s, Vesc0s, Sr0s,
     &                GMeff, GMStoP, GMpart1,GMpart2, Eps2,Eps2s,Eps2p, 
     &                Omega23, Dini, Vini, Trmin, DTrBp, DTrBs, Eps,
     &                RtrBp, RtrBs, BgasP, BgasS, BgasP2, BgasS2,
     &                YBgasP2, YBgasS2, CTrBgasP, CtrBgasS, Vesc0gasP,
     &                Vesc0gasS, Sr0gasP, Sr0gasS, DTrBgasP, DTrBgasS,
     &                RTrBgasP, RTrBgasS, Eps2gasP, Eps2gasS   


      Rs2    = Xs(1)*Xs(1) + Xs(2)*Xs(2) + Xs(3)*Xs(3) ! dist. of S from CM

      GMcorr = GMeff * Young76c(Drel/Bp)     ! correction of GMeff due to
                                             !   Plummer-->de Vaucouleurs
      DO 10 I = 1, 3


        As     = ACCPLUMM (GMcorr, Xs(I), Rs2, Eps2)
        Vs (I) = Vs(I) + As                              ! update velocity
        Xs (I) = Xs(I) + Vs (I)                          ! update position

        Xp (I) = - GMStoP * Xs (I)                       ! adjust primary
        Vp (I) = - GMStoP * Vs (I)                       ! adjust primary

 10   CONTINUE

      END

*************************************************************************
      SUBROUTINE EqMovG0P
*************************************************************************
! like EqMovG1P but only for initialization of leapfrog at t=0: 
!   integrates velocity of the two galaxies half-step backward 
                                   
      IMPLICIT REAL*8 (A-H,O-Z)

      COMMON /GXIES/  Xp(3),Xs(3),Vp(3),Vs(3)
      COMMON /GXIES2/ GMp, Bp, Bp2, YBp2, CTrBp, Sigma0p, Vesc0p, Sr0p,
     &                GMs, Bs, Bs2, YBs2, CTrBs, Sigma0s, Vesc0s, Sr0s,
     &                GMeff, GMStoP, GMpart1,GMpart2, Eps2,Eps2s,Eps2p, 
     &                Omega23, Dini, Vini, Trmin, DTrBp, DTrBs, Eps,
     &                RtrBp, RtrBs, BgasP, BgasS, BgasP2, BgasS2,
     &                YBgasP2, YBgasS2, CTrBgasP, CtrBgasS, Vesc0gasP,
     &                Vesc0gasS, Sr0gasP, Sr0gasS, DTrBgasP, DTrBgasS,
     &                RTrBgasP, RTrBgasS, Eps2gasP, Eps2gasS   


      Rs2    = Xs(1)*Xs(1) + Xs(2)*Xs(2) + Xs(3)*Xs(3)

      DO 10 I = 1, 3

        As     = ACCPLUMM (GMeff, Xs(I), Rs2, Eps2)
        Vs (I) = Vs(I) - 0.5d0 * As                      ! 

 10   CONTINUE

      END

*************************************************************************
      SUBROUTINE EqMovG0H
*************************************************************************
! like EqMovG1H but only for initialization of leapfrog at t=0: 
!   integrates velocity of the two galaxies half-step backward 
                                   
      IMPLICIT REAL*8 (A-H,O-Z)

      COMMON /GXIES/  Xp(3),Xs(3),Vp(3),Vs(3)
      COMMON /GXIES2/ GMp, Bp, Bp2, YBp2, CTrBp, Sigma0p, Vesc0p, Sr0p,
     &                GMs, Bs, Bs2, YBs2, CTrBs, Sigma0s, Vesc0s, Sr0s,
     &                GMeff, GMStoP, GMpart1,GMpart2, Eps2,Eps2s,Eps2p, 
     &                Omega23, Dini, Vini, Trmin, DTrBp, DTrBs, Eps,
     &                RtrBp, RtrBs, BgasP, BgasS, BgasP2, BgasS2,
     &                YBgasP2, YBgasS2, CTrBgasP, CtrBgasS, Vesc0gasP,
     &                Vesc0gasS, Sr0gasP, Sr0gasS, DTrBgasP, DTrBgasS,
     &                RTrBgasP, RTrBgasS, Eps2gasP, Eps2gasS   

      Rs   = DSQRT(Xs(1)*Xs(1) + Xs(2)*Xs(2) + Xs(3)*Xs(3))
      
      DO 10 I = 1, 3
        
        As     = ACCHERNQ (GMeff, Xs(I), Rs, Eps)
        Vs (I) = Vs(I) - 0.5d0 * As

 10   CONTINUE

      END
      
*************************************************************************
      SUBROUTINE EqMovG0V
*************************************************************************
! like EqMovG0P but Primary's potential = de Vaucouleurs
                                   
      IMPLICIT REAL*8 (A-H,O-Z)

      COMMON /GXIES/  Xp(3),Xs(3),Vp(3),Vs(3)
      COMMON /GXIES2/ GMp, Bp, Bp2, YBp2, CTrBp, Sigma0p, Vesc0p, Sr0p,
     &                GMs, Bs, Bs2, YBs2, CTrBs, Sigma0s, Vesc0s, Sr0s,
     &                GMeff, GMStoP, GMpart1,GMpart2, Eps2,Eps2s,Eps2p, 
     &                Omega23, Dini, Vini, Trmin, DTrBp, DTrBs, Eps,
     &                RtrBp, RtrBs, BgasP, BgasS, BgasP2, BgasS2,
     &                YBgasP2, YBgasS2, CTrBgasP, CtrBgasS, Vesc0gasP,
     &                Vesc0gasS, Sr0gasP, Sr0gasS, DTrBgasP, DTrBgasS,
     &                RTrBgasP, RTrBgasS, Eps2gasP, Eps2gasS   


      Rs2    = Xs(1)*Xs(1) + Xs(2)*Xs(2) + Xs(3)*Xs(3)

      GMcorr = GMeff * Young76c(Separation()/Bp) ! correction of GMeff due to
                                                 !   Plummer-->de Vaucouleurs

      DO 10 I = 1, 3

        As     = ACCPLUMM (GMcorr, Xs(I), Rs2, Eps2)
        Vs (I) = Vs(I) - 0.5d0 * As                      ! 

 10   CONTINUE

      END



*************************************************************************
      SUBROUTINE EqMovT1P
*************************************************************************
! integrates eq. of motion of test particles via leapfrog
! Primary's potential = Plummer
                                   
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'merge.inp'

      COMMON /GXIES/  Xp(3),Xs(3),Vp(3),Vs(3)
      COMMON /GXIES2/ GMp, Bp, Bp2, YBp2, CTrBp, Sigma0p, Vesc0p, Sr0p,
     &                GMs, Bs, Bs2, YBs2, CTrBs, Sigma0s, Vesc0s, Sr0s,
     &                GMeff, GMStoP, GMpart1,GMpart2, Eps2,Eps2s,Eps2p, 
     &                Omega23, Dini, Vini, Trmin, DTrBp, DTrBs, Eps,
     &                RtrBp, RtrBs, BgasP, BgasS, BgasP2, BgasS2,
     &                YBgasP2, YBgasS2, CTrBgasP, CtrBgasS, Vesc0gasP,
     &                Vesc0gasS, Sr0gasP, Sr0gasS, DTrBgasP, DTrBgasS,
     &                RTrBgasP, RTrBgasS, Eps2gasP, Eps2gasS 
      COMMON /TestP/  X(NX), Y(NX), Z(NX), VX(NX), VY(NX), VZ(NX)
      COMMON /TestP1/ N1, N2, NALL        
 
      DO 10 I = 1, NALL

        Rp2    = (X(I)-Xp(1))**2 + (Y(I)-Xp(2))**2 + (Z(I)-Xp(3))**2
        Rs2    = (X(I)-Xs(1))**2 + (Y(I)-Xs(2))**2 + (Z(I)-Xs(3))**2

        Ax     =   ACCPLUMM (GMs, X(I)-Xs(1), Rs2, Eps2s)
     &           + ACCPLUMM (GMp, X(I)-Xp(1), Rp2, Eps2p)
        Ay     =   ACCPLUMM (GMs, Y(I)-Xs(2), Rs2, Eps2s)
     &           + ACCPLUMM (GMp, Y(I)-Xp(2), Rp2, Eps2p)
        Az     =   ACCPLUMM (GMs, Z(I)-Xs(3), Rs2, Eps2s)
     &           + ACCPLUMM (GMp, Z(I)-Xp(3), Rp2, Eps2p)

        VX (I) = VX (I) + Ax                      
        VY (I) = VY (I) + Ay                      
        VZ (I) = VZ (I) + Az                      

        X  (I) = X (I) + VX(I)                      
        Y  (I) = Y (I) + VY(I)                   
        Z  (I) = Z (I) + VZ(I)                   

 10   CONTINUE

      END
      
*************************************************************************
      SUBROUTINE EqMovT1Pprint
*************************************************************************
! integrates eq. of motion of test particles via leapfrog
! Primary's potential = Plummer
                                   
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'merge.inp'

      COMMON /GXIES/  Xp(3),Xs(3),Vp(3),Vs(3)
      COMMON /GXIES2/ GMp, Bp, Bp2, YBp2, CTrBp, Sigma0p, Vesc0p, Sr0p,
     &                GMs, Bs, Bs2, YBs2, CTrBs, Sigma0s, Vesc0s, Sr0s,
     &                GMeff, GMStoP, GMpart1,GMpart2, Eps2,Eps2s,Eps2p, 
     &                Omega23, Dini, Vini, Trmin, DTrBp, DTrBs, Eps,
     &                RtrBp, RtrBs, BgasP, BgasS, BgasP2, BgasS2,
     &                YBgasP2, YBgasS2, CTrBgasP, CtrBgasS, Vesc0gasP,
     &                Vesc0gasS, Sr0gasP, Sr0gasS, DTrBgasP, DTrBgasS,
     &                RTrBgasP, RTrBgasS, Eps2gasP, Eps2gasS 
      COMMON /TestP/  X(NX), Y(NX), Z(NX), VX(NX), VY(NX), VZ(NX)
      COMMON /TestP1/ N1, N2, NALL
      
      
!      OPEN(80,FILE='positions_stars_plum.dat') 
! 800  FORMAT(f9.5, f9.5, f9.5, f9.5, f9.5, f9.5, f9.5)          
 
      DO 10 I = 1, NALL

        Rp2    = (X(I)-Xp(1))**2 + (Y(I)-Xp(2))**2 + (Z(I)-Xp(3))**2
        Rs2    = (X(I)-Xs(1))**2 + (Y(I)-Xs(2))**2 + (Z(I)-Xs(3))**2
        R      = DSQRT(X(I)**2 + Y(I)**2 + Z(I)**2) !distance from center of picture
        

        Ax     =   ACCPLUMM (GMs, X(I)-Xs(1), Rs2, Eps2s)
     &           + ACCPLUMM (GMp, X(I)-Xp(1), Rp2, Eps2p)
        Ay     =   ACCPLUMM (GMs, Y(I)-Xs(2), Rs2, Eps2s)
     &           + ACCPLUMM (GMp, Y(I)-Xp(2), Rp2, Eps2p)
        Az     =   ACCPLUMM (GMs, Z(I)-Xs(3), Rs2, Eps2s)
     &           + ACCPLUMM (GMp, Z(I)-Xp(3), Rp2, Eps2p)

        VX (I) = VX (I) + Ax                      
        VY (I) = VY (I) + Ay                      
        VZ (I) = VZ (I) + Az                      

        X  (I) = X (I) + VX(I)                      
        Y  (I) = Y (I) + VY(I)                   
        Z  (I) = Z (I) + VZ(I)     
        
!        WRITE(80,800) R, X(I), Y(I), Z(I), VX(I), VY(I), VZ(I)

        print*, R, X(I), Y(I), Z(I),  VX(I), VY(I), VZ(I)
        
 10   CONTINUE

      END
      
*************************************************************************
      SUBROUTINE EqMovTgas1P
*************************************************************************
! integrates eq. of motion of gas particles via leapfrog
! Primary's potential = Plummer
                                   
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'merge.inp'

      COMMON /GXIES/  Xp(3),Xs(3),Vp(3),Vs(3)
      COMMON /GXIES2/ GMp, Bp, Bp2, YBp2, CTrBp, Sigma0p, Vesc0p, Sr0p,
     &                GMs, Bs, Bs2, YBs2, CTrBs, Sigma0s, Vesc0s, Sr0s,
     &                GMeff, GMStoP, GMpart1,GMpart2, Eps2,Eps2s,Eps2p, 
     &                Omega23, Dini, Vini, Trmin, DTrBp, DTrBs, Eps,
     &                RtrBp, RtrBs, BgasP, BgasS, BgasP2, BgasS2,
     &                YBgasP2, YBgasS2, CTrBgasP, CtrBgasS, Vesc0gasP,
     &                Vesc0gasS, Sr0gasP, Sr0gasS, DTrBgasP, DTrBgasS,
     &                RTrBgasP, RTrBgasS, Eps2gasP, Eps2gasS 
      COMMON /Gas/    Xg(NX), Yg(NX), Zg(NX), VXg(NX), VYg(NX), VZg(NX),
     &                Xhalf(NX), Yhalf(NX), Zhalf(NX), Mesh(NX)
      COMMON /Gas1/   N1g, N2g, NALLg        
 
      DO 10 I = 1, NALLg

        Rp2    = (Xg(I)-Xp(1))**2 + (Yg(I)-Xp(2))**2 + (Zg(I)-Xp(3))**2
        Rs2    = (Xg(I)-Xs(1))**2 + (Yg(I)-Xs(2))**2 + (Zg(I)-Xs(3))**2
        R      = DSQRT(Xg(I)**2 + Yg(I)**2 + Zg(I)**2)

        Ax     =   ACCPLUMM (GMs, Xg(I)-Xs(1), Rs2, Eps2gasS)
     &           + ACCPLUMM (GMp, Xg(I)-Xp(1), Rp2, Eps2gasP)
        Ay     =   ACCPLUMM (GMs, Yg(I)-Xs(2), Rs2, Eps2gasS)
     &           + ACCPLUMM (GMp, Yg(I)-Xp(2), Rp2, Eps2gasP)
        Az     =   ACCPLUMM (GMs, Zg(I)-Xs(3), Rs2, Eps2gasS)
     &           + ACCPLUMM (GMp, Zg(I)-Xp(3), Rp2, Eps2gasP)

        VXg (I) = VXg (I) + Ax                      
        VYg (I) = VYg (I) + Ay                      
        VZg (I) = VZg (I) + Az                      

        Xtmp = Xg(I)   ! temporary positions saved for calculating 
        Ytmp = Yg(I)   ! position at mid-steps (same time as velocity)
        Ztmp = Zg(I)
        
        Xg  (I) = Xg (I) + VXg(I)                      
        Yg  (I) = Yg (I) + VYg(I)                   
        Zg  (I) = Zg (I) + VZg(I)                   
        
        Xhalf(I) = (Xtmp + Xg(I)) / 2.d0 ! positions at mid-step (same time
        Yhalf(I) = (Ytmp + Yg(I)) / 2.d0 ! as velocities) to calculate encounters 
        Zhalf(I) = (Ztmp + Zg(I)) / 2.d0 ! according to sticky particles scheme

c        PRINT*, R, Xg(I), Yg(I), Zg(I), VXg(I), VYg(I), VZg(I)
 10   CONTINUE

      END

*************************************************************************
      SUBROUTINE EqMovTgas1Pprint
*************************************************************************
! integrates eq. of motion of gas particles via leapfrog
! Primary's potential = Plummer
                                   
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'merge.inp'

      COMMON /GXIES/  Xp(3),Xs(3),Vp(3),Vs(3)
      COMMON /GXIES2/ GMp, Bp, Bp2, YBp2, CTrBp, Sigma0p, Vesc0p, Sr0p,
     &                GMs, Bs, Bs2, YBs2, CTrBs, Sigma0s, Vesc0s, Sr0s,
     &                GMeff, GMStoP, GMpart1,GMpart2, Eps2,Eps2s,Eps2p, 
     &                Omega23, Dini, Vini, Trmin, DTrBp, DTrBs, Eps,
     &                RtrBp, RtrBs, BgasP, BgasS, BgasP2, BgasS2,
     &                YBgasP2, YBgasS2, CTrBgasP, CtrBgasS, Vesc0gasP,
     &                Vesc0gasS, Sr0gasP, Sr0gasS, DTrBgasP, DTrBgasS,
     &                RTrBgasP, RTrBgasS, Eps2gasP, Eps2gasS 
      COMMON /Gas/    Xg(NX), Yg(NX), Zg(NX), VXg(NX), VYg(NX), VZg(NX),
     &                Xhalf(NX), Yhalf(NX), Zhalf(NX), Mesh(NX)
      COMMON /Gas1/   N1g, N2g, NALLg        
 
      DO 10 I = 1, NALLg

        Rp2    = (Xg(I)-Xp(1))**2 + (Yg(I)-Xp(2))**2 + (Zg(I)-Xp(3))**2
        Rs2    = (Xg(I)-Xs(1))**2 + (Yg(I)-Xs(2))**2 + (Zg(I)-Xs(3))**2
        R      = DSQRT(Xg(I)**2 + Yg(I)**2 + Zg(I)**2)

        Ax     =   ACCPLUMM (GMs, Xg(I)-Xs(1), Rs2, Eps2gasS)
     &           + ACCPLUMM (GMp, Xg(I)-Xp(1), Rp2, Eps2gasP)
        Ay     =   ACCPLUMM (GMs, Yg(I)-Xs(2), Rs2, Eps2gasS)
     &           + ACCPLUMM (GMp, Yg(I)-Xp(2), Rp2, Eps2gasP)
        Az     =   ACCPLUMM (GMs, Zg(I)-Xs(3), Rs2, Eps2gasS)
     &           + ACCPLUMM (GMp, Zg(I)-Xp(3), Rp2, Eps2gasP)

        VXg (I) = VXg (I) + Ax                      
        VYg (I) = VYg (I) + Ay                      
        VZg (I) = VZg (I) + Az                      

        Xtmp = Xg(I)   ! temporary positions saved for calculating 
        Ytmp = Yg(I)   ! position at mid-steps (same time as velocity)
        Ztmp = Zg(I)
        
        Xg  (I) = Xg (I) + VXg(I)                      
        Yg  (I) = Yg (I) + VYg(I)                   
        Zg  (I) = Zg (I) + VZg(I)                   
        
        Xhalf(I) = (Xtmp + Xg(I)) / 2.d0 ! positions at mid-step (same time
        Yhalf(I) = (Ytmp + Yg(I)) / 2.d0 ! as velocities) to calculate encounters 
        Zhalf(I) = (Ztmp + Zg(I)) / 2.d0 ! according to sticky particles scheme

        PRINT*, R, Xg(I), Yg(I), Zg(I), VXg(I), VYg(I), VZg(I)
 10   CONTINUE

      END

*************************************************************************
      SUBROUTINE EqMovT1Hprint
*************************************************************************
! integrates eq. of motion of test particles via leapfrog
! Primary's potential = Hernquist
                                   
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'merge.inp'

      COMMON /GXIES/  Xp(3),Xs(3),Vp(3),Vs(3)
      COMMON /GXIES2/ GMp, Bp, Bp2, YBp2, CTrBp, Sigma0p, Vesc0p, Sr0p,
     &                GMs, Bs, Bs2, YBs2, CTrBs, Sigma0s, Vesc0s, Sr0s,
     &                GMeff, GMStoP, GMpart1,GMpart2, Eps2,Eps2s,Eps2p, 
     &                Omega23, Dini, Vini, Trmin, DTrBp, DTrBs, Eps,
     &                RtrBp, RtrBs, BgasP, BgasS, BgasP2, BgasS2,
     &                YBgasP2, YBgasS2, CTrBgasP, CtrBgasS, Vesc0gasP,
     &                Vesc0gasS, Sr0gasP, Sr0gasS, DTrBgasP, DTrBgasS,
     &                RTrBgasP, RTrBgasS, Eps2gasP, Eps2gasS   
      COMMON /TestP/  X(NX), Y(NX), Z(NX), VX(NX), VY(NX), VZ(NX)
      COMMON /TestP1/ N1, N2, NALL 
      
      
!      OPEN(80,FILE='positions_stars_hern.dat') 
! 800  FORMAT(f9.5, f9.5, f9.5, f9.5, f9.5, f9.5, f9.5)      
 
      DO 10 I = 1, NALL
      
        Rp = DSQRT((X(I)-Xp(1))**2 + (Y(I)-Xp(2))**2 + (Z(I)-Xp(3))**2)
        Rs = DSQRT((X(I)-Xs(1))**2 + (Y(I)-Xs(2))**2 + (Z(I)-Xs(3))**2)
        R      = DSQRT(X(I)**2 + Y(I)**2 + Z(I)**2) !distance from center of picture
        

        Ax     =   ACCHERNQ (GMs, X(I)-Xs(1), Rs, Bs)
     &           + ACCHERNQ (GMp, X(I)-Xp(1), Rp, Bp)
        Ay     =   ACCHERNQ (GMs, Y(I)-Xs(2), Rs, Bs)
     &           + ACCHERNQ (GMp, Y(I)-Xp(2), Rp, Bp)
        Az     =   ACCHERNQ (GMs, Z(I)-Xs(3), Rs, Bs)
     &           + ACCHERNQ (GMp, Z(I)-Xp(3), Rp, Bp)

        VX (I) = VX (I) + Ax                      
        VY (I) = VY (I) + Ay                      
        VZ (I) = VZ (I) + Az                      

        X  (I) = X (I) + VX(I)                      
        Y  (I) = Y (I) + VY(I)                   
        Z  (I) = Z (I) + VZ(I)
        
!        WRITE(80,800) R, X(I), Y(I), Z(I), VX(I), VY(I), VZ(I)
        print*, R, X(I), Y(I), Z(I), VX(I), VY(I), VZ(I)

 10   CONTINUE

      END
*************************************************************************
      SUBROUTINE EqMovT1H
*************************************************************************
! integrates eq. of motion of test particles via leapfrog
! Primary's potential = Hernquist
                                   
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'merge.inp'

      COMMON /GXIES/  Xp(3),Xs(3),Vp(3),Vs(3)
      COMMON /GXIES2/ GMp, Bp, Bp2, YBp2, CTrBp, Sigma0p, Vesc0p, Sr0p,
     &                GMs, Bs, Bs2, YBs2, CTrBs, Sigma0s, Vesc0s, Sr0s,
     &                GMeff, GMStoP, GMpart1,GMpart2, Eps2,Eps2s,Eps2p, 
     &                Omega23, Dini, Vini, Trmin, DTrBp, DTrBs, Eps,
     &                RtrBp, RtrBs, BgasP, BgasS, BgasP2, BgasS2,
     &                YBgasP2, YBgasS2, CTrBgasP, CtrBgasS, Vesc0gasP,
     &                Vesc0gasS, Sr0gasP, Sr0gasS, DTrBgasP, DTrBgasS,
     &                RTrBgasP, RTrBgasS, Eps2gasP, Eps2gasS   
      COMMON /TestP/  X(NX), Y(NX), Z(NX), VX(NX), VY(NX), VZ(NX)
      COMMON /TestP1/ N1, N2, NALL        
 
      DO 10 I = 1, NALL
      
        Rp = DSQRT((X(I)-Xp(1))**2 + (Y(I)-Xp(2))**2 + (Z(I)-Xp(3))**2)
        Rs = DSQRT((X(I)-Xs(1))**2 + (Y(I)-Xs(2))**2 + (Z(I)-Xs(3))**2)

        Ax     =   ACCHERNQ (GMs, X(I)-Xs(1), Rs, Bs)
     &           + ACCHERNQ (GMp, X(I)-Xp(1), Rp, Bp)
        Ay     =   ACCHERNQ (GMs, Y(I)-Xs(2), Rs, Bs)
     &           + ACCHERNQ (GMp, Y(I)-Xp(2), Rp, Bp)
        Az     =   ACCHERNQ (GMs, Z(I)-Xs(3), Rs, Bs)
     &           + ACCHERNQ (GMp, Z(I)-Xp(3), Rp, Bp)

        VX (I) = VX (I) + Ax                      
        VY (I) = VY (I) + Ay                      
        VZ (I) = VZ (I) + Az                      

        X  (I) = X (I) + VX(I)                      
        Y  (I) = Y (I) + VY(I)                   
        Z  (I) = Z (I) + VZ(I)

 10   CONTINUE

      END
           
*************************************************************************
      SUBROUTINE EqMovTgas1H
*************************************************************************
! integrates eq. of motion for gas particles via leapfrog
! Primary's potential = Hernquist
                                   
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'merge.inp'

      COMMON /GXIES/  Xp(3),Xs(3),Vp(3),Vs(3)
      COMMON /GXIES2/ GMp, Bp, Bp2, YBp2, CTrBp, Sigma0p, Vesc0p, Sr0p,
     &                GMs, Bs, Bs2, YBs2, CTrBs, Sigma0s, Vesc0s, Sr0s,
     &                GMeff, GMStoP, GMpart1,GMpart2, Eps2,Eps2s,Eps2p, 
     &                Omega23, Dini, Vini, Trmin, DTrBp, DTrBs, Eps,
     &                RtrBp, RtrBs, BgasP, BgasS, BgasP2, BgasS2,
     &                YBgasP2, YBgasS2, CTrBgasP, CtrBgasS, Vesc0gasP,
     &                Vesc0gasS, Sr0gasP, Sr0gasS, DTrBgasP, DTrBgasS,
     &                RTrBgasP, RTrBgasS, Eps2gasP, Eps2gasS 
      COMMON /Gas/    Xg(NX), Yg(NX), Zg(NX), VXg(NX), VYg(NX), VZg(NX),
     &                Xhalf(NX), Yhalf(NX), Zhalf(NX), Mesh(NX)
      COMMON /Gas1/   N1g, N2g, NALLg        
 
      DO 10 I = 1, NALLg

        Rp = DSQRT((Xg(I)-Xp(1))**2 + (Yg(I)-Xp(2))**2
     &               + (Zg(I)-Xp(3))**2)
        Rs = DSQRT((Xg(I)-Xs(1))**2 + (Yg(I)-Xs(2))**2
     &               + (Zg(I)-Xs(3))**2)
        R  = DSQRT(Xg(I)**2 + Yg(I)**2 + Zg(I)**2)

        Ax     =   ACCHERNQ (GMs, Xg(I)-Xs(1), Rs, BgasS)
     &           + ACCHERNQ (GMp, Xg(I)-Xp(1), Rp, BgasP)
        Ay     =   ACCHERNQ (GMs, Yg(I)-Xs(2), Rs, BgasS)
     &           + ACCHERNQ (GMp, Yg(I)-Xp(2), Rp, BgasP)
        Az     =   ACCHERNQ (GMs, Zg(I)-Xs(3), Rs, BgasS)
     &           + ACCHERNQ (GMp, Zg(I)-Xp(3), Rp, BgasP)

        VXg(I) = VXg(I) + Ax                      
        VYg(I) = VYg(I) + Ay                      
        VZg(I) = VZg(I) + Az                      

        Xtmp = Xg(I)   ! temporary positions saved for calculating 
        Ytmp = Yg(I)   ! position at mid-steps (same time as velocity)
        Ztmp = Zg(I)
        
        Xg(I) = Xg(I) + VXg(I)                      
        Yg(I) = Yg(I) + VYg(I)                   
        Zg(I) = Zg(I) + VZg(I)                   
        
        Xhalf(I) = (Xtmp + Xg(I)) / 2.d0 ! positions at mid-step (same time
        Yhalf(I) = (Ytmp + Yg(I)) / 2.d0 ! as velocities) to calculate encounters 
        Zhalf(I) = (Ztmp + Zg(I)) / 2.d0 ! according to sticky particles scheme

c        PRINT*, R, Xg(I), Yg(I), Zg(I), VXg(I), VYg(I), VZg(I)
 10   CONTINUE

      END

*************************************************************************
      SUBROUTINE EqMovTgas1Hprint
*************************************************************************
! integrates eq. of motion for gas particles via leapfrog
! Primary's potential = Hernquist
                                   
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'merge.inp'

      COMMON /GXIES/  Xp(3),Xs(3),Vp(3),Vs(3)
      COMMON /GXIES2/ GMp, Bp, Bp2, YBp2, CTrBp, Sigma0p, Vesc0p, Sr0p,
     &                GMs, Bs, Bs2, YBs2, CTrBs, Sigma0s, Vesc0s, Sr0s,
     &                GMeff, GMStoP, GMpart1,GMpart2, Eps2,Eps2s,Eps2p, 
     &                Omega23, Dini, Vini, Trmin, DTrBp, DTrBs, Eps,
     &                RtrBp, RtrBs, BgasP, BgasS, BgasP2, BgasS2,
     &                YBgasP2, YBgasS2, CTrBgasP, CtrBgasS, Vesc0gasP,
     &                Vesc0gasS, Sr0gasP, Sr0gasS, DTrBgasP, DTrBgasS,
     &                RTrBgasP, RTrBgasS, Eps2gasP, Eps2gasS 
      COMMON /Gas/    Xg(NX), Yg(NX), Zg(NX), VXg(NX), VYg(NX), VZg(NX),
     &                Xhalf(NX), Yhalf(NX), Zhalf(NX), Mesh(NX)
      COMMON /Gas1/   N1g, N2g, NALLg        
 
      DO 10 I = 1, NALLg

        Rp = DSQRT((Xg(I)-Xp(1))**2 + (Yg(I)-Xp(2))**2
     &               + (Zg(I)-Xp(3))**2)
        Rs = DSQRT((Xg(I)-Xs(1))**2 + (Yg(I)-Xs(2))**2
     &               + (Zg(I)-Xs(3))**2)
        R  = DSQRT(Xg(I)**2 + Yg(I)**2 + Zg(I)**2)

        Ax     =   ACCHERNQ (GMs, Xg(I)-Xs(1), Rs, BgasS)
     &           + ACCHERNQ (GMp, Xg(I)-Xp(1), Rp, BgasP)
        Ay     =   ACCHERNQ (GMs, Yg(I)-Xs(2), Rs, BgasS)
     &           + ACCHERNQ (GMp, Yg(I)-Xp(2), Rp, BgasP)
        Az     =   ACCHERNQ (GMs, Zg(I)-Xs(3), Rs, BgasS)
     &           + ACCHERNQ (GMp, Zg(I)-Xp(3), Rp, BgasP)

        VXg(I) = VXg(I) + Ax                      
        VYg(I) = VYg(I) + Ay                      
        VZg(I) = VZg(I) + Az                      

        Xtmp = Xg(I)   ! temporary positions saved for calculating 
        Ytmp = Yg(I)   ! position at mid-steps (same time as velocity)
        Ztmp = Zg(I)
        
        Xg(I) = Xg(I) + VXg(I)                      
        Yg(I) = Yg(I) + VYg(I)                   
        Zg(I) = Zg(I) + VZg(I)                   
        
        Xhalf(I) = (Xtmp + Xg(I)) / 2.d0 ! positions at mid-step (same time
        Yhalf(I) = (Ytmp + Yg(I)) / 2.d0 ! as velocities) to calculate encounters 
        Zhalf(I) = (Ztmp + Zg(I)) / 2.d0 ! according to sticky particles scheme

        PRINT*, R, Xg(I), Yg(I), Zg(I), VXg(I), VYg(I), VZg(I)
 10   CONTINUE

      END

*************************************************************************
      SUBROUTINE EqMovT1V
*************************************************************************
! integrates eq. of motion of test particles via leapfrog
! Primary's potential = de Vaucouleurs
                                   
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'merge.inp'

      COMMON /GXIES/  Xp(3),Xs(3),Vp(3),Vs(3)
      COMMON /GXIES2/ GMp, Bp, Bp2, YBp2, CTrBp, Sigma0p, Vesc0p, Sr0p,
     &                GMs, Bs, Bs2, YBs2, CTrBs, Sigma0s, Vesc0s, Sr0s,
     &                GMeff, GMStoP, GMpart1,GMpart2, Eps2,Eps2s,Eps2p, 
     &                Omega23, Dini, Vini, Trmin, DTrBp, DTrBs, Eps,
     &                RtrBp, RtrBs, BgasP, BgasS, BgasP2, BgasS2,
     &                YBgasP2, YBgasS2, CTrBgasP, CtrBgasS, Vesc0gasP,
     &                Vesc0gasS, Sr0gasP, Sr0gasS, DTrBgasP, DTrBgasS,
     &                RTrBgasP, RTrBgasS, Eps2gasP, Eps2gasS   
      COMMON /TestP/  X(NX), Y(NX), Z(NX), VX(NX), VY(NX), VZ(NX)
      COMMON /TestP1/ N1, N2, NALL        
 
      DO 10 I = 1, NALL

        Rp2    = (X(I)-Xp(1))**2 + (Y(I)-Xp(2))**2 + (Z(I)-Xp(3))**2
        Rs2    = (X(I)-Xs(1))**2 + (Y(I)-Xs(2))**2 + (Z(I)-Xs(3))**2

        GMpCorr = GMp * Young76c(DSQRT(Rp2)/Bp)  ! correction of GMp due to
                                                 !   Plummer-->de Vaucouleurs

        Ax     =   ACCPLUMM (GMs,     X(I)-Xs(1), Rs2, Eps2s)
     &           + ACCPLUMM (GMpCorr, X(I)-Xp(1), Rp2, Eps2p)
        Ay     =   ACCPLUMM (GMs,     Y(I)-Xs(2), Rs2, Eps2s)
     &           + ACCPLUMM (GMpCorr, Y(I)-Xp(2), Rp2, Eps2p)
        Az     =   ACCPLUMM (GMs,     Z(I)-Xs(3), Rs2, Eps2s)
     &           + ACCPLUMM (GMpCorr, Z(I)-Xp(3), Rp2, Eps2p)

        VX (I) = VX (I) + Ax                      
        VY (I) = VY (I) + Ay                      
        VZ (I) = VZ (I) + Az                      

        X  (I) = X (I) + VX(I)                      
        Y  (I) = Y (I) + VY(I)                   
        Z  (I) = Z (I) + VZ(I)                   

 10   CONTINUE

      END


*************************************************************************
      SUBROUTINE EqMovT0P
*************************************************************************
! like EqMovT1P but only for initialization of leapfrog at t=0: 
!   integrates velocity of test particles half-step backward 

                                   
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'merge.inp'

      COMMON /GXIES/  Xp(3),Xs(3),Vp(3),Vs(3)
      COMMON /GXIES2/ GMp, Bp, Bp2, YBp2, CTrBp, Sigma0p, Vesc0p, Sr0p,
     &                GMs, Bs, Bs2, YBs2, CTrBs, Sigma0s, Vesc0s, Sr0s,
     &                GMeff, GMStoP, GMpart1,GMpart2, Eps2,Eps2s,Eps2p, 
     &                Omega23, Dini, Vini, Trmin, DTrBp, DTrBs, Eps,
     &                RtrBp, RtrBs, BgasP, BgasS, BgasP2, BgasS2,
     &                YBgasP2, YBgasS2, CTrBgasP, CtrBgasS, Vesc0gasP,
     &                Vesc0gasS, Sr0gasP, Sr0gasS, DTrBgasP, DTrBgasS,
     &                RTrBgasP, RTrBgasS, Eps2gasP, Eps2gasS   
      COMMON /TestP/  X(NX), Y(NX), Z(NX), VX(NX), VY(NX), VZ(NX)
      COMMON /TestP1/ N1, N2, NALL        

      DO 10 I = 1, NALL

        Rp2  = (X(I)-Xp(1))**2 + (Y(I)-Xp(2))**2 + (Z(I)-Xp(3))**2
        Rs2  = (X(I)-Xs(1))**2 + (Y(I)-Xs(2))**2 + (Z(I)-Xs(3))**2

        Ax     =   ACCPLUMM (GMs, X(I)-Xs(1), Rs2, Eps2s)
     &           + ACCPLUMM (GMp, X(I)-Xp(1), Rp2, Eps2p)
        Ay     =   ACCPLUMM (GMs, Y(I)-Xs(2), Rs2, Eps2s)
     &           + ACCPLUMM (GMp, Y(I)-Xp(2), Rp2, Eps2p)
        Az     =   ACCPLUMM (GMs, Z(I)-Xs(3), Rs2, Eps2s)
     &           + ACCPLUMM (GMp, Z(I)-Xp(3), Rp2, Eps2p)

        VX (I) = VX (I) - 0.5d0 * Ax                      
        VY (I) = VY (I) - 0.5d0 * Ay                      
        VZ (I) = VZ (I) - 0.5d0 * Az                      

 10   CONTINUE

      END
      
*************************************************************************
      SUBROUTINE EqMovTgas0P
*************************************************************************
! like EqMovTgas1P but only for initialization of leapfrog at t=0: 
!   integrates velocity of test particles half-step backward 
                                   
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'merge.inp'

      COMMON /GXIES/  Xp(3),Xs(3),Vp(3),Vs(3)
      COMMON /GXIES2/ GMp, Bp, Bp2, YBp2, CTrBp, Sigma0p, Vesc0p, Sr0p,
     &                GMs, Bs, Bs2, YBs2, CTrBs, Sigma0s, Vesc0s, Sr0s,
     &                GMeff, GMStoP, GMpart1,GMpart2, Eps2,Eps2s,Eps2p, 
     &                Omega23, Dini, Vini, Trmin, DTrBp, DTrBs, Eps,
     &                RtrBp, RtrBs, BgasP, BgasS, BgasP2, BgasS2,
     &                YBgasP2, YBgasS2, CTrBgasP, CtrBgasS, Vesc0gasP,
     &                Vesc0gasS, Sr0gasP, Sr0gasS, DTrBgasP, DTrBgasS,
     &                RTrBgasP, RTrBgasS, Eps2gasP, Eps2gasS 
      COMMON /Gas/    Xg(NX), Yg(NX), Zg(NX), VXg(NX), VYg(NX), VZg(NX),
     &                Xhalf(NX), Yhalf(NX), Zhalf(NX), Mesh(NX)
      COMMON /Gas1/   N1g, N2g, NALLg        
 
      DO 10 I = 1, NALLg

        Rp2    = (Xg(I)-Xp(1))**2 + (Yg(I)-Xp(2))**2 + (Zg(I)-Xp(3))**2
        Rs2    = (Xg(I)-Xs(1))**2 + (Yg(I)-Xs(2))**2 + (Zg(I)-Xs(3))**2

        Ax     =   ACCPLUMM (GMs, Xg(I)-Xs(1), Rs2, Eps2gasS)
     &           + ACCPLUMM (GMp, Xg(I)-Xp(1), Rp2, Eps2gasP)
        Ay     =   ACCPLUMM (GMs, Yg(I)-Xs(2), Rs2, Eps2gasS)
     &           + ACCPLUMM (GMp, Yg(I)-Xp(2), Rp2, Eps2gasP)
        Az     =   ACCPLUMM (GMs, Zg(I)-Xs(3), Rs2, Eps2gasS)
     &           + ACCPLUMM (GMp, Zg(I)-Xp(3), Rp2, Eps2gasP)

        VXg (I) = VXg (I) - 0.5d0 * Ax                      
        VYg (I) = VYg (I) - 0.5d0 * Ay                      
        VZg (I) = VZg (I) - 0.5d0 * Az  
        
        !print*, Xg(I), Yg(I), Zg(I), Vxg(I), vyg(i), vzg(i)                    

 10   CONTINUE

      END

      
*************************************************************************
      SUBROUTINE EqMovT0H
*************************************************************************
! like EqMovT1H but only for initialization of leapfrog at t=0: 
!   integrates velocity of test particles half-step backward 

                                   
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'merge.inp'

      COMMON /GXIES/  Xp(3),Xs(3),Vp(3),Vs(3)
      COMMON /GXIES2/ GMp, Bp, Bp2, YBp2, CTrBp, Sigma0p, Vesc0p, Sr0p,
     &                GMs, Bs, Bs2, YBs2, CTrBs, Sigma0s, Vesc0s, Sr0s,
     &                GMeff, GMStoP, GMpart1,GMpart2, Eps2,Eps2s,Eps2p, 
     &                Omega23, Dini, Vini, Trmin, DTrBp, DTrBs, Eps,
     &                RtrBp, RtrBs, BgasP, BgasS, BgasP2, BgasS2,
     &                YBgasP2, YBgasS2, CTrBgasP, CtrBgasS, Vesc0gasP,
     &                Vesc0gasS, Sr0gasP, Sr0gasS, DTrBgasP, DTrBgasS,
     &                RTrBgasP, RTrBgasS, Eps2gasP, Eps2gasS   
      COMMON /TestP/  X(NX), Y(NX), Z(NX), VX(NX), VY(NX), VZ(NX)
      COMMON /TestP1/ N1, N2, NALL        

      DO 10 I = 1, NALL
      
        Rp  = DSQRT((X(I)-Xp(1))**2 + (Y(I)-Xp(2))**2 + (Z(I)-Xp(3))**2)
        Rs  = DSQRT((X(I)-Xs(1))**2 + (Y(I)-Xs(2))**2 + (Z(I)-Xs(3))**2)

        Ax     =   ACCHERNQ (GMs, X(I)-Xs(1), Rs, Bs)
     &           + ACCHERNQ (GMp, X(I)-Xp(1), Rp, Bp)
        Ay     =   ACCHERNQ (GMs, Y(I)-Xs(2), Rs, Bs)
     &           + ACCHERNQ (GMp, Y(I)-Xp(2), Rp, Bp)
        Az     =   ACCHERNQ (GMs, Z(I)-Xs(3), Rs, Bs)
     &           + ACCHERNQ (GMp, Z(I)-Xp(3), Rp, Bp)

        VX (I) = VX (I) - 0.5d0 * Ax                      
        VY (I) = VY (I) - 0.5d0 * Ay                      
        VZ (I) = VZ (I) - 0.5d0 * Az                      

 10   CONTINUE

      END
      
*************************************************************************
      SUBROUTINE EqMovTgas0H
*************************************************************************
! like EqMovTgas1H but only for initialization of leapfrog at t=0: 
!   integrates velocity of test particles half-step backward 
                                   
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'merge.inp'

      COMMON /GXIES/  Xp(3),Xs(3),Vp(3),Vs(3)
      COMMON /GXIES2/ GMp, Bp, Bp2, YBp2, CTrBp, Sigma0p, Vesc0p, Sr0p,
     &                GMs, Bs, Bs2, YBs2, CTrBs, Sigma0s, Vesc0s, Sr0s,
     &                GMeff, GMStoP, GMpart1,GMpart2, Eps2,Eps2s,Eps2p, 
     &                Omega23, Dini, Vini, Trmin, DTrBp, DTrBs, Eps,
     &                RtrBp, RtrBs, BgasP, BgasS, BgasP2, BgasS2,
     &                YBgasP2, YBgasS2, CTrBgasP, CtrBgasS, Vesc0gasP,
     &                Vesc0gasS, Sr0gasP, Sr0gasS, DTrBgasP, DTrBgasS,
     &                RTrBgasP, RTrBgasS, Eps2gasP, Eps2gasS 
      COMMON /Gas/    Xg(NX), Yg(NX), Zg(NX), VXg(NX), VYg(NX), VZg(NX),
     &                Xhalf(NX), Yhalf(NX), Zhalf(NX), Mesh(NX)
      COMMON /Gas1/   N1g, N2g, NALLg        
 
      DO 10 I = 1, NALLg

        Rp = DSQRT((Xg(I)-Xp(1))**2 + (Yg(I)-Xp(2))**2
     &        + (Zg(I)-Xp(3))**2)
        Rs = DSQRT((Xg(I)-Xs(1))**2 + (Yg(I)-Xs(2))**2
     &        + (Zg(I)-Xs(3))**2)

        Ax     =   ACCHERNQ (GMs, Xg(I)-Xs(1), Rs, BgasS)
     &           + ACCHERNQ (GMp, Xg(I)-Xp(1), Rp, BgasP)
        Ay     =   ACCHERNQ (GMs, Yg(I)-Xs(2), Rs, BgasS)
     &           + ACCHERNQ (GMp, Yg(I)-Xp(2), Rp, BgasP)
        Az     =   ACCHERNQ (GMs, Zg(I)-Xs(3), Rs, BgasS)
     &           + ACCHERNQ (GMp, Zg(I)-Xp(3), Rp, BgasP)

        VXg (I) = VXg (I) - 0.5d0 * Ax                      
        VYg (I) = VYg (I) - 0.5d0 * Ay                      
        VZg (I) = VZg (I) - 0.5d0 * Az
        
        

 10   CONTINUE

      END

*************************************************************************
      SUBROUTINE EqMovT0V
*************************************************************************
! like EqMovT0P but Primary's potential = deVacoulers
                                   
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'merge.inp'

      COMMON /GXIES/  Xp(3),Xs(3),Vp(3),Vs(3)
      COMMON /GXIES2/ GMp, Bp, Bp2, YBp2, CTrBp, Sigma0p, Vesc0p, Sr0p,
     &                GMs, Bs, Bs2, YBs2, CTrBs, Sigma0s, Vesc0s, Sr0s,
     &                GMeff, GMStoP, GMpart1,GMpart2, Eps2,Eps2s,Eps2p, 
     &                Omega23, Dini, Vini, Trmin, DTrBp, DTrBs, Eps,
     &                RtrBp, RtrBs, BgasP, BgasS, BgasP2, BgasS2,
     &                YBgasP2, YBgasS2, CTrBgasP, CtrBgasS, Vesc0gasP,
     &                Vesc0gasS, Sr0gasP, Sr0gasS, DTrBgasP, DTrBgasS,
     &                RTrBgasP, RTrBgasS, Eps2gasP, Eps2gasS   
      COMMON /TestP/  X(NX), Y(NX), Z(NX), VX(NX), VY(NX), VZ(NX)
      COMMON /TestP1/ N1, N2, NALL        

      DO 10 I = 1, NALL

        Rp2    = (X(I)-Xp(1))**2 + (Y(I)-Xp(2))**2 + (Z(I)-Xp(3))**2
        Rs2    = (X(I)-Xs(1))**2 + (Y(I)-Xs(2))**2 + (Z(I)-Xs(3))**2

        GMpCorr = GMp * Young76c(DSQRT(Rp2)/Bp)  ! correction of GMp due to
                                                 !   Plummer-->de Vaucouleurs

        Ax     =   ACCPLUMM (GMs,     X(I)-Xs(1), Rs2, Eps2s)
     &           + ACCPLUMM (GMpCorr, X(I)-Xp(1), Rp2, Eps2p)
        Ay     =   ACCPLUMM (GMs,     Y(I)-Xs(2), Rs2, Eps2s)
     &           + ACCPLUMM (GMpCorr, Y(I)-Xp(2), Rp2, Eps2p)
        Az     =   ACCPLUMM (GMs,     Z(I)-Xs(3), Rs2, Eps2s)
     &           + ACCPLUMM (GMpCorr, Z(I)-Xp(3), Rp2, Eps2p)

        VX (I) = VX (I) - 0.5d0 * Ax                      
        VY (I) = VY (I) - 0.5d0 * Ay                      
        VZ (I) = VZ (I) - 0.5d0 * Az                      

 10   CONTINUE

      END

************************************************************************
      SUBROUTINE StickyParticles
************************************************************************
! test for close gas particles and change their relative velocities
! accordingly to sticky particles scheme

      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'merge.inp'
      
      integer:: COL=0 
      COMMON /Gas/    Xg(NX), Yg(NX), Zg(NX), VXg(NX), VYg(NX), VZg(NX),
     &                Xhalf(NX), Yhalf(NX), Zhalf(NX), Mesh(NX)
      COMMON /Gas1/   N1g, N2g, NALLg 
       
      ! Parameters fo sticky particles scheme
      BetaT = 1.    ! changes tangential part of velocity (1 to conserve angular momentum) 
      BetaR = 0.5d0 ! changes radial part of velocity (1 for ellastic collision)
      Rcl = 0.02d0     ! diameter of cloud in kpc
      
      Dcl = 2.0*Rcl
      
c      MeshNr = 0     !Nr. indexing meshgrid for particles
c      DO 10 I = 1, NALLg
c      
c        DO 20 J = -100, 100
c          
c          DO 30 K = -100, 100
c            IF ((Xg(I) .GT. J) .AND. (Xg(I) .LT. J+1) .AND.
c     &         (Yg(I) .GT. K) .AND. (Yg(I) .LT. K+1)) THEN
c             Mesh(I) = MeshNr
c            ENDIF
c          MeshNr = MeshNr + 1  
c 30       CONTINUE
c        MeshNR = MeshNr + 1
c 20     CONTINUE
c      MeshNr = MeshNr + 1
c 10   CONTINUE
 

      DO 40 I = 1, NALLg-1
       
         DO 50 J = I+1, NALLg
            
            D = DSQRT  ((Xhalf(I)- Xhalf(J))**2.d0 
     &                + (Yhalf(I)- Yhalf(J))**2.d0
     &                + (Zhalf(I)- Zhalf(J))**2.d0)
            DotV=VXg(I)*VXg(J)+VYg(I)*VYg(J)+VZg(I)*VZg(J)    !dot product of two velocity vectors
            IF (D .LE. Dcl .AND. DotV .LE. 0) THEN
               WX = VXg(I) - VXg(J)
               WY = VYg(I) - VYg(J)
               WZ = VZg(I) - VZg(J)
               VXtmp = VXg(I)   ! temporary velocities saved for computation
               VYtmp = VYg(I)   ! of velocities after encounter
               VZtmp = VZg(I)
               RX = Xhalf(J) - Xhalf(I)
               RY = Yhalf(J) - Yhalf(I)
               RZ = Zhalf(J) - Zhalf(I)
               R2 = RX**2.d0 + RY**2.d0 + RZ**2.d0
               VXg(I)=0.5d0*(BetaT+1)*WX - 0.5d0*(BetaT-BetaR)
     &                 *(RX*WX)*RX/R2 + VXg(J)
               VYg(I)=0.5d0*(BetaT+1)*WY - 0.5d0*(BetaT-BetaR)
     &                 *(RY*WY)*RY/R2 + VYg(J)
               VZg(I)=0.5d0*(BetaT+1)*WZ - 0.5d0*(BetaT-BetaR)
     &                 *(RZ*WZ)*RZ/R2 + VZg(J)
               VXg(J)=0.5d0*(BetaT+1)*WX + 0.5d0*(BetaT-BetaR)
     &                 *(RX*WX)*RX/R2 + VXtmp
               VYg(J)=0.5d0*(BetaT+1)*WY + 0.5d0*(BetaT-BetaR)
     &                 *(RY*WY)*RY/R2 + VYtmp
               VZg(J)=0.5d0*(BetaT+1)*WZ + 0.5d0*(BetaT-BetaR)
     &                 *(RZ*WZ)*RZ/R2 + VZtmp
               COL = COL + 1
               EXIT
            ENDIF
 50      CONTINUE
 40   CONTINUE
 
      OPEN(80,FILE='kolizie.dat') 
      WRITE(80,*) COL
      END
 
*************************************************************************
      SUBROUTINE Disrupt (D, V)
*************************************************************************
! test galaxy separation and if subcritical, then catastrophically 
!      disrupt the secondary (setting its mass to zero)
! relative velocity computed for information 
                                   
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'merge.inp'

      COMMON /GXIES/  Xp(3),Xs(3),Vp(3),Vs(3)
      COMMON /GXIES2/ GMp, Bp, Bp2, YBp2, CTrBp, Sigma0p, Vesc0p, Sr0p,
     &                GMs, Bs, Bs2, YBs2, CTrBs, Sigma0s, Vesc0s, Sr0s,
     &                GMeff, GMStoP, GMpart1,GMpart2, Eps2,Eps2s,Eps2p, 
     &                Omega23, Dini, Vini, Trmin, DTrBp, DTrBs, Eps,
     &                RtrBp, RtrBs, BgasP, BgasS, BgasP2, BgasS2,
     &                YBgasP2, YBgasS2, CTrBgasP, CtrBgasS, Vesc0gasP,
     &                Vesc0gasS, Sr0gasP, Sr0gasS, DTrBgasP, DTrBgasS,
     &                RTrBgasP, RTrBgasS, Eps2gasP, Eps2gasS   

      DATA IDisrupt, Rdisrupt /0, 0.1d0/
      SAVE IDisrupt, Rdisrupt
      
      IF (IDisrupt .EQ. 0) THEN  ! secondary still undisrupted, test separation
        
        D = Separation()                                  

        V2 = (Vs(1)-Vp(1))**2 + (Vs(2)-Vp(2))**2 + (Vs(3)-Vp(3))**2  ! v_rel^2
        V  = DSQRT (V2)     ! relative velocity

!       IF (D .LT. Rdisrupt * Bp) THEN  ! separation subcritical -> disruption
        IF (D .LT. Rdisrupt * Bp .OR. Xs(1). LT. 0.d0) THEN  
          GMs    = 0.d0                 ! set secondary's mass to 0        
          GMstoP = 0.

          DO 10 I = 1, 3
            Xp(I) = 0.d0;   Vp(I) = 0.d0  ! place primary to the center-of-mass
 10       CONTINUE
  
          IDisrupt = 1
          D   = 0.d0
          V   = 0.d0
        ENDIF

      ENDIF

      END

*************************************************************************
      REAL*8 FUNCTION SEPARATION()
*************************************************************************
! computes separation of the primary and secondary
                                   
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /GXIES/  Xp(3),Xs(3),Vp(3),Vs(3)

      D2 = (Xs(1)-Xp(1))**2 + (Xs(2)-Xp(2))**2 + (Xs(3)-Xp(3))**2 
      Separation = DSQRT (D2)     ! secondary-to-primary separation

      END
      

*************************************************************************
      REAL*8 FUNCTION ACCPLUMM (GM, X, R2, B2)
*************************************************************************
! computes acceleration for a Plummer potential
                                   
      IMPLICIT REAL*8 (A-H,O-Z)

      ACCPLUMM  = - GM * X / (R2 + B2)**1.5d0  ! acceleration

      END

*************************************************************************
      REAL*8 FUNCTION ACCHERNQ (GM, X, R, B)
*************************************************************************
! like ACCPLUMM but Primary potential = Hernquist
                                   
      IMPLICIT REAL*8 (A-H,O-Z)

      ACCHERNQ  = - GM * X  / (R*(R + B)**2.d0) ! acceleration

      END

*********************y***************************************************
      SUBROUTINE SetUnits
************************************************************************

      IMPLICIT REAL*8 (A-H,O-Z)

! G=1
! distance unit = 1 kpc
! time unit equal to the chosen time step
! DTMyr = time step in Myr
! velocity unit: VUNIT = 1 kpc/TUNIT = 977.88/DTMyr [km/s]
! mass unit: GMUNIT = 2.2233 * 10^11/(DTMyr^2) [MSun]
! VUNIT1,YVUNIT1 convert vel. from 100m/s to code unit & vice-versa during I/O

      COMMON /UNITS/ TUNIT,VUNIT,YVUNIT,GMUNIT,GMUNIT6

      TUNIT   = 0.01d0             ! time step in Myr = code time unit
!     TUNIT   = 1.0d0             ! time step in Myr = code time unit

      VUNIT   = 977.88D0 / TUNIT  ! velocity unit in km/s 
      YVUNIT  = 1./VUNIT

      GMUNIT  = 2.2233D+11 / (TUNIT**2)   ! mass unit in Mo
      GMUNIT6 = GMUNIT*(1.D-6) ! converts surface density to Mo/pc^2

      END


******************************************
      REAL*8 FUNCTION GASDEV( )
******************************************
* Gaussian random number generator using Box-Muller transformation
* Returns a normally distributed deviate with zero mean and unit variance
      
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*4 RAND
      DATA ISET /0/
      SAVE ISET,GSET

c      CALL SRAND(ISEED1)

      IF (ISET.EQ.0) THEN
1       V1=2.*RAND(0)-1. 
        V2=2.*RAND(0)-1.
	  R=V1**2+V2**2
	  IF (R.GE.1) GO TO 1
      FAC=DSQRT(-2.*LOG(R)/R)
	  GSET=V1*FAC
	  GASDEV=V2*FAC
	  ISET=1
        ELSE
          GASDEV=GSET
	  ISET=0
      ENDIF

      END		

******************************************************************
      SUBROUTINE PALETT(TYPE, CONTRA, BRIGHT)
******************************************************************
C Set a "palette" of colors in the range of color indices used by
C PGIMAG.
C-----------------------------------------------------------------------
      INTEGER TYPE
      REAL CONTRA, BRIGHT
C
      REAL GL(2), GR(2), GG(2), GB(2)
      REAL RL(9), RR(9), RG(9), RB(9)
      REAL HL(5), HR(5), HG(5), HB(5)
      REAL WL(10), WR(10), WG(10), WB(10)
      REAL AL(20), AR(20), AG(20), AB(20)
C
      DATA GL /0.0, 1.0/
      DATA GR /0.0, 1.0/
      DATA GG /0.0, 1.0/
      DATA GB /0.0, 1.0/
C
      DATA RL /-0.5, 0.0, 0.17, 0.33, 0.50, 0.67, 0.83, 1.0, 1.7/
      DATA RR / 0.0, 0.0,  0.0,  0.0,  0.6,  1.0,  1.0, 1.0, 1.0/
      DATA RG / 0.0, 0.0,  0.0,  1.0,  1.0,  1.0,  0.6, 0.0, 1.0/
      DATA RB / 0.0, 0.3,  0.8,  1.0,  0.3,  0.0,  0.0, 0.0, 1.0/
C
      DATA HL /0.0, 0.2, 0.4, 0.6, 1.0/
      DATA HR /0.0, 0.5, 1.0, 1.0, 1.0/
      DATA HG /0.0, 0.0, 0.5, 1.0, 1.0/
      DATA HB /0.0, 0.0, 0.0, 0.3, 1.0/
C
      DATA WL /0.0, 0.5, 0.5, 0.7, 0.7, 0.85, 0.85, 0.95, 0.95, 1.0/
      DATA WR /0.0, 1.0, 0.0, 0.0, 0.3,  0.8,  0.3,  1.0,  1.0, 1.0/
      DATA WG /0.0, 0.5, 0.4, 1.0, 0.0,  0.0,  0.2,  0.7,  1.0, 1.0/
      DATA WB /0.0, 0.0, 0.0, 0.0, 0.4,  1.0,  0.0,  0.0, 0.95, 1.0/
C
      DATA AL /0.0, 0.1, 0.1, 0.2, 0.2, 0.3, 0.3, 0.4, 0.4, 0.5,
     :         0.5, 0.6, 0.6, 0.7, 0.7, 0.8, 0.8, 0.9, 0.9, 1.0/
      DATA AR /0.0, 0.0, 0.3, 0.3, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0,
     :         0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0/
      DATA AG /0.0, 0.0, 0.3, 0.3, 0.0, 0.0, 0.0, 0.0, 0.8, 0.8,
     :         0.6, 0.6, 1.0, 1.0, 1.0, 1.0, 0.8, 0.8, 0.0, 0.0/
      DATA AB /0.0, 0.0, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.9, 0.9,
     :         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/
C
      IF (TYPE.EQ.1) THEN
C        -- gray scale
         CALL PGCTAB(GL, GR, GG, GB, 2, CONTRA, BRIGHT)
      ELSE IF (TYPE.EQ.2) THEN
C        -- rainbow
         CALL PGCTAB(RL, RR, RG, RB, 9, CONTRA, BRIGHT)
      ELSE IF (TYPE.EQ.3) THEN
C        -- heat
         CALL PGCTAB(HL, HR, HG, HB, 5, CONTRA, BRIGHT)
      ELSE IF (TYPE.EQ.4) THEN
C        -- weird IRAF
         CALL PGCTAB(WL, WR, WG, WB, 10, CONTRA, BRIGHT)
      ELSE IF (TYPE.EQ.5) THEN
C        -- AIPS
         CALL PGCTAB(AL, AR, AG, AB, 20, CONTRA, BRIGHT)
      END IF
      END

****************************************************************************
      SUBROUTINE FILENAME3 (I, IFILE, FileType, GifFile)
****************************************************************************
*Gives names and opens output files for 2D density and velocity maps
      CHARACTER IFILE*3
      CHARACTER FileNo*5, FileType*3, GifFile*9
      PARAMETER(M=48)

      IX    = I
      I4    = IX / 10000
      IX    = IX - I4*10000
      I3    = IX  / 1000
      IX    = IX - I3*1000
      I2    = IX / 100
      IX    = IX - I2*100
      I1    = IX / 10
      I0    = IX - I1*10    

      FileNo=CHAR(I4+M)//CHAR(I3+M)//CHAR(I2+M)//CHAR(I1+M)//CHAR(I0+M)

      IF (FileType .EQ. 'rho') THEN   ! open files for ascii data
        OPEN(20, FILE = IFILE//'s'//FileNo//'-rho')    ! stars - density
        OPEN(21, FILE = IFILE//'g'//FileNo//'-rho')    ! gas
!        OPEN(22, FILE = IFILE//'s'//FileNo//'-vrm')    ! stars - mean rad vel
!        OPEN(23, FILE = IFILE//'g'//FileNo//'-vrm')    ! gas
      ELSEIF (FileType .EQ. 'gif') THEN   ! assign name to a gif image
        GifFile = IFILE//'-'//FileNo         
      ENDIF

      RETURN
      END


*************************************************************************
      SUBROUTINE MPRODUCT (A, B, C, N)
*************************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)

      DIMENSION A(N,N), B(N,N), C(N,N)

      DO 10 I=1,N
        DO 10 J=1,N
          C(I,J) = 0.d0
          DO 20 K=1,N
            C(I,J) = C(I,J) + A(I,K) * B(K,J)
 20       CONTINUE  
 10   CONTINUE    

      END

******************************************
      REAL*8 FUNCTION YOUNG76b (S) 
******************************************
      IMPLICIT REAL*8 (A-H,O-Z)

      OPEN (12, FILE = 'young1976.dat')

      YOUNG76b = 0.d0    ! unnecessary initilization to avoid warning messages
                         !   during compilation 

      X1 = 0.;   Y1 = 0.  

      DO 10 I = 1, 124            ! no of lines in young1976.txt
        READ(12,*) X2, Y2
        IF (S .LE. X2) THEN       ! simpler condition than in Young76 
          YOUNG76b = Y1 + (Y2-Y1)/(X2-X1)*(S-X1) 
          GOTO 11
        ELSE 
          X1 = X2;  Y1 = Y2
        ENDIF
 10   CONTINUE
  
 11   CLOSE (12)

      END      

******************************************
      REAL*8 FUNCTION YOUNG76c (S) 
******************************************
* unlike Young76 & Young76b, doesn't read data from the table
*   but from arrays saved in the memory 

      IMPLICIT REAL*8 (A-H,O-Z)

      PARAMETER (NY = 124)
      COMMON /YOUNG/ YOUNG1(NY), YOUNG2(NY)
      SAVE /YOUNG/

      YOUNG76c = 0.    ! unnecessary initilization to avoid warning messages
                         !   during compilation 

      X1 = 0.;   Y1 = 0.  

      DO 10 I = 1, 124            ! no of lines in young1976.txt
        X2 = YOUNG1(I)
        Y2 = YOUNG2(I)  

        IF (S .LE. X2) THEN       ! simpler condition than in Young76
          YOUNG76c = Y1 + (Y2-Y1)/(X2-X1)*(S-X1) 
          GOTO 11
        ELSE 
          X1 = X2;  Y1 = Y2
        ENDIF
 10   CONTINUE
  
 11   RETURN

      END      


******************************************
      SUBROUTINE YOUNG76c0 
******************************************
* initializes arrays used by Young76c

      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NY = 124)
      COMMON /YOUNG/ YOUNG1(NY), YOUNG2(NY)
      SAVE /YOUNG/

      OPEN (12, FILE = 'young1976.dat')

      DO 10 I = 1, 124            ! no of lines in young1976.txt
        READ(12,*) YOUNG1(I), YOUNG2(I)
 10   CONTINUE
  
      CLOSE (12)

      END      
