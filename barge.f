! --- barge.f
!     Program for response analysis of 2d rectangular barge.
!
!     computation for added mass at infinite frequency has been added.
!     time dependence of exp(i omega t)
!
!     last upadated on 2023/03/29
!
      PROGRAM BARGE
      IMPLICIT NONE
      INTEGER N_MAX
      DOUBLE PRECISION PI
      DOUBLE COMPLEX ZI
      PARAMETER (N_MAX=30, ZI=(0.D0,1.D0))
      INTEGER nn , i , ipvts(0:N_MAX) , ipvta(0:N_MAX) , nf , j , info
      DOUBLE PRECISION rho , k , h , d , b , k_1 , b22_far , b33_far , 
     &                 b44_far , b24_far , b42_far , f , g , omega , 
     &                 df , zg , cc , cg , m , rx , gm , bb , kg , 
     &                 fdrift_fixed , fdrift_motion
      DOUBLE PRECISION k2(0:N_MAX) , n1(0:N_MAX) , n2(0:N_MAX) , 
     &                 wkn(N_MAX) , g2(0:N_MAX) , g3(0:N_MAX) , 
     &                 g4(0:N_MAX) , g7(0:N_MAX) , c(0:N_MAX,0:N_MAX) , 
     &                 fs(0:N_MAX,0:N_MAX) , fa(0:N_MAX,0:N_MAX) ,
     &                 err(5)
      DOUBLE PRECISION a22inf , a33inf , a44inf , a24inf , a42inf ,
     &                 k3(0:N_MAX) , n3(0:N_MAX) , 
     &                 rmats(0:N_MAX,0:N_MAX) , rmata(0:N_MAX,0:N_MAX) ,
     &                 rg1(0:N_MAX) , rg2(0:N_MAX) , rg3(0:N_MAX) ,
     &                 rg4(0:N_MAX) , rg5(0:N_MAX) , rg6(0:N_MAX) ,
     &                 rg7(0:N_MAX) , ra2(0:N_MAX) , ra3(0:N_MAX) ,
     &                 ra4(0:N_MAX) , rba(0:N_MAX,2) , rb2(0:N_MAX) ,
     &                 rb3(0:N_MAX) , rb4(0:N_MAX)
      DOUBLE COMPLEX r , t , a2p , a3p , a4p , a2m , a3m , a4m , a33 , 
     &               a44 , a22 , a24 , a42 , f2_far , f3_far , f4_far , 
     &               f2 , f3 , f4 , r_motion , t_motion
      DOUBLE COMPLEX k1(0:N_MAX) , mats(0:N_MAX,0:N_MAX) , 
     &               mata(0:N_MAX,0:N_MAX) , as(0:N_MAX) , aa(0:N_MAX) ,
     &               bs(0:N_MAX) , ba(0:N_MAX) , a2(0:N_MAX) , 
     &               a3(0:N_MAX) , a4(0:N_MAX) , b2(0:N_MAX) , 
     &               b3(0:N_MAX) , b4(0:N_MAX) , g1(0:N_MAX) , 
     &               g5(0:N_MAX) , g6(0:N_MAX) ,
     &               zbs(0:N_MAX,2) , zba(0:N_MAX,3) , 
     &               za(3,3) , zf(3)
      PI = 4.D0*atan(1.D0)
! --- start of input data ----------------------------------------------
      PRINT * , 'input fluid density, rho(Mg/m^3)'
      READ * , rho
      PRINT * , 'input gravity acceleration, g(m/s^2)'
      READ * , g
      PRINT * , 'input full width of the barge, B(m)'
      READ * , bb
      PRINT * , 'input draft of the barge, d(m)'
      READ * , d
      PRINT * , 'input water depth, h(m)'
      READ * , h
      PRINT * , 'input height of the center of gravity, KG(m)'
      READ * , kg
      PRINT * , 'input radius of gyration, rx(m)'
      READ * , rx
      PRINT * , 'input interval of frequencies, df(Hz)'
      READ * , df
      PRINT * , 'input number of frequencies, nf'
      READ * , nf
! --- end of input data ------------------------------------------------
      OPEN (9,FILE='barge.dat')
      WRITE (9,'(1p,e12.5,a)') rho , ' ! fluid density, rho(Mg/m^3)'
      WRITE (9,'(1p,e12.5,a)') g , ' ! gravity acceleration, g(m/s^2)'
      WRITE (9,'(1p,e12.5,a)') bb, ' ! full width of the barge, B(m)'
      WRITE (9,'(1p,e12.5,a)') d , ' ! draft of the barge, d(m)'
      WRITE (9,'(1p,e12.5,a)') h , ' ! water depth, h(m)'
      WRITE (9,'(1p,e12.5,a)') kg ,' ! center of gravity, KG(m)'
      WRITE (9,'(1p,e12.5,a)') rx ,' ! radius of gyration, rx(m)'
      WRITE (9,'(1p,e12.5,a)') df , 
     &                         ' ! interval of frequencies, df(Hz)'
      WRITE (9,'(i12,a)') nf , ' ! number of frequencies, nf'
      CLOSE (9)
!
      OPEN (10,FILE='amass.csv')
      OPEN (11,FILE='rdamp.csv')
      OPEN (12,FILE='rdamp_far.csv')
      OPEN (13,FILE='force.csv')
      OPEN (14,FILE='force_far.csv')
      OPEN (15,FILE='rao.csv')
      OPEN (16,FILE='rt.csv')
      OPEN (17,FILE='drift.csv')
      OPEN (18,FILE='wnumber.csv')
      OPEN (19,FILE='kochin.csv')
      OPEN (20,FILE='check.csv')
      OPEN (21,FILE='amass_inf.csv') ! added mass at infinite frequency
!
      nn = N_MAX
      b = bb/2
      zg = kg-d
!
      WRITE (10,99001)
      WRITE (11,99002)
      WRITE (12,99002)
      WRITE (13,99003)
      WRITE (14,99003)
      WRITE (15,99004)
      WRITE (16,99005)
      WRITE (17,99006)
      WRITE (18,99007)
      WRITE (19,99009)
      WRITE (20,99010)
      WRITE (21,99011)
!
      m = rho*bb*d
      gm = d/2 + bb**2/(12*d) - kg ! GM=KB+BM-KG
!
      DO i = 1 , nf
         f = df*i
         omega = 2*pi*f
         k = omega**2/g
         CALL COEF1(k,k_1,h,d,nn,k1,k2,n1,n2,wkn)
         CALL COEF2(h,d,nn,N_MAX,k1,k2,n1,n2,c)
         CALL CALFSA(h,b,nn,N_MAX,k2,c,fs,fa)
         CALL MATRX(h,b,nn,N_MAX,k1,fs,fa,mats,mata,as,aa)
         CALL CALGM(h,d,b,zg,nn,k1,k2,n1,n2,g1,g2,g3,g4,g5,g6,g7)
         CALL VECTR(h,b,nn,N_MAX,k2,c,g1,g2,g3,g4,g5,g6,a2,a3,a4)
         DO j = 0 , nn
            zbs(j,1) = as(j)
            zbs(j,2) = a3(j)
         ENDDO
         CALL ZGESV(nn+1,2,mats,N_MAX+1,ipvts,zbs,N_MAX+1,info)
         if (info .ne. 0) stop 'error at ZGESV 1'
         DO j = 0 , nn
            as(j) = zbs(j,1)
            a3(j) = zbs(j,2)
         ENDDO
         DO j = 0 , nn
            zba(j,1) = aa(j)
            zba(j,2) = a2(j)
            zba(j,3) = a4(j)
         ENDDO
         CALL ZGESV(nn+1,3,mata,N_MAX+1,ipvta,zba,N_MAX+1,info)
         if (info .ne. 0) stop 'error at ZGESV 2'
         DO j = 0 , nn
            aa(j) = zba(j,1)
            a2(j) = zba(j,2)
            a4(j) = zba(j,3)
         ENDDO
         CALL CALB(k1,b,nn,N_MAX,c,g2,g4,as,aa,a2,a3,a4,bs,ba,b2,b3,b4)
         CALL FORCE(rho,g,h,b,nn,k_1,n1,k2,n2,g1,g6,g7,aa,ba,bs,f2,f3,
     &              f4)
         CALL ADMAS(rho,h,d,b,nn,k2,n2,g1,g6,g7,a2,a4,b2,b3,b4,a22,a33,
     &              a44,a24,a42)
         r = (as(0)+aa(0))*EXP(k1(0)*b)
         t = (as(0)-aa(0))*EXP(k1(0)*b)
         err(1) = abs(r)**2 + abs(t)**2 -1.d0
         a2p = -k1(0)*h*SINH(k_1*h)*EXP(k1(0)*b)/n1(0)*a2(0)/omega
         a3p = -k1(0)*h*SINH(k_1*h)*EXP(k1(0)*b)/n1(0)*a3(0)/omega
         a4p = -k1(0)*h**2*SINH(k_1*h)*EXP(k1(0)*b)/n1(0)*a4(0)/omega
         cc = omega/k_1
         cg = 0.5D0*cc*(1.D0+2*k_1*h/sinh(2*k_1*h))
         f2_far = -2*rho*g*a2p*cg
         f3_far = -2*rho*g*a3p*cg
         f4_far = -2*rho*g*a4p*cg
         a2m = -a2p
         a3m =  a3p
         a4m = -a4p
         b22_far = rho*g*cg*(a2m*CONJG(a2m)+a2p*CONJG(a2p))
         b33_far = rho*g*cg*(a3m*CONJG(a3m)+a3p*CONJG(a3p))
         b44_far = rho*g*cg*(a4m*CONJG(a4m)+a4p*CONJG(a4p))
         b24_far = rho*g*cg*(a2m*CONJG(a4m)+a2p*CONJG(a4p))
         b42_far = rho*g*cg*(a4m*CONJG(a2m)+a4p*CONJG(a2p))
         err(2) = abs(a2p-r*CONJG(a2p)-t*CONJG(a2m))
         err(3) = abs(a3p-r*CONJG(a3p)-t*CONJG(a3m))
         err(4) = abs(a4p-r*CONJG(a4p)-t*CONJG(a4m))
         za(1,1) = -omega**2*(m+a22)
         za(1,2) = 0.0D0
         za(1,3) = -omega**2*a24
         za(2,1) = 0.0D0
         za(2,2) = -omega**2*(m+a33) + rho*g*bb
         za(2,3) = 0.0D0
         za(3,1) = -omega**2*a42
         za(3,2) = 0.0D0
         za(3,3) = -omega**2*(m*rx**2+a44) + m*g*gm
!
         WRITE (10,'(1p,e11.4,5(",",e11.4))') f , DBLE(a22) , DBLE(a33)
     &                           , DBLE(a44) , DBLE(a24) , DBLE(a42)
         WRITE (11,'(1p,e11.4,5(",",e11.4))') f , -omega*DIMAG(a22) ,
     &                           -omega*DIMAG(a33) , -omega*DIMAG(a44)
     &                         , -omega*DIMAG(a24) , -omega*DIMAG(a42)
         WRITE (12,'(1p,e11.4,5(",",e11.4))') f , b22_far , b33_far ,
     &                             b44_far , b24_far , b42_far
         WRITE (13,'(1p,e11.4,6(",",e11.4))') f , f2 , f3 , f4
         WRITE (14,'(1p,e11.4,6(",",e11.4))') f , f2_far , f3_far , 
     &                             f4_far
!
         zf(1) = f2
         zf(2) = f3
         zf(3) = f4
         CALL ZGESV(3,1,za,3,ipvts,zf,3,info)
         if (info .ne. 0) stop 'error at ZGESV 2'
         WRITE (15,'(1p,e11.4,6(",",e11.4))') f , zf(1) , zf(2) , zf(3)
!
         r_motion = r + ZI*omega*(zf(1)*a2p+zf(2)*a3p+zf(3)*a4p)
         t_motion = t + ZI*omega*(zf(1)*a2m+zf(2)*a3m+zf(3)*a4m)
         err(5) = abs(r_motion)**2 + abs(t_motion)**2 -1.d0
!
         WRITE (16,'(1p,e11.4,8(",",e11.4))') f , r , t , r_motion , 
     &                             t_motion
         WRITE (20,'(1p,e11.4,5(",",e11.4))') f , (err(j), j = 1 , 5)
!
         fdrift_fixed  = -0.5D0*rho*g*cg/cc*(1.D0+abs(r)**2-abs(t)**2)
         fdrift_motion = -0.5D0*rho*g*cg/cc*
     &                   (1.D0+abs(r_motion)**2-abs(t_motion)**2)
         WRITE (17,'(1p,e11.4,2(",",e11.4))') f , fdrift_fixed , 
     &                                        fdrift_motion
         WRITE (18,'(1p,e11.4,4(",",e11.4))') f , omega , 1.D0/f , 
     &                                        k*bb , k_1*bb
         WRITE (19,'(1p,e11.4,12(",",e11.4))') f , a2p , a2m , a3p , 
     &                                        a3m , a4p , a4m
      ENDDO

!     calculation of added mass at infinite frequency.
      CALL COEF3(h,d,nn,N_MAX,k2,k3,n2,n3,c)
      CALL CALFSA(h,b,nn,N_MAX,k2,c,fs,fa)
      CALL RMATRX(h,b,nn,N_MAX,k3,fs,fa,rmats,rmata)
      CALL RCALGM(h,d,b,zg,nn,k3,k2,n3,n2,rg1,rg2,rg3,rg4,rg5,rg6,rg7)
      CALL RVECTR(h,b,nn,N_MAX,k2,c,rg1,rg2,rg3,rg4,rg5,rg6,ra2,ra3,ra4)
      CALL DGESV(nn+1,1,rmats,N_MAX+1,ipvts,ra3,N_MAX+1,info)
      if (info .ne. 0) stop 'error at DGESV 1'
      DO j = 0 , nn
         rba(j,1) = ra2(j)
         rba(j,2) = ra4(j)
      ENDDO
      CALL DGESV(nn+1,2,rmata,N_MAX+1,ipvta,rba,N_MAX+1,info)
      if (info .ne. 0) stop 'error at DGESV 2'
      DO j = 0 , nn
         ra2(j) = rba(j,1)
         ra4(j) = rba(j,2)
      ENDDO
      CALL CALRB(b,nn,N_MAX,c,rg2,rg4,ra2,ra3,ra4,rb2,rb3,rb4)
      CALL RADMAS(rho,h,d,b,nn,k2,n2,rg1,rg6,rg7,ra2,ra4,rb2,rb3,rb4,
     &            a22inf,a33inf,a44inf,a24inf,a42inf)
      WRITE (21,'(1p,e11.4,4(",",e11.4))') a22inf , a33inf
     &                           , a44inf , a24inf , a42inf
!
99001 FORMAT('f(Hz),A22(Mg/m),A33(Mg/m),A44(Mgm),A24(Mg),A42(Mg)')
99002 FORMAT('f(Hz),B22(Mg/ms),B33(Mg/ms),B44(Mgm/s),B24(Mg/s),'
     &       'B42(Mg/s)')
99003 FORMAT('f(Hz),f2r(kN/m2),f2i(kN/m2),f3r(kN/m2),f3i(kN/m2),'
     &      ,'f4r(kN/m),f4i(kN/m)')
99004 FORMAT('f(Hz),H2r(-),H2i(-),H3r(-),H3i(-),H4r(rad/m),H4i(rad/m)')
99005 FORMAT('f(Hz),Rr_fixed,Ri_fixed,Tr_fixed,Ti_fixed,Rr_motion,'
     &      ,'Ri_motion,Tr_motion,Ti_motion')
99006 FORMAT('f(Hz),Fdy_fixed(kN/m3),Fdy_motion(kN/m3)')
99007 FORMAT('f(Hz),omega(rad/s),T(s),KB(-),kB(-)')
99009 FORMAT('f(Hz),A2+(-),A2-(-),A3+(-),A3-(-),A4+(-),A4-(-)')
99010 FORMAT('f(Hz),Err1(-),Err2(-),Err3(-),Err4(-),Err5(-)')
99011 FORMAT('A22(Mg/m),A33(Mg/m),A44(Mgm),A24(Mg),A42(Mg)')
      END
! --- subroutine calb
      SUBROUTINE CALB(K1,B,Nn,N_max,C,G2,G4,As,Aa,A2,A3,A4,Bs,Ba,B2,B3,
     &                B4)
      IMPLICIT NONE
      INTEGER Nn , N_max , n , l
      DOUBLE PRECISION B , C(0:N_max,0:Nn) , G2(0:Nn) , G4(0:Nn)
      DOUBLE COMPLEX K1(0:Nn) , A2(0:Nn) , A3(0:Nn) , A4(0:Nn) , 
     &               B2(0:Nn) , B3(0:Nn) , B4(0:Nn) , As(0:Nn) , 
     &               Aa(0:Nn) , Bs(0:Nn) , Ba(0:Nn)
      DO n = 0 , Nn
         Bs(n) = 0.5D0*EXP(K1(0)*B)*C(n,0)
         Ba(n) = 0.5D0*EXP(K1(0)*B)*C(n,0)
         B2(n) = 0.D0
         B3(n) = -G2(n)
         B4(n) = -G4(n)
         DO l = 0 , Nn
            Bs(n) = Bs(n) + As(l)*C(n,l)
            Ba(n) = Ba(n) + Aa(l)*C(n,l)
            B2(n) = B2(n) + A2(l)*C(n,l)
            B3(n) = B3(n) + A3(l)*C(n,l)
            B4(n) = B4(n) + A4(l)*C(n,l)
         ENDDO
      ENDDO
      END
! --- subroutine calrb
      SUBROUTINE CALRB(B,Nn,N_max,C,G2,G4,A2,A3,A4,B2,B3,B4)
      IMPLICIT NONE
      INTEGER Nn , N_max , n , l
      DOUBLE PRECISION B , C(0:N_max,0:Nn) , G2(0:Nn) , G4(0:Nn)
      DOUBLE PRECISION A2(0:Nn) , A3(0:Nn) , A4(0:Nn) , 
     &               B2(0:Nn) , B3(0:Nn) , B4(0:Nn)
      DO n = 0 , Nn
         B2(n) = 0.D0
         B3(n) = -G2(n)
         B4(n) = -G4(n)
         DO l = 0 , Nn
            B2(n) = B2(n) + A2(l)*C(n,l)
            B3(n) = B3(n) + A3(l)*C(n,l)
            B4(n) = B4(n) + A4(l)*C(n,l)
         ENDDO
      ENDDO
      END
! --- subroutine force
      SUBROUTINE FORCE(Rho,G,H,B,Nn,K_1,N1,K2,N2,G1,G6,G7,Aa,Ba,Bs,F2,
     &                 F3,F4)
      IMPLICIT NONE
      INTEGER Nn , n
      DOUBLE PRECISION Rho , G , H , B , K_1 , N1(0:Nn) , K2(0:Nn) , 
     &                 N2(0:Nn) , G7(0:Nn)
      DOUBLE COMPLEX ZI , G1(0:Nn) , G6(0:Nn) , Aa(0:Nn) , Ba(0:Nn) , 
     &               Bs(0:Nn) , F2 , F3 , F4
      PARAMETER (ZI=(0.D0,1.D0))
      F2 = 0.5D0*EXP(ZI*K_1*B)*G1(0) + Aa(0)*G1(0)
      F3 = B*Bs(0)/N2(0)
      F4 = 0.5D0*H**2*EXP(ZI*K_1*B)*G6(0) + H**2*Aa(0)*G6(0) + Ba(0)
     &     *B**2/(3.D0*N2(0))
      DO n = 1 , Nn
         F2 = F2 + Aa(n)*G1(n)
         F3 = F3 + Bs(n)*(-1)**n*TANH(K2(n)*B)/(K2(n)*N2(n))
         F4 = F4 + H**2*Aa(n)*G6(n) + (-1)**n*Ba(n)/N2(n)*G7(n)
      ENDDO
      F2 = -2*Rho*G*N1(0)*H/COSH(K_1*H)*F2
      F3 = 2*Rho*G*N1(0)/COSH(K_1*H)*F3
      F4 = 2*Rho*G*N1(0)/COSH(K_1*H)*F4
      END
! --- subroutine admas
      SUBROUTINE ADMAS(Rho,H,D,B,Nn,K2,N2,G1,G6,G7,A2,A4,B2,B3,B4,A22,
     &                 A33,A44,A24,A42)
      IMPLICIT NONE
      INTEGER Nn , n
      DOUBLE PRECISION Rho , H , D , B , K2(0:Nn) , N2(0:Nn) , G7(0:Nn)
      DOUBLE COMPLEX G1(0:Nn) , G6(0:Nn) , A2(0:Nn) , A4(0:Nn) , 
     &               B2(0:Nn) , B3(0:Nn) , B4(0:Nn) , A22 , A33 , A44 , 
     &               A24 , A42
      A22 = A2(0)*G1(0)
      A33 = H*B*B3(0)/N2(0)
      A44 = H**2*A4(0)*G6(0) + B**2*B4(0)/(3.D0*N2(0))
      A24 = A4(0)*G1(0)
      A42 = H**2*A2(0)*G6(0) + B**2*B2(0)/(3.D0*N2(0))
      DO n = 1 , Nn
         A22 = A22 + A2(n)*G1(n)
         A33 = A33 + H*(-1)**n*TANH(K2(n)*B)*B3(n)/(N2(n)*K2(n))
         A44 = A44 + H**2*A4(n)*G6(n) + (-1)**n*B4(n)/N2(n)*G7(n)
         A24 = A24 + A4(n)*G1(n)
         A42 = A42 + H**2*A2(n)*G6(n) + (-1)**n*B2(n)/N2(n)*G7(n)
      ENDDO
      A22 = -Rho*2*H**2*A22
      A33 = Rho*2*(A33+((H-D)**2*B-B**3/3.D0)/(2*(H-D)))
      A44 = Rho*2*H**2*(A44+B**2*(B*(H-D)**2-B**3/5.D0)/(6*H**2*(H-D)))
      A24 = -Rho*2*H**3*A24
      A42 = Rho*2*H*A42
      END
! --- subroutine radmas
      SUBROUTINE RADMAS(Rho,H,D,B,Nn,K2,N2,G1,G6,G7,A2,A4,B2,B3,B4,A22,
     &                 A33,A44,A24,A42)
      IMPLICIT NONE
      INTEGER Nn , n
      DOUBLE PRECISION Rho , H , D , B , K2(0:Nn) , N2(0:Nn) , G7(0:Nn)
      DOUBLE PRECISION G1(0:Nn) , G6(0:Nn) , A2(0:Nn) , A4(0:Nn) , 
     &               B2(0:Nn) , B3(0:Nn) , B4(0:Nn) , A22 , A33 , A44 , 
     &               A24 , A42
      A22 = A2(0)*G1(0)
      A33 = H*B*B3(0)/N2(0)
      A44 = H**2*A4(0)*G6(0) + B**2*B4(0)/(3.D0*N2(0))
      A24 = A4(0)*G1(0)
      A42 = H**2*A2(0)*G6(0) + B**2*B2(0)/(3.D0*N2(0))
      DO n = 1 , Nn
         A22 = A22 + A2(n)*G1(n)
         A33 = A33 + H*(-1)**n*TANH(K2(n)*B)*B3(n)/(N2(n)*K2(n))
         A44 = A44 + H**2*A4(n)*G6(n) + (-1)**n*B4(n)/N2(n)*G7(n)
         A24 = A24 + A4(n)*G1(n)
         A42 = A42 + H**2*A2(n)*G6(n) + (-1)**n*B2(n)/N2(n)*G7(n)
      ENDDO
      A22 = -Rho*2*H**2*A22
      A33 = Rho*2*(A33+((H-D)**2*B-B**3/3.D0)/(2*(H-D)))
      A44 = Rho*2*H**2*(A44+B**2*(B*(H-D)**2-B**3/5.D0)/(6*H**2*(H-D)))
      A24 = -Rho*2*H**3*A24
      A42 = Rho*2*H*A42
      END
! --- subroutine vectr
      SUBROUTINE VECTR(H,B,Nn,N_max,K2,C,G1,G2,G3,G4,G5,G6,A2,A3,A4)
      IMPLICIT NONE
      INTEGER Nn , N_max , l , m
      DOUBLE PRECISION H , B , K2(0:Nn) , C(0:N_max,0:Nn) , G2(0:Nn) , 
     &                 G3(0:Nn) , G4(0:Nn)
      DOUBLE COMPLEX A3(0:Nn) , zsum3 , A4(0:Nn) , zsum4 , G5(0:Nn) , 
     &               G6(0:Nn) , G1(0:Nn) , A2(0:Nn)
      DO m = 0 , Nn
         zsum3 = 0.D0
         zsum4 = G4(0)*H/B*C(0,m)
         DO l = 1 , Nn
            zsum3 = zsum3 + G2(l)*K2(l)*H*TANH(K2(l)*B)*C(l,m)
            zsum4 = zsum4 + G4(l)*K2(l)*H/TANH(K2(l)*B)*C(l,m)
         ENDDO
         A2(m) = -G1(m)
         A3(m) = zsum3 + G3(m)
         A4(m) = zsum4 - G5(m) + G6(m)
      ENDDO
      END
! --- subroutine vectr
      SUBROUTINE RVECTR(H,B,Nn,N_max,K2,C,G1,G2,G3,G4,G5,G6,A2,A3,A4)
      IMPLICIT NONE
      INTEGER Nn , N_max , l , m
      DOUBLE PRECISION H , B , K2(0:Nn) , C(0:N_max,0:Nn) , G2(0:Nn) , 
     &                 G3(0:Nn) , G4(0:Nn)
      DOUBLE PRECISION A3(0:Nn) , zsum3 , A4(0:Nn) , zsum4 , G5(0:Nn) ,
     &               G6(0:Nn) , G1(0:Nn) , A2(0:Nn)
      DO m = 0 , Nn
         zsum3 = 0.D0
         zsum4 = G4(0)*H/B*C(0,m)
         DO l = 1 , Nn
            zsum3 = zsum3 + G2(l)*K2(l)*H*TANH(K2(l)*B)*C(l,m)
            zsum4 = zsum4 + G4(l)*K2(l)*H/TANH(K2(l)*B)*C(l,m)
         ENDDO
         A2(m) = -G1(m)
         A3(m) = zsum3 + G3(m)
         A4(m) = zsum4 - G5(m) + G6(m)
      ENDDO
      END
! --- subroutine matrx
      SUBROUTINE MATRX(H,B,Nn,N_max,K1,Fs,Fa,Mats,Mata,As,Aa)
      IMPLICIT NONE
      INTEGER Nn , N_max
      DOUBLE PRECISION H , B , Fs(0:N_max,0:Nn) , Fa(0:N_max,0:Nn)
      DOUBLE COMPLEX K1(0:Nn) , Mats(0:N_max,0:Nn) , Mata(0:N_max,0:Nn)
     &               , As(0:Nn) , Aa(0:Nn)
      INTEGER m , l
      DOUBLE PRECISION DELT
      EXTERNAL DELT
      DO m = 0 , Nn
         As(m) = 0.5D0*EXP(K1(0)*B)*(DELT(m,0)*K1(m)*H-Fs(m,0))
         Aa(m) = 0.5D0*EXP(K1(0)*B)*(DELT(m,0)*K1(m)*H-Fa(m,0))
         DO l = 0 , Nn
            Mats(m,l) = DELT(m,l)*K1(m)*H + Fs(m,l)
            Mata(m,l) = DELT(m,l)*K1(m)*H + Fa(m,l)
         ENDDO
      ENDDO
      END
! --- subroutine rmatrx
      SUBROUTINE RMATRX(H,B,Nn,N_max,K1,Fs,Fa,Mats,Mata)
      IMPLICIT NONE
      INTEGER Nn , N_max
      DOUBLE PRECISION H , B , Fs(0:N_max,0:Nn) , Fa(0:N_max,0:Nn)
      DOUBLE PRECISION K1(0:Nn) , Mats(0:N_max,0:Nn) , 
     &                 Mata(0:N_max,0:Nn)
      INTEGER m , l
      DOUBLE PRECISION DELT
      EXTERNAL DELT
      DO m = 0 , Nn
         DO l = 0 , Nn
            Mats(m,l) = DELT(m,l)*K1(m)*H + Fs(m,l)
            Mata(m,l) = DELT(m,l)*K1(m)*H + Fa(m,l)
         ENDDO
      ENDDO
      END
! --- subroutine calgm
      SUBROUTINE CALGM(H,D,B,Zg,Nn,K1,K2,N1,N2,G1,G2,G3,G4,G5,G6,G7)
      IMPLICIT NONE
      INTEGER Nn , n
      DOUBLE PRECISION H , D , B , Zg , K2(0:Nn) , N1(0:Nn) , N2(0:Nn) ,
     &                 G2(0:Nn) , G3(0:Nn) , G4(0:Nn) , G7(0:Nn)
      DOUBLE COMPLEX K1(0:Nn) , G5(0:Nn) , G6(0:Nn) , G1(0:Nn)
      G2(0) = 1.D0/(2*H**2*(H-D)*N2(0))*((H-D)**3/3.D0-B**2*(H-D))
      G4(0) = B/(6*H**3*(H-D)*N2(0))*((H-D)**3-B**2*(H-D))
      G7(0) = 0.D0
      DO n = 1 , Nn
         G2(n) = 1.D0/(2*H**2*(H-D)*N2(n))
     &           *(1.D0/K2(n)*((H-D)**2-2.D0/K2(n)**2)*SIN(K2(n)*(H-D))
     &           +2*(H-D)/K2(n)**2*COS(K2(n)*(H-D))-B**2/K2(n)
     &           *SIN(K2(n)*(H-D)))
         G4(n) = B/(2*H**3*(H-D)*N2(n))
     &           *(1.D0/K2(n)*((H-D)**2-2.D0/K2(n)**2)*SIN(K2(n)*(H-D))
     &           +2*(H-D)/K2(n)**2*COS(K2(n)*(H-D))-B**2/3/K2(n)
     &           *SIN(K2(n)*(H-D)))
         G7(n) = B**2*(1.D0/TANH(K2(n)*B)/(K2(n)*B)-1.D0/(K2(n)*B)**2)
      ENDDO
      DO n = 0 , Nn
         G1(n) = (SIN(K1(n)*H)-SIN(K1(n)*(H-D)))/(N1(n)*K1(n)*H)
         G3(n) = B/(H*(H-D))/(N1(n)*K1(n))*SIN(K1(n)*(H-D))
         G5(n) = 1.D0/(2*H**2*(H-D)*N1(n))
     &           *(1/K1(n)*((H-D)**2-2/K1(n)**2)*SIN(K1(n)*(H-D))
     &           +2*(H-D)/K1(n)**2*COS(K1(n)*(H-D))-B**2/K1(n)
     &           *SIN(K1(n)*(H-D)))
         G6(n) = 1.D0/(H**2*N1(n))
     &           *(1.D0/K1(n)**2*COS(K1(n)*H)+D/K1(n)*SIN(K1(n)*(H-D))
     &           -1.D0/K1(n)**2*COS(K1(n)*(H-D))) - Zg/H*G1(n)
      ENDDO
      END
! --- subroutine rcalgm
      SUBROUTINE RCALGM(H,D,B,Zg,Nn,K1,K2,N1,N2,G1,G2,G3,G4,G5,G6,G7)
      IMPLICIT NONE
      INTEGER Nn , n
      DOUBLE PRECISION H , D , B , Zg , K2(0:Nn) , N1(0:Nn) , N2(0:Nn) ,
     &                 G2(0:Nn) , G3(0:Nn) , G4(0:Nn) , G7(0:Nn)
      DOUBLE PRECISION K1(0:Nn) , G5(0:Nn) , G6(0:Nn) , G1(0:Nn)
      G2(0) = 1.D0/(2*H**2*(H-D)*N2(0))*((H-D)**3/3.D0-B**2*(H-D))
      G4(0) = B/(6*H**3*(H-D)*N2(0))*((H-D)**3-B**2*(H-D))
      G7(0) = 0.D0
      DO n = 1 , Nn
         G2(n) = 1.D0/(2*H**2*(H-D)*N2(n))
     &           *(1.D0/K2(n)*((H-D)**2-2.D0/K2(n)**2)*SIN(K2(n)*(H-D))
     &           +2*(H-D)/K2(n)**2*COS(K2(n)*(H-D))-B**2/K2(n)
     &           *SIN(K2(n)*(H-D)))
         G4(n) = B/(2*H**3*(H-D)*N2(n))
     &           *(1.D0/K2(n)*((H-D)**2-2.D0/K2(n)**2)*SIN(K2(n)*(H-D))
     &           +2*(H-D)/K2(n)**2*COS(K2(n)*(H-D))-B**2/3/K2(n)
     &           *SIN(K2(n)*(H-D)))
         G7(n) = B**2*(1.D0/TANH(K2(n)*B)/(K2(n)*B)-1.D0/(K2(n)*B)**2)
      ENDDO
      DO n = 0 , Nn
         G1(n) = (SIN(K1(n)*H)-SIN(K1(n)*(H-D)))/(N1(n)*K1(n)*H)
         G3(n) = B/(H*(H-D))/(N1(n)*K1(n))*SIN(K1(n)*(H-D))
         G5(n) = 1.D0/(2*H**2*(H-D)*N1(n))
     &           *(1/K1(n)*((H-D)**2-2/K1(n)**2)*SIN(K1(n)*(H-D))
     &           +2*(H-D)/K1(n)**2*COS(K1(n)*(H-D))-B**2/K1(n)
     &           *SIN(K1(n)*(H-D)))
         G6(n) = 1.D0/(H**2*N1(n))
     &           *(1.D0/K1(n)**2*COS(K1(n)*H)+D/K1(n)*SIN(K1(n)*(H-D))
     &           -1.D0/K1(n)**2*COS(K1(n)*(H-D))) - Zg/H*G1(n)
      ENDDO
      END
! --- subroutine calfsa
      SUBROUTINE CALFSA(H,B,Nn,N_max,K2,C,Fs,Fa)
      IMPLICIT NONE
      INTEGER Nn , N_max
      DOUBLE PRECISION H , B , K2(0:Nn) , C(0:N_max,0:Nn) , 
     &                 Fs(0:N_max,0:Nn) , Fa(0:N_max,0:Nn) , zsum2 , 
     &                 zsum4
      INTEGER m , l , n
      DO m = 0 , Nn
         DO n = 0 , Nn
            zsum2 = 0.D0
            zsum4 = 0.D0
            DO l = 1 , Nn
               zsum2 = zsum2 + K2(l)*H*TANH(K2(l)*B)*C(l,m)*C(l,n)
               zsum4 = zsum4 + K2(l)*H/TANH(K2(l)*B)*C(l,m)*C(l,n)
            ENDDO
            Fs(m,n) = zsum2
            Fa(m,n) = zsum4 + H/B*C(0,m)*C(0,n)
         ENDDO
      ENDDO
      END
! --- subroutine coef2
      SUBROUTINE COEF2(H,D,Nn,N_max,K1,K2,N1,N2,C)
      IMPLICIT NONE
      INTEGER Nn , N_max
      DOUBLE PRECISION H , D , K2(0:Nn) , N1(0:Nn) , N2(0:Nn) , 
     &                 C(0:N_max,0:Nn)
      DOUBLE COMPLEX K1(0:Nn)
      INTEGER m , n
      DO m = 0 , Nn
         DO n = 0 , Nn
            C(m,n) = 1.D0/(2*H*N2(m)*N1(n))
     &               *(1.D0/(K2(m)-K1(n))*SIN((K2(m)-K1(n))*(H-D))
     &               +1.D0/(K2(m)+K1(n))*SIN((K2(m)+K1(n))*(H-D)))
         ENDDO
      ENDDO
      END
! --- subroutine coef3
      SUBROUTINE COEF3(H,D,Nn,N_max,K2,K3,N2,N3,C)
      IMPLICIT NONE
      INTEGER Nn , N_max
      DOUBLE PRECISION pi
      DOUBLE PRECISION H , D , K2(0:Nn) , N2(0:Nn) , K3(0:Nn) , 
     &                 N3(0:Nn) , C(0:N_max,0:Nn) , EPSL
      INTEGER m , n
      EXTERNAL EPSL
      pi = 4.D0*atan(1.D0)

      DO n = 0 , Nn
         K3(n) = (2*n+1)*pi/(2*H)
         N3(n) = 1.D0/SQRT(2.D0)
      ENDDO

      DO m = 0 , Nn
         DO n = 0 , Nn
            if(K2(m).eq.K3(n)) then
               C(m,n) = 1.D0/(2*H*N2(m)*N3(n))
     &            *((H-D)
     &            +1.D0/(K2(m)+K3(n))*SIN((K2(m)+K3(n))*(H-D)))
            else
               C(m,n) = 1.D0/(2*H*N2(m)*N3(n))
     &            *(1.D0/(K2(m)-K3(n))*SIN((K2(m)-K3(n))*(H-D))
     &            +1.D0/(K2(m)+K3(n))*SIN((K2(m)+K3(n))*(H-D)))
            endif
         ENDDO
      ENDDO
      END
! --- subroutine coef1
      SUBROUTINE COEF1(K,K_1,H,D,Nn,K1,K2,N1,N2,Wkn)
      IMPLICIT NONE
      DOUBLE PRECISION PI
      DOUBLE COMPLEX ZI
      PARAMETER (ZI=(0.D0,1.D0))
      INTEGER Nn , n
      DOUBLE PRECISION K , K_1 , H , D , wk0 , K2(0:Nn) , 
     &                 N1(0:Nn) , N2(0:Nn) , Wkn(Nn) , EPSL
      DOUBLE COMPLEX K1(0:Nn)
      EXTERNAL EPSL
      PI = 4.D0*atan(1.D0)
      CALL WNUMBER(K*H,wk0)
      K_1 = wk0/H
      K1(0) = ZI*K_1
      N1(0) = SQRT((1.D0+SINH(2*wk0)/(2*wk0))/2)
      CALL EVANESD(K*H,Wkn,Nn)
      DO n = 1 , Nn
         K1(n) = Wkn(n)/H
         N1(n) = SQRT((1.D0+SIN(2*Wkn(n))/(2*Wkn(n)))/2)
      ENDDO
      DO n = 0 , Nn
         K2(n) = n*pi/(H-D)
         N2(n) = SQRT((1.D0-D/H)/EPSL(n))
      ENDDO
      END
! --- function delt
      DOUBLE PRECISION FUNCTION DELT(N,M)
      IMPLICIT NONE
      INTEGER N , M
      IF ( N.EQ.M ) THEN
         DELT = 1.D0
      ELSE
         DELT = 0.D0
      ENDIF
      END
! --- function epsl
      DOUBLE PRECISION FUNCTION EPSL(M)
      IMPLICIT NONE
      INTEGER M
      IF ( M.EQ.0 ) THEN
         EPSL = 1.D0
      ELSEIF ( M.GT.0 ) THEN
         EPSL = 2.D0
      ELSE
         EPSL = 0.D0
         PRINT * , 'error at epsl: m=' , M
      ENDIF
      END
! --- subroutine elanesd
      SUBROUTINE EVANESD(X,Wkn,N)
!
!     calculate evanescent wave number.
!
!     input  ... x       (n.d.)  (=omega*omega*h/g=k0*h*tanh(k0*h))
!     output ... wkn(n)  (n.d.)  (=kn*h)
!
!     programmed based on A.O.R., vol.12, pp.14-18, 1990.
!
      IMPLICIT NONE
      INTEGER N , i
      DOUBLE PRECISION pi
      DOUBLE PRECISION X , Wkn(N) , u , d , EVANESDF
      pi = 4.D0*atan(1.D0)
      IF ( X.LT.0.D0 .OR. N.LT.1 ) STOP 'error at evanes'
! -------------------------
      u = 3.0D0*X/(7.0D0+3.0D0*X)
      d = 0.0159D0 + u*(0.1032D0+u*(4.3152D0-u*2.8768D0))
      d = EVANESDF(X,d,1)
      Wkn(1) = pi - d
      IF ( N.LT.2 ) RETURN
! -------------------------
      DO i = 2 , N
         d = d - pi*X/(X*X+pi*(i+1)*(pi*i-d))
         d = EVANESDF(X,d,i)
         Wkn(i) = pi*i - d
      ENDDO
      END
! --- function evalesdf
      DOUBLE PRECISION FUNCTION EVANESDF(X,D,N)
      IMPLICIT NONE
      DOUBLE PRECISION X , D
      INTEGER N
      DOUBLE PRECISION pi
      DOUBLE PRECISION y , f0 , f1 , f2
      pi = 4.D0*atan(1.D0)
      y = pi*N - D
      f0 = DATAN(X/y) - D
      f1 = X/(X*X+y*y) - 1.0D0
      f2 = 2.0D0*X*y/(X*X+y*y)**2
      EVANESDF = D - f0/f1*(1.0D0+0.5D0*f0*f2/(f1*f1))
      END
! --- subroutine wnumber
      SUBROUTINE WNUMBER(X,Wk0)
      IMPLICIT NONE
      DOUBLE PRECISION X , a , p
      INTEGER i
      DOUBLE PRECISION Wk0 , b(6) , c(9) , sum
      DATA b/0.000000122D0 , 0.073250017D0 , -0.009899981D0 , 
     &     0.002640863D0 , -0.000829239D0 , -0.000176411D0/
      DATA c/1.D0 , -0.33333372D0 , -0.01109668D0 , 0.01726435D0 , 
     &     0.01325580D0 , -0.00116594D0 , 0.00829006D0 , -0.01252603D0 , 
     &     0.00404923D0/
      SAVE b , c
      IF ( X.LT.0.0D0 ) STOP 'error at wnumber'
      IF ( X.LE.2.0D0 ) THEN
         a = 0.5D0*X
         p = 1.0D0
         sum = c(1)
         DO i = 2 , 9
            p = p*a
            sum = sum + c(i)*p
         ENDDO
         Wk0 = SQRT(X)/sum
      ELSE
         a = 0.5D0*X*EXP(4.0D0-2.0D0*X)
         p = 1.0D0
         sum = b(1)
         DO i = 2 , 6
            p = p*a
            sum = sum + b(i)*p
         ENDDO
         Wk0 = X + sum
      ENDIF
      END
