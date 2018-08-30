define(REAL,`double precision')dnl
define(INTEGER,`integer')dnl
c     HELPER FUNCTIONS
c     -- dimension independent --

ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                   c
c     Interpolates cell average to cell center      c
c                                                   c
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine weno_interp_c_to_c_k(Q,I,k)
      implicit none
      INTEGER k
      REAL Q(-k:k)
      REAL I
      REAL c(0:k,0:(k-1))
      REAL wp(0:(k-1))
      REAL wn(0:(k-1))
      REAL s(0:1)

      INTEGER r,j
      REAL i_vals(0:(k-1))

      REAL si(0:(k-1))
      REAL eps
      REAL alpha(0:(k-1))
      REAL tot

      REAL wpf(0:(k-1))
      REAL wnf(0:(k-1))

      call calculate_weights(k,c,wp,wn,s)
      do r=0,k-1
        i_vals(r) = 0.d0
        do j=0,k-1
          i_vals(r) = i_vals(r)+c(r+1,j)*Q(j-r)
        enddo
      enddo
      select case(k)
        case(2)
          si(0) = (Q(1)-Q(0))**2
          si(1) = (Q(0)-Q(-1))**2
        case(3)
          si(0) = 13.d0/12.d0*(Q(0)-2.d0*Q(1)+Q(2))**2
     &          +0.25d0*(3.d0*Q(0)-4.d0*Q(1)+Q(2))**2
          si(1) = 13.d0/12.d0*(Q(-1)-2.d0*Q(0)+Q(1))**2
     &          +0.25d0*(Q(-1)-Q(1))**2
          si(2) = 13.d0/12.d0*(Q(-2)-2.d0*Q(-1)+Q(0))**2
     &          +0.25d0*(Q(-2)-4.d0*Q(-1)+3.d0*Q(0))**2
      end select

      tot = 0.d0
      eps = 1.0d-6
      do r=0,k-1
        alpha(r) = wp(r)/(eps+si(r))**2
        tot = tot+alpha(r)
      enddo
      tot = 1.d0/tot
      do r=0,k-1
        wpf(r) = tot*alpha(r)
      enddo
      do r=0,k-1
        alpha(r) = wn(r)/(eps+si(r))**2
        tot = tot+alpha(r)
      enddo
      tot = 1.d0/tot
      do r=0,k-1
        wnf(r) = tot*alpha(r)
      enddo
      I = 0.d0
      do r=0,k-1
        I = I + s(0)*wpf(r)*i_vals(r)
     &        - s(1)*wnf(r)*i_vals(r)
      enddo
      endsubroutine

ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                   c
c     Interpolates cell side to cell center         c
c                                                   c
ccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine weno_interp_s_to_c_k(Q,I,k)
      implicit none
      INTEGER k
      REAL Q((-k+1):k)
      REAL I

      REAL si(0:(k-1))
      REAL qi(0:(k-1))
      REAL wi(0:(k-1))
      REAL inv_tot
      REAL eps
      INTEGER j

      eps = 1.d0-6

      select case (k)
        case (2)
          qi(0) = -0.125d0*Q(-1)+0.75d0*Q(0)+0.325d0*Q(1)
          qi(1) = 0.375d0*Q(0)+0.75d0*Q(1)-0.125d0*Q(2)

          si(0) = 13.d0/12.d0*Q(-1)**2
     &       -13.d0/3.d0*Q(-1)*Q(0)
     &       +16.d0/3.d0*Q(0)**2
     &       +13.d0/6.d0*Q(-1)*Q(1)
     &       -19.d0/3.d0*Q(0)*Q(1)
     &       +25.d0/12.d0*Q(1)**2
          si(1) = 25.d0/12.d0*Q(0)**2
     &       -19.d0/3.d0*Q(0)*Q(1)
     &       +16.d0/3.d0*Q(1)**2
     &       +13.d0/6.d0*Q(0)*Q(2)
     &       -13.d0/3.d0*Q(1)*Q(2)
     &       +13.d0/12.d0*Q(2)**2

          wi(0) = 0.5d0/(eps+si(0))**2
          wi(1) = 0.5d0/(eps+si(1))**2
        case (3)
          qi(0) = 0.0625d0*Q(-2)-0.3125d0*Q(-1)
     &         +0.9375d0*Q(0)+0.3125d0*Q(1)
          qi(1) = -0.0625d0*Q(-1)+0.5625d0*Q(0)
     &         +0.5625d0*Q(1)-0.0625d0*Q(2)
          qi(2) = 0.3125d0*Q(0)+0.9375d0*Q(1)
     &         -0.3125d0*Q(2)+0.0625d0*Q(3)

          si(0) = 61.d0/45.d0*Q(-2)**2-553.d0/60.d0*Q(-2)*Q(-1)
     &      +248.d0/15.d0*Q(-1)**2+103.d0/10.d0*Q(-2)*Q(0)
     &      -2309.d0/60.d0*Q(-1)*Q(0)+721.d0/30.d0*Q(0)**2
     &      -683.d0/180.d0*Q(-2)*Q(1)+407.d0/90.d0*Q(1)**2
     &      +439.d0/30.d0*Q(-1)*Q(1)-1193.d0/60.d0*Q(0)*Q(1)
         si(1) = 61.d0/45.d0*Q(-1)**2-141.d0/20.d0*Q(-1)*Q(0)
     &      +331.d0/30.d0*Q(0)**2+179.d0/30.d0*Q(-1)*Q(1)
     &      -1259.d0/60.d0*Q(0)*Q(1)+331.d0/30.d0*Q(1)**2
     &      -293.d0/180.d0*Q(-1)*Q(2)+61.d0/45.d0*Q(2)**2
     &      +179.d0/30.d0*Q(0)*Q(2)-141.d0/20.d0*Q(1)*Q(2)
         si(2) = 407.d0/90.d0*Q(0)**2-1193.d0/60.d0*Q(0)*Q(1)
     &      +721.d0/30.d0*Q(1)**2+439.d0/30.d0*Q(0)*Q(2)
     &      -2309.d0/60.d0*Q(1)*Q(2)+248.d0/15.d0*Q(2)**2
     &      -683.d0/180.d0*Q(0)*Q(3)+61.d0/45.d0*Q(3)**2
     &      +103.d0/10.d0*Q(1)*Q(3)-553.d0/60.d0*Q(2)*Q(3)

         wi(0) = 0.1875d0/(eps+si(0))**2
         wi(1) = 0.625d0/(eps+si(1))**2
         wi(2) = 0.1875d0/(eps+si(2))**2
      endselect
      inv_tot = 0.d0
      do j=0,k-1
        inv_tot = inv_tot+wi(j)
      enddo
      I = 0.d0
      do j=0,k-1
        wi(j) = inv_tot*wi(j)
        I = I + wi(j)*qi(j)
      enddo
      end subroutine

ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                   c
c     Differentiates cell center quantities         c
c     to cell center                                c
c                                                   c
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine weno_der(q, dqdx, dx, k)
      implicit none
      INTEGER k
      REAL dqdx
      REAL dx
      REAL q(-k:k)

      REAL dq(0:k)
      REAL si(0:k)
      REAL w(0:k)
      REAL inv_tot
      REAL inv_dx
      REAL eps
      INTEGER j

      eps = 1.d-6
      select case (k)
        case (2)
          dq(0) = 0.5d0*(q(-1)-4.d0*q(-1)+3.d0*q(0))
          dq(1) = 0.5d0*(q(1)-q(-1))
          dq(2) = 0.5d0*(-3.d0*q(0)+4.d0*q(1)-q(2))

          si(0) = 25.d0/12.d0*q(-2)**2-31.d0/3.d0*q(-2)*q(-1)
     &          +40.d0/3.d0*q(-1)**2+37.d0/6.d0*q(-2)*q(0)
     &          -49.d0/3.d0*q(-1)*q(0)
     &          +61.d0/12.d0*q(0)**2
          si(1) = 13.d0/12.d0*q(-1)**2-13.d0/3.d0*q(-1)*q(0)
     &          +16.d0/3.d0*q(0)**2+13.d0/6.d0*q(-1)*q(1)
     &          -19.d0/3.d0*q(0)*q(1)+25.d0/12.d0*q(1)**2
          si(2) = 25.d0/12.d0*q(0)**2-19.d0/3.d0*q(0)*q(1)
     &          +16.d0/3.d0*q(1)**2+13.d0/6.d0*q(0)*q(2)
     &          -13.d0/3.d0*q(1)*q(2)+13.d0/12.d0*q(2)**2

          w(0) = 1.d0/(6.d0*(eps+si(0))**2)
          w(1) = 2.d0/(3.d0*(eps+si(1))**2)
          w(2) = 1.d0/(6.d0*(eps+si(2))**2)
        case(3)
          dq(0) = 1.d0/6.d0*(-2.d0*q(-3)+9.d0*q(-2)
     &          -18*q(-1)+11*q(0))
          dq(1) = 1.d0/6.d0*(q(-2)-6.d0*q(-1)
     &          +3.d0*q(0)+2.d0*q(1))
          dq(2) = 1.d0/6.d0*(-2.d0*q(-1)-3.d0*q(0)
     &          +6.d0*q(1)-q(2))
          dq(3) = 1.d0/6.d0*(-11.d0*q(0)+18.d0*q(1)
     &          -9.d0*q(2)+2.d0*q(3))

          si(0) = 407.d0/90.d0*q(-3)**2-1943.d0/60.d0*q(-3)*q(-2)
     &          +878.d0/15.d0*q(-2)**2+1189.d0/30.d0*q(-3)*q(-1)
     &          -8699.d0/60.d0*q(-2)*q(-1)+1373.d0/15.d0*q(-1)**2
     &          -2933.d0/180.d0*q(-3)*q(0)+603.d0/10.d0*q(-2)*q(0)
     &          -4663.d0/60.d0*q(-1)*q(0)+1517.d0/90.d0*q(0)**2
          si(1) = 61.d0/45.d0*q(-2)**2-553.d0/60.d0*q(-2)*q(-1)
     &          +248.d0/15.d0*q(-1)**2+103.d0/10.d0*q(-2)*q(0)
     &          -2309.d0/60.d0*q(-1)*q(0)+721.d0/30.d0*q(0)**2
     &          -683.d0/180.d0*q(-2)*q(1)+439.d0/30.d0*q(-1)*q(1)
     &          -1193.d0/60.d0*q(0)*q(1)+407.d0/90.d0*q(1)**2
          si(2) = 61.d0/45.d0*q(-1)**2-141.d0/20.d0*q(-1)*q(0)
     &          +331.d0/30.d0*q(0)**2+179.d0/30.d0*q(-1)*q(1)
     &          -1259.d0/60.d0*q(0)*q(1)+331.d0/30.d0*q(1)**2
     &          -293.d0/180.d0*q(-1)*q(2)+179.d0/30.d0*q(0)*q(2)
     &          -141.d0/20.d0*q(1)*q(2)+61.d0/45.d0*q(2)**2
          si(3) = 407.d0/90.d0*q(0)**2-1193.d0/60.d0*q(0)*q(1)
     &          +721.d0/30.d0*q(1)**2+439.d0/30.d0*q(0)*q(2)
     &          -2309.d0/60.d0*q(1)*q(2)+248.d0/15.d0*q(2)**2
     &          -683.d0/180.d0*q(0)*q(3)+103.d0/10.d0*q(1)*q(3)
     &          -553.d0/60.d0*q(2)*q(3)+61.d0/45.d0*q(3)**2

          w(0) = 1.d0/(20.d0*(1.0d-6+si(0))**2)
          w(1) = 9.d0/(20.d0*(1.0d-6+si(1))**2)
          w(2) = 9.d0/(20.d0*(1.0d-6+si(2))**2)
          w(3) = 1.d0/(20.d0*(1.0d-6+si(3))**2)
      endselect
      inv_tot = 0.d0
      do j=0,k
        inv_tot = inv_tot+w(j)
      enddo
      dqdx = 0.d0
      do j=0,k
        w(j) = inv_tot*w(j)
        dqdx = dqdx + w(j)*dq(j)
      enddo
      dqdx = dqdx/dx
      end subroutine

      subroutine calculate_weights(k, c, wp, wn, s)
      INTEGER k
      REAL c(0:k,0:k-1)
      REAL wp(0:k-1)
      REAL wn(0:k-1)
      REAL s(0:1)

      select case (k)
        case (3)
          c(0,0) = 2.958333333333325d0
          c(0,1) = -2.916666666666655d0
          c(0,2) = 0.958333333333329d0
          c(1,0) = 0.958333333333329d0
          c(1,1) = 0.083333333333333d0
          c(1,2) = -0.041666666666666d0
          c(2,0) = -0.041666666666d0
          c(2,1) = 1.083333333333337d0
          c(2,2) = -0.041666666666666d0
          c(3,0) = -0.041666666666666d0
          c(3,1) = 0.083333333333333d0
          c(3,2) = 0.958333333333337d0
          wp(0) = 9.d0/80.d0
          wp(2) = 9.d0/80.d0
          wp(1) = 49.d0/20.d0
          wn(0) = 9.d0/40.d0
          wn(1) = 49.d0/40.d0
          wn(2) = 9.d0/40.d0
          s(0) = 214.d0/80.d0
          s(1) = 67.d0/40.d0
      endselect
      endsubroutine
