define(NDIM,2)dnl
define(REAL,`double precision')dnl
define(INTEGER,`integer')dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim2d.i)dnl
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       Computes r = grad u
c
c       where u is a vector valued face centered quantity
c       and r = is a cell centered quantity
c       du contains the gradient of the x component of u
c       dv contains the gradient of the y component of u
c       uses adaptive WENO stencils
c       interpolation order is given by k_i
c       differentiation order is given by k_d
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine weno_grad_s_to_c(u_data_0, u_data_1, u_gcw,
     &                            du, dv, d_gcw,
     &                            dx, alpha,
     &                            ilower0, iupper0, ilower1, iupper1,
     &                            k_i, k_d)

      implicit none

      INTEGER ilower0, iupper0
      INTEGER ilower1, iupper1
      INTEGER k_i, k_d

      REAL alpha
      REAL dx(0:1)

      INTEGER u_gcw
      REAL u_data_0(FACE2d0(ilower,iupper,u_gcw))
      REAL u_data_1(FACE2d1(ilower,iupper,u_gcw))

      INTEGER d_gcw
      REAL du(CELL2d(ilower,iupper,d_gcw),0:1)
      REAL dv(CELL2d(ilower,iupper,d_gcw),0:1)

      INTEGER u_i_gcw
      INTEGER st_w
      REAL u_val(CELL2d(ilower,iupper,(k_d+1)),0:1)

      REAL inv_dx, inv_dy
      REAL um, up, u0
      REAL u_i(0:k_i)
      REAL w_i(0:k_i)
      REAL s_i(0:k_i)
      REAL u_d(0:k_d+1)
      REAL w_d(0:k_d+1)
      REAL s_d(0:k_d+1)
      REAL w0,w1,w2
      REAL s0,s1,s2
      REAL inv_tot
      INTEGER i0, i1, j, k

      st_w = k_d+1
      inv_dx = 1.d0/dx(0)
      inv_dy = 1.d0/dx(1)
!     reconstruct to cell sides for u
      do i0 = ilower0-st_w, iupper0+st_w
        do i1 = ilower1-st_w, iupper1+st_w
        if(k_i == 2) then
          u_i(0) = -0.125d0*u_data_0(i0-1,i1)
     &       +0.75d0*u_data_0(i0,i1)
     &       +0.375d0*u_data_0(i0+1,i1)
          u_i(1) = 0.375d0*u_data_0(i0,i1)
     &       +0.75d0*u_data_0(i0+1,i1)
     &       -0.125d0*u_data_0(i0+2,i1)

          s_i(0) = 13.d0/12.d0*u_data_0(i0-1,i1)**2
     &       -13.d0/3.d0*u_data_0(i0-1,i1)*u_data_0(i0,i1)
     &       +16.d0/3.d0*u_data_0(i0,i1)**2
     &       +13.d0/6.d0*u_data_0(i0-1,i1)*u_data_0(i0+1,i1)
     &       -19.d0/3.d0*u_data_0(i0,i1)*u_data_0(i0+1,i1)
     &       +25.d0/12.d0*u_data_0(i0+1,i1)**2
          s_i(1) = 25.d0/12.d0*u_data_0(i0,i1)**2
     &       -19.d0/3.d0*u_data_0(i0,i1)*u_data_0(i0+1,i1)
     &       +16.d0/3.d0*u_data_0(i0+1,i1)**2
     &       +13.d0/6.d0*u_data_0(i0,i1)*u_data_0(i0+2,i1)
     &       -13.d0/3.d0*u_data_0(i0+1,i1)*u_data_0(i0+2,i1)
     &       +13.d0/12.d0*u_data_0(i0+2,i1)**2

          w_i(0) = 1.d0/(2.d0*(1.d-6+s_i(0))**2)
          w_i(1) = 1.d0/(2.d0*(1.d-6+s_i(1))**2)
          inv_tot = 1.d0/(w_i(0)+w_i(1))
          w_i(0) = w_i(0)*inv_tot
          w_i(1) = w_i(1)*inv_tot
          u_val(i0,i1,0) = w_i(0)*u_i(0)+w_i(1)*u_i(1)
        elseif (k_i == 3) then
          u_i(0) = 0.0625d0*u_data_0(i0-2,i1)
     &   -0.3125d0*u_data_0(i0-1,i1)
     &   +0.9375d0*u_data_0(i0,i1)
     &   +0.3125d0*u_data_0(i0+1,i1)
          u_i(1) = -0.0625d0*u_data_0(i0-1,i1)
     &   +0.5625d0*u_data_0(i0,i1)
     &   +0.5625d0*u_data_0(i0+1,i1)
     &   -0.0625d0*u_data_0(i0+2,i1)
          u_i(2) = 0.3125d0*u_data_0(i0,i1)
     &   +0.9375d0*u_data_0(i0+1,i1)
     &   -0.3125d0*u_data_0(i0+2,i1)
     &   +0.0625d0*u_data_0(i0+3,i1)

         s_i(0) = 61.d0/45.d0*u_data_0(i0-2,i1)**2
     &      -553.d0/60.d0*u_data_0(i0-2,i1)*u_data_0(i0-1,i1)
     &      +248.d0/15.d0*u_data_0(i0-1,i1)**2
     &      +103.d0/10.d0*u_data_0(i0-2,i1)*u_data_0(i0,i1)
     &      -2309.d0/60.d0*u_data_0(i0-1,i1)*u_data_0(i0,i1)
     &      +721.d0/30.d0*u_data_0(i0,i1)**2
     &      -683.d0/180.d0*u_data_0(i0-2,i1)*u_data_0(i0+1,i1)
     &      +407.d0/90.d0*u_data_0(i0+1,i1)**2
     &      +439.d0/30.d0*u_data_0(i0-1,i1)*u_data_0(i0+1,i1)
     &      -1193.d0/60.d0*u_data_0(i0,i1)*u_data_0(i0+1,i1)
         s_i(1) = 61.d0/45.d0*u_data_0(i0-1,i1)**2
     &      -141.d0/20.d0*u_data_0(i0-1,i1)*u_data_0(i0,i1)
     &      +331.d0/30.d0*u_data_0(i0,i1)**2
     &      +179.d0/30.d0*u_data_0(i0-1,i1)*u_data_0(i0+1,i1)
     &      -1259.d0/60.d0*u_data_0(i0,i1)*u_data_0(i0+1,i1)
     &      +331.d0/30.d0*u_data_0(i0+1,i1)**2
     &      -293.d0/180.d0*u_data_0(i0-1,i1)*u_data_0(i0+2,i1)
     &      +61.d0/45.d0*u_data_0(i0+2,i1)**2
     &      +179.d0/30.d0*u_data_0(i0,i1)*u_data_0(i0+2,i1)
     &      -141.d0/20.d0*u_data_0(i0+1,i1)*u_data_0(i0+2,i1)
         s_i(2) = 407.d0/90.d0*u_data_0(i0,i1)**2
     &      -1193.d0/60.d0*u_data_0(i0,i1)*u_data_0(i0+1,i1)
     &      +721.d0/30.d0*u_data_0(i0+1,i1)**2
     &      +439.d0/30.d0*u_data_0(i0,i1)*u_data_0(i0+2,i1)
     &      -2309.d0/60.d0*u_data_0(i0+1,i1)*u_data_0(i0+2,i1)
     &      +248.d0/15.d0*u_data_0(i0+2,i1)**2
     &      -683.d0/180.d0*u_data_0(i0,i1)*u_data_0(i0+3,i1)
     &      +61.d0/45.d0*u_data_0(i0+3,i1)**2
     &      +103.d0/10.d0*u_data_0(i0+1,i1)*u_data_0(i0+3,i1)
     &      -553.d0/60.d0*u_data_0(i0+2,i1)*u_data_0(i0+3,i1)

          w_i(0) = 0.1875d0/((1.d-6+s_i(0))**2)
          w_i(1) = 0.625d0/((1.d-6+s_i(1))**2)
          w_i(2) = 0.1875d0/((1.d-6+s_i(2))**2)

          inv_tot = 1.d0/(w_i(0)+w_i(1)+w_i(2))

          w_i(0) = w_i(0)*inv_tot
          w_i(1) = w_i(1)*inv_tot
          w_i(2) = w_i(2)*inv_tot
          u_val(i0,i1,0) = w_i(0)*u_i(0)+w_i(1)*u_i(1)+w_i(2)*u_i(2)
        endif
        enddo
      enddo
!     reconstruct to cell sides for v
      do i1 = ilower1-st_w, iupper1+st_w
        do i0 = ilower0-st_w, iupper0+st_w
        if(k_i == 2) then
          u_i(0) = -0.125d0*u_data_1(i1-1,i0)
     &       +0.75d0*u_data_1(i1,i0)
     &       +0.375d0*u_data_1(i1+1,i0)
          u_i(1) = 0.375d0*u_data_1(i1,i0)
     &       +0.75d0*u_data_1(i1+1,i0)
     &       -0.125d0*u_data_1(i1+2,i0)

          s_i(0) = 13.d0/12.d0*u_data_1(i1-1,i0)**2
     &       -13.d0/3.d0*u_data_1(i1-1,i0)*u_data_1(i1,i0)
     &       +16.d0/3.d0*u_data_1(i1,i0)**2
     &       +13.d0/6.d0*u_data_1(i1-1,i0)*u_data_1(i1+1,i0)
     &       -19.d0/3.d0*u_data_1(i1,i0)*u_data_1(i1+1,i0)
     &       +25.d0/12.d0*u_data_1(i1+1,i0)**2
          s_i(1) = 25.d0/12.d0*u_data_1(i1,i0)**2
     &       -19.d0/3.d0*u_data_1(i1,i0)*u_data_1(i1+1,i0)
     &       +16.d0/3.d0*u_data_1(i1+1,i0)**2
     &       +13.d0/6.d0*u_data_1(i1,i0)*u_data_1(i1+2,i0)
     &       -13.d0/3.d0*u_data_1(i1+1,i0)*u_data_1(i1+2,i0)
     &       +13.d0/12.d0*u_data_1(i1+2,i0)**2

          w_i(0) = 1.d0/(2.d0*(1.d-6+s_i(0))**2)
          w_i(1) = 1.d0/(2.d0*(1.d-6+s_i(1))**2)
          inv_tot = 1.d0/(w_i(0)+w_i(1))
          w_i(0) = w_i(0)*inv_tot
          w_i(1) = w_i(1)*inv_tot
          u_val(i0,i1,0) = w_i(0)*u_i(0)+w_i(1)*u_i(1)
        elseif (k_i == 3) then
          u_i(0) = 0.0625d0*u_data_1(i1-2,i0)
     &   -0.3125d0*u_data_1(i1-1,i0)
     &   +0.9375d0*u_data_1(i1,i0)
     &   +0.3125d0*u_data_1(i1+1,i0)
          u_i(1) = -0.0625d0*u_data_1(i1-1,i0)
     &   +0.5625d0*u_data_1(i1,i0)
     &   +0.5625d0*u_data_1(i1+1,i0)
     &   -0.0625d0*u_data_1(i1+2,i0)
          u_i(2) = 0.3125d0*u_data_1(i1,i0)
     &   +0.9375d0*u_data_1(i1+1,i0)
     &   -0.3125d0*u_data_1(i1+2,i0)
     &   +0.0625d0*u_data_1(i1+3,i0)

         s_i(0) = 61.d0/45.d0*u_data_1(i1-2,i0)**2
     &      -553.d0/60.d0*u_data_1(i1-2,i0)*u_data_1(i1-1,i0)
     &      +248.d0/15.d0*u_data_1(i1-1,i0)**2
     &      +103.d0/10.d0*u_data_1(i1-2,i0)*u_data_1(i1,i0)
     &      -2309.d0/60.d0*u_data_1(i1-1,i0)*u_data_1(i1,i0)
     &      +721.d0/30.d0*u_data_1(i1,i0)**2
     &      -683.d0/180.d0*u_data_1(i1-2,i0)*u_data_1(i1+1,i0)
     &      +407.d0/90.d0*u_data_1(i1+1,i0)**2
     &      +439.d0/30.d0*u_data_1(i1-1,i0)*u_data_1(i1+1,i0)
     &      -1193.d0/60.d0*u_data_1(i1,i0)*u_data_1(i1+1,i0)
         s_i(1) = 61.d0/45.d0*u_data_1(i1-1,i0)**2
     &      -141.d0/20.d0*u_data_1(i1-1,i0)*u_data_1(i1,i0)
     &      +331.d0/30.d0*u_data_1(i1,i0)**2
     &      +179.d0/30.d0*u_data_1(i1-1,i0)*u_data_1(i1+1,i0)
     &      -1259.d0/60.d0*u_data_1(i1,i0)*u_data_1(i1+1,i0)
     &      +331.d0/30.d0*u_data_1(i1+1,i0)**2
     &      -293.d0/180.d0*u_data_1(i1-1,i0)*u_data_1(i1+2,i0)
     &      +61.d0/45.d0*u_data_1(i1+2,i0)**2
     &      +179.d0/30.d0*u_data_1(i1,i0)*u_data_1(i1+2,i0)
     &      -141.d0/20.d0*u_data_1(i1+1,i0)*u_data_1(i1+2,i0)
         s_i(2) = 407.d0/90.d0*u_data_1(i1,i0)**2
     &      -1193.d0/60.d0*u_data_1(i1,i0)*u_data_1(i1+1,i0)
     &      +721.d0/30.d0*u_data_1(i1+1,i0)**2
     &      +439.d0/30.d0*u_data_1(i1,i0)*u_data_1(i1+2,i0)
     &      -2309.d0/60.d0*u_data_1(i1+1,i0)*u_data_1(i1+2,i0)
     &      +248.d0/15.d0*u_data_1(i1+2,i0)**2
     &      -683.d0/180.d0*u_data_1(i1,i0)*u_data_1(i1+3,i0)
     &      +61.d0/45.d0*u_data_1(i1+3,i0)**2
     &      +103.d0/10.d0*u_data_1(i1+1,i0)*u_data_1(i1+3,i0)
     &      -553.d0/60.d0*u_data_1(i1+2,i0)*u_data_1(i1+3,i0)

          w_i(0) = 0.1875d0/((1.d-6+s_i(0))**2)
          w_i(1) = 0.625d0/((1.d-6+s_i(1))**2)
          w_i(2) = 0.1875d0/((1.d-6+s_i(2))**2)

          inv_tot = 1.d0/(w_i(0)+w_i(1)+w_i(2))

          w_i(0) = w_i(0)*inv_tot
          w_i(1) = w_i(1)*inv_tot
          w_i(2) = w_i(2)*inv_tot
          u_val(i0,i1,1) = w_i(0)*u_i(0)+w_i(1)*u_i(1)+w_i(2)*u_i(2)
        endif
        enddo
      enddo

      do j=0,1
        do i0 = ilower0, iupper0
          do i1 = ilower1, iupper1
          if(k_d == 2) then
            u_d(0) = 0.5d0*inv_dx*(u_val(i0-2,i1,j)
     &            -4.d0*u_val(i0-1,i1,j)+3.d0*u_val(i0,i1,j))
            u_d(1) = 0.5d0*inv_dx*(u_val(i0+1,i1,j)-u_val(i0-1,i1,j))
            u_d(2) = 0.5d0*inv_dx*(-3.d0*u_val(i0,i1,j)
     &            +4.d0*u_val(i0+1,i1,j)-u_val(i0+2,i1,j))

            s_d(0) = 25.d0/12.d0*u_val(i0-2,i1,j)**2
     &             -31.d0/3.d0*u_val(i0-2,i1,j)*u_val(i0-1,i1,j)
     &             +40.d0/3.d0*u_val(i0-1,i1,j)**2
     &             +37.d0/6.d0*u_val(i0-2,i1,j)*u_val(i0,i1,j)
     &             -49.d0/3.d0*u_val(i0-1,i1,j)*u_val(i0,i1,j)
     &             +61.d0/12.d0*u_val(i0,i1,j)**2
            s_d(1) = 13.d0/12.d0*u_val(i0-1,i1,j)**2
     &             -13.d0/3.d0*u_val(i0-1,i1,j)*u_val(i0,i1,j)
     &             +16.d0/3.d0*u_val(i0,i1,j)**2
     &             +13.d0/6.d0*u_val(i0-1,i1,j)*u_val(i0+1,i1,j)
     &             -19.d0/3.d0*u_val(i0,i1,j)*u_val(i0+1,i1,j)
     &             +25.d0/12.d0*u_val(i0+1,i1,j)**2
            s_d(2) = 25.d0/12.d0*u_val(i0,i1,j)**2
     &             -19.d0/3.d0*u_val(i0,i1,j)*u_val(i0+1,i1,j)
     &             +16.d0/3.d0*u_val(i0+1,i1,j)**2
     &             +13.d0/6.d0*u_val(i0,i1,j)*u_val(i0+2,i1,j)
     &             -13.d0/3.d0*u_val(i0+1,i1,j)*u_val(i0+2,i1,j)
     &             +13.d0/12.d0*u_val(i0+2,i1,j)**2

            w_d(0) = 1.d0/(6.d0*(1.d-6+s_d(0))**2)
            w_d(1) = 2.d0/(3.d0*(1.d-6+s_d(1))**2)
            w_d(2) = 1.d0/(6.d0*(1.d-6+s_d(2))**2)
          elseif(k_d == 3) then
            u_d(0) = inv_dx/6.d0*(-2.d0*u_val(i0-3,i1,j)
     &             +9.d0*u_val(i0-2,i1,j)-18*u_val(i0-1,i1,j)
     &             +11*u_val(i0,i1,j))
            u_d(1) = inv_dx/6.d0*(u_val(i0-2,i1,j)
     &             -6.d0*u_val(i0-1,i1,j)+3.d0*u_val(i0,i1,j)
     &             +2.d0*u_val(i0+1,i1,j))
            u_d(2) = inv_dx/6.d0*(-2.d0*u_val(i0-1,i1,j)
     &             -3.d0*u_val(i0,i1,j)+6.d0*u_val(i0+1,i1,j)
     &             -u_val(i0+2,i1,j))
            u_d(3) = inv_dx/6.d0*(-11.d0*u_val(i0,i1,j)
     &             +18.d0*u_val(i0+1,i1,j)-9.d0*u_val(i0+2,i1,j)
     &             +2.d0*u_val(i0+3,i1,j))

            s_d(0) = 407.d0/90.d0*u_val(i0-3,i1,j)**2
     &             -1943.d0/60.d0*u_val(i0-3,i1,j)*u_val(i0-2,i1,j)
     &             +878.d0/15.d0*u_val(i0-2,i1,j)**2
     &             +1189.d0/30.d0*u_val(i0-3,i1,j)*u_val(i0-1,i1,j)
     &             -8699.d0/60.d0*u_val(i0-2,i1,j)*u_val(i0-1,i1,j)
     &             +1373.d0/15.d0*u_val(i0-1,i1,j)**2
     &             -2933.d0/180.d0*u_val(i0-3,i1,j)*u_val(i0,i1,j)
     &             +603.d0/10.d0*u_val(i0-2,i1,j)*u_val(i0,i1,j)
     &             -4663.d0/60.d0*u_val(i0-1,i1,j)*u_val(i0,i1,j)
     &             +1517.d0/90.d0*u_val(i0,i1,j)**2
            s_d(1) = 61.d0/45.d0*u_val(i0-2,i1,j)**2
     &             -553.d0/60.d0*u_val(i0-2,i1,j)*u_val(i0-1,i1,j)
     &             +248.d0/15.d0*u_val(i0-1,i1,j)**2
     &             +103.d0/10.d0*u_val(i0-2,i1,j)*u_val(i0,i1,j)
     &             -2309.d0/60.d0*u_val(i0-1,i1,j)*u_val(i0,i1,j)
     &             +721.d0/30.d0*u_val(i0,i1,j)**2
     &             -683.d0/180.d0*u_val(i0-2,i1,j)*u_val(i0+1,i1,j)
     &             +439.d0/30.d0*u_val(i0-1,i1,j)*u_val(i0+1,i1,j)
     &             -1193.d0/60.d0*u_val(i0,i1,j)*u_val(i0+1,i1,j)
     &             +407.d0/90.d0*u_val(i0+1,i1,j)**2
            s_d(2) = 61.d0/45.d0*u_val(i0-1,i1,j)**2
     &             -141.d0/20.d0*u_val(i0-1,i1,j)*u_val(i0,i1,j)
     &             +331.d0/30.d0*u_val(i0,i1,j)**2
     &             +179.d0/30.d0*u_val(i0-1,i1,j)*u_val(i0+1,i1,j)
     &             -1259.d0/60.d0*u_val(i0,i1,j)*u_val(i0+1,i1,j)
     &             +331.d0/30.d0*u_val(i0+1,i1,j)**2
     &             -293.d0/180.d0*u_val(i0-1,i1,j)*u_val(i0+2,i1,j)
     &             +179.d0/30.d0*u_val(i0,i1,j)*u_val(i0+2,i1,j)
     &             -141.d0/20.d0*u_val(i0+1,i1,j)*u_val(i0+2,i1,j)
     &             +61.d0/45.d0*u_val(i0+2,i1,j)**2
            s_d(3) = 407.d0/90.d0*u_val(i0,i1,j)**2
     &             -1193.d0/60.d0*u_val(i0,i1,j)*u_val(i0+1,i1,j)
     &             +721.d0/30.d0*u_val(i0+1,i1,j)**2
     &             +439.d0/30.d0*u_val(i0,i1,j)*u_val(i0+2,i1,j)
     &             -2309.d0/60.d0*u_val(i0+1,i1,j)*u_val(i0+2,i1,j)
     &             +248.d0/15.d0*u_val(i0+2,i1,j)**2
     &             -683.d0/180.d0*u_val(i0,i1,j)*u_val(i0+3,i1,j)
     &             +103.d0/10.d0*u_val(i0+1,i1,j)*u_val(i0+3,i1,j)
     &             -553.d0/60.d0*u_val(i0+2,i1,j)*u_val(i0+3,i1,j)
     &             +61.d0/45.d0*u_val(i0+3,i1,j)**2

            w_d(0) = 1.d0/(20.d0*(1.0d-6+s_d(0))**2)
            w_d(1) = 9.d0/(20.d0*(1.0d-6+s_d(1))**2)
            w_d(2) = 9.d0/(20.d0*(1.0d-6+s_d(2))**2)
            w_d(3) = 1.d0/(20.d0*(1.0d-6+s_d(3))**2)
          endif

            inv_tot = 0.d0
            do k=0,k_d
              inv_tot = inv_tot + w_d(k)
            enddo
            inv_tot = 1.d0/inv_tot

            if(j == 0) then
              du(i0,i1,0) = 0.d0
              do k=0,k_d
                w_d(k) = inv_tot*w_d(k)
                du(i0,i1,0) = du(i0,i1,0) + w_d(k)*u_d(k)
              enddo
            else
              dv(i0,i1,0) = 0.d0
              do k=0,k_d
                w_d(k) = inv_tot*w_d(k)
                dv(i0,i1,0) = dv(i0,i1,0) + w_d(k)*u_d(k)
              enddo
            endif

          if(k_d == 2) then
            u_d(0) = 0.5d0*inv_dx*(u_val(i0,i1-2,j)
     &            -4.d0*u_val(i0,i1-1,j)+3.d0*u_val(i0,i1,j))
            u_d(1) = 0.5d0*inv_dx*(u_val(i0,i1+1,j)-u_val(i0,i1-1,j))
            u_d(2) = 0.5d0*inv_dx*(-3.d0*u_val(i0,i1,j)
     &            +4.d0*u_val(i0,i1+1,j)-u_val(i0,i1+2,j))

            s_d(0) = 25.d0/12.d0*u_val(i0,i1-2,j)**2
     &             -31.d0/3.d0*u_val(i0,i1-2,j)*u_val(i0,i1-1,j)
     &             +40.d0/3.d0*u_val(i0,i1-1,j)**2
     &             +37.d0/6.d0*u_val(i0,i1-2,j)*u_val(i0,i1,j)
     &             -49.d0/3.d0*u_val(i0,i1-1,j)*u_val(i0,i1,j)
     &             +61.d0/12.d0*u_val(i0,i1,j)**2
            s_d(1) = 13.d0/12.d0*u_val(i0,i1-1,j)**2
     &             -13.d0/3.d0*u_val(i0,i1-1,j)*u_val(i0,i1,j)
     &             +16.d0/3.d0*u_val(i0,i1,j)**2
     &             +13.d0/6.d0*u_val(i0,i1-1,j)*u_val(i0,i1+1,j)
     &             -19.d0/3.d0*u_val(i0,i1,j)*u_val(i0,i1+1,j)
     &             +25.d0/12.d0*u_val(i0,i1+1,j)**2
            s_d(2) = 25.d0/12.d0*u_val(i0,i1,j)**2
     &             -19.d0/3.d0*u_val(i0,i1,j)*u_val(i0,i1+1,j)
     &             +16.d0/3.d0*u_val(i0,i1+1,j)**2
     &             +13.d0/6.d0*u_val(i0,i1,j)*u_val(i0,i1+2,j)
     &             -13.d0/3.d0*u_val(i0,i1+1,j)*u_val(i0,i1+2,j)
     &             +13.d0/12.d0*u_val(i0,i1+2,j)**2

            w_d(0) = 1.d0/(6.d0*(1.d-6+s_d(0))**2)
            w_d(1) = 2.d0/(3.d0*(1.d-6+s_d(1))**2)
            w_d(2) = 1.d0/(6.d0*(1.d-6+s_d(2))**2)
          elseif(k_d == 3) then
            u_d(0) = inv_dx/6.d0*(-2.d0*u_val(i0,i1-3,j)
     &             +9.d0*u_val(i0,i1-2,j)-18.d0*u_val(i0,i1-1,j)
     &             +11.d0*u_val(i0,i1,j))
            u_d(1) = inv_dx/6.d0*(u_val(i0,i1-2,j)
     &             -6.d0*u_val(i0,i1-1,j)+3.d0*u_val(i0,i1,j)
     &             +2.d0*u_val(i0,i1+1,j))
            u_d(2) = inv_dx/6.d0*(-2.d0*u_val(i0,i1-1,j)
     &             -3.d0*u_val(i0,i1,j)+6.d0*u_val(i0,i1+1,j)
     &             -u_val(i0,i1+2,j))
            u_d(3) = inv_dx/6.d0*(-11.d0*u_val(i0,i1,j)
     &             +18.d0*u_val(i0,i1+1,j)-9.d0*u_val(i0,i1+2,j)
     &             +2.d0*u_val(i0,i1+3,j))

            s_d(0) = 407.d0/90.d0*u_val(i0,i1-3,j)**2
     &             -1943.d0/60.d0*u_val(i0,i1-3,j)*u_val(i0,i1-2,j)
     &             +878.d0/15.d0*u_val(i0,i1-2,j)**2
     &             +1189.d0/30.d0*u_val(i0,i1-3,j)*u_val(i0,i1-1,j)
     &             -8699.d0/60.d0*u_val(i0,i1-2,j)*u_val(i0,i1-1,j)
     &             +1373.d0/15.d0*u_val(i0,i1-1,j)**2
     &             -2933.d0/180.d0*u_val(i0,i1-3,j)*u_val(i0,i1,j)
     &             +603.d0/10.d0*u_val(i0,i1-2,j)*u_val(i0,i1,j)
     &             -4663.d0/60.d0*u_val(i0,i1-1,j)*u_val(i0,i1,j)
     &             +1517.d0/90.d0*u_val(i0,i1,j)**2
            s_d(1) = 61.d0/45.d0*u_val(i0,i1-2,j)**2
     &             -553.d0/60.d0*u_val(i0,i1-2,j)*u_val(i0,i1-1,j)
     &             +248.d0/15.d0*u_val(i0,i1-1,j)**2
     &             +103.d0/10.d0*u_val(i0,i1-2,j)*u_val(i0,i1,j)
     &             -2309.d0/60.d0*u_val(i0,i1-1,j)*u_val(i0,i1,j)
     &             +721.d0/30.d0*u_val(i0,i1,j)**2
     &             -683.d0/180.d0*u_val(i0,i1-2,j)*u_val(i0,i1+1,j)
     &             +439.d0/30.d0*u_val(i0,i1-1,j)*u_val(i0,i1+1,j)
     &             -1193.d0/60.d0*u_val(i0,i1,j)*u_val(i0,i1+1,j)
     &             +407.d0/90.d0*u_val(i0,i1+1,j)**2
            s_d(2) = 61.d0/45.d0*u_val(i0,i1-1,j)**2
     &             -141.d0/20.d0*u_val(i0,i1-1,j)*u_val(i0,i1,j)
     &             +331.d0/30.d0*u_val(i0,i1,j)**2
     &             +179.d0/30.d0*u_val(i0,i1-1,j)*u_val(i0,i1+1,j)
     &             -1259.d0/60.d0*u_val(i0,i1,j)*u_val(i0,i1+1,j)
     &             +331.d0/30.d0*u_val(i0,i1+1,j)**2
     &             -293.d0/180.d0*u_val(i0,i1-1,j)*u_val(i0,i1+2,j)
     &             +179.d0/30.d0*u_val(i0,i1,j)*u_val(i0,i1+2,j)
     &             -141.d0/20.d0*u_val(i0,i1+1,j)*u_val(i0,i1+2,j)
     &             +61.d0/45.d0*u_val(i0,i1+2,j)**2
            s_d(3) = 407.d0/90.d0*u_val(i0,i1,j)**2
     &             -1193.d0/60.d0*u_val(i0,i1,j)*u_val(i0,i1+1,j)
     &             +721.d0/30.d0*u_val(i0,i1+1,j)**2
     &             +439.d0/30.d0*u_val(i0,i1,j)*u_val(i0,i1+2,j)
     &             -2309.d0/60.d0*u_val(i0,i1+1,j)*u_val(i0,i1+2,j)
     &             +248.d0/15.d0*u_val(i0,i1+2,j)**2
     &             -683.d0/180.d0*u_val(i0,i1,j)*u_val(i0,i1+3,j)
     &             +103.d0/10.d0*u_val(i0,i1+1,j)*u_val(i0,i1+3,j)
     &             -553.d0/60.d0*u_val(i0,i1+2,j)*u_val(i0,i1+3,j)
     &             +61.d0/45.d0*u_val(i0,i1+3,j)**2

            w_d(0) = 1.d0/(20.d0*(1.0d-6+s_d(0))**2)
            w_d(1) = 9.d0/(20.d0*(1.0d-6+s_d(1))**2)
            w_d(2) = 9.d0/(20.d0*(1.0d-6+s_d(2))**2)
            w_d(3) = 1.d0/(20.d0*(1.0d-6+s_d(3))**2)
          endif

            inv_tot = 0.d0
            do k=0,k_d
              inv_tot = inv_tot + w_d(k)
            enddo
            inv_tot = 1.d0/inv_tot

            if(j == 0) then
              du(i0,i1,1) = 0.d0
              do k=0,k_d
                w_d(k) = inv_tot*w_d(k)
                du(i0,i1,1) = du(i0,i1,1) + w_d(k)*u_d(k)
              enddo
            else
              dv(i0,i1,1) = 0.d0
              do k=0,k_d
                w_d(k) = inv_tot*w_d(k)
                dv(i0,i1,1) = dv(i0,i1,1) + w_d(k)*u_d(k)
              enddo
            endif
          enddo
        enddo
      enddo
      end subroutine
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       Computes r = grad q
c
c       where d is vector valued cell centered
c       and q is scalar valued cell centered
c       weno centered differences based on
c       speeds stored in u (vector valued cell centered)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine compute_grad_weno_2d(
     &                     q_data, q_gcw,
     &                     r_data, r_gcw,
     &                     u_data, u_gcw,
     &                     ilower0, ilower1,
     &                     iupper0, iupper1,
     &                     k, dx,
     &                     interp_weights,
     &                     smooth_weights)
      implicit none
      INTEGER k
      INTEGER ilower0, iupper0
      INTEGER ilower1, iupper1

      INTEGER q_gcw
      REAL q_data(CELL2d(ilower,iupper,q_gcw))

      INTEGER u_gcw
      REAL u_data(CELL2d(ilower,iupper,u_gcw),0:1)

      INTEGER r_gcw
      REAL r_data(CELL2d(ilower,iupper,r_gcw),0:1)

      REAL interp_weights(0:(k-1),0:(k-1))
      REAL smooth_weights(0:k-1)

      REAL interp_values(0:k-1)
      REAL smooth_id(0:k-1)
      REAL weights(0:k-1)
      REAL s_data_0(FACE2d0(ilower,iupper,0))
      REAL s_data_1(FACE2d1(ilower,iupper,0))

      INTEGER i0, i1
      INTEGER j, r, step

      REAL dx(0:1)
      REAL eps, total, alpha
      REAL u_ij, v_ij
      eps = 1.0d-6

      call reconstruct_data_on_patch_upwind_2d(q_data, q_gcw,
     &         s_data_0, s_data_1, 0,
     &         u_data, u_gcw,
     &         ilower0, ilower1, iupper0, iupper1,
     &         interp_weights,
     &         smooth_weights, k)

       do i0 = ilower0,iupper0
         do i1 = ilower1,iupper1
           r_data(i0,i1,0) = (s_data_0(i0+1,i1)-s_data_0(i0,i1))/dx(0)
           r_data(i0,i1,1) = (s_data_1(i1+1,i0)-s_data_1(i1,i0))/dx(1)
         enddo
       enddo
       end subroutine

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       Computes r = grad q
c
c       where d is vector valued cell centered
c       and q is scalar valued cell centered
c       using an adaptive stencil based on
c       smoothness indicators of interpolating polynomials
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine compute_grad_adapt_2d(dx,
     &         q_data, q_gcw, r_data, r_gcw,
     &         ilower0, iupper0, ilower1, iupper1, k)

      implicit none
      INTEGER ilower0, iupper0
      INTEGER ilower1, iupper1
      INTEGER k

      REAL dx(0:1)
c
c      Q Data
c
      INTEGER q_gcw
      REAL q_data(CELL2d(ilower,iupper,q_gcw))
c
c      Return Data
c
      INTEGER r_gcw
      REAL r_data(CELL2d(ilower,iupper,r_gcw),0:1)

      INTEGER i0, i1
      INTEGER j
      REAL alpha(0:(k-1))
      REAL dq_dx(0:(k-1)), dq_dy(0:(k-1))
      REAL si(0:(k-1))
      REAL u_ij, v_ij
      REAL eps, inv_tot
      REAL scale_x, scale_y
      scale_x = 1.0/(dx(0))
      scale_y = 1.0/(dx(1))
      eps = 1.0d-6

      do i1 = ilower1, iupper1
        do i0 = ilower0, iupper0
         if (k == 3) then
c         Do WENO interpolation in x derivative
         dq_dx(0) = (3.d0*q_data(i0,i1)-4.d0*q_data(i0-1,i1)+
     &              q_data(i0-2,i1))*scale_x*0.5d0
         dq_dx(1) = (q_data(i0+1,i1)-q_data(i0-1,i1))
     &             *scale_x*0.5d0
         dq_dx(2) = (-3.d0*q_data(i0,i1)+4.d0*q_data(i0+1,i1)-
     &              q_data(i0+2,i1))*scale_x*0.5d0
c         Calculate smoothness indicators
         si(0) = 13.d0/12.d0*
     &      (q_data(i0-2,i1)-2.d0*q_data(i0-1,i1)+q_data(i0,i1))**2
     &      +0.25d0*
     &      (q_data(i0-2,i1)-4.d0*q_data(i0-1,i1)+3.d0*q_data(i0,i1))**2
         si(1) = 13.d0/12.d0*
     &      (q_data(i0-1,i1)-2.d0*q_data(i0,i1)+q_data(i0+1,i1))**2
     &      +0.25d0*
     &      (q_data(i0-1,i1)-q_data(i0+1,i1))**2
         si(2) = 13.d0/12.d0*
     &      (q_data(i0,i1)-2.d0*q_data(i0+1,i1)+q_data(i0+2,i1))**2
     &      +0.25d0*
     &      (3.d0*q_data(i0,i1)-4.d0*q_data(i0+1,i1)+q_data(i0+2,i1))**2
         alpha(0) = 1.d0/(6.d0*(eps+si(0))**1)
         alpha(1) = 2.d0/(3.d0*(eps+si(1))**1)
         alpha(2) = 1.d0/(6.d0*(eps+si(2))**1)
         inv_tot = 1.d0/(alpha(0) + alpha(1) + alpha(2))
         alpha(0) = alpha(0)*inv_tot
         alpha(1) = alpha(1)*inv_tot
         alpha(2) = alpha(2)*inv_tot
         r_data(i0,i1,0) = (alpha(0)*dq_dx(0)
     &                 +alpha(1)*dq_dx(1)
     &                 +alpha(2)*dq_dx(2))
c        Do WENO interpolation in y derivative
         dq_dy(0) = (3.d0*q_data(i0,i1)-4.d0*q_data(i0,i1-1)+
     &              q_data(i0,i1-2))*scale_y*0.5d0
         dq_dy(1) = (q_data(i0,i1+1)-q_data(i0,i1-1))
     &             *scale_y*0.5d0
         dq_dy(2) = (-3.d0*q_data(i0,i1)+4.d0*q_data(i0,i1+1)-
     &             q_data(i0,i1+2))*scale_y*0.5d0
c         Calculate smoothness indicators
         si(0) = 13.d0/12.d0*
     &      (q_data(i0,i1-2)-2.d0*q_data(i0,i1-1)+q_data(i0,i1))**2
     &      +0.25d0*
     &      (q_data(i0,i1-2)-4.d0*q_data(i0,i1-1)+3.d0*q_data(i0,i1))**2
         si(1) = 13.d0/12.d0*
     &      (q_data(i0,i1-1)-2.d0*q_data(i0,i1)+q_data(i0,i1+1))**2
     &      +0.25d0*
     &      (q_data(i0,i1-1)-q_data(i0,i1+1))**2
         si(2) = 13.d0/12.d0*
     &      (q_data(i0,i1)-2.d0*q_data(i0,i1+1)+q_data(i0,i1+2))**2
     &      +0.25d0*
     &      (3.d0*q_data(i0,i1)-4.d0*q_data(i0,i1+1)+q_data(i0,i1+2))**2
         alpha(0) = 1.d0/(6.d0*(eps+si(0))**2)
         alpha(1) = 2.d0/(3.d0*(eps+si(1))**2)
         alpha(2) = 1.d0/(6.d0*(eps+si(2))**2)
         inv_tot = 1.d0/(alpha(0) + alpha(1) + alpha(2))
         alpha(0) = alpha(0)*inv_tot
         alpha(1) = alpha(1)*inv_tot
         alpha(2) = alpha(2)*inv_tot
         r_data(i0,i1,1) =
     &                 (alpha(0)*dq_dy(0)
     &                 +alpha(1)*dq_dy(1)
     &                 +alpha(2)*dq_dy(2))
         elseif (k == 2) then
           dq_dx(0) = scale_x*(q_data(i0,i1)-q_data(i0-1,i1))
           dq_dx(1) = scale_x*(q_data(i0+1,i1)-q_data(i0,i1))

           si(0) = (q_data(i0-1,i1)-q_data(i0,i1))**2
           si(1) = (q_data(i0,i1)-q_data(i0+1,i1))**2

           alpha(0) = 0.5d0/((eps+si(0))**1)
           alpha(1) = 0.5d0/((eps+si(1))**1)
           inv_tot = 1.d0/(alpha(0)+alpha(1))
           alpha(0) = alpha(0)*inv_tot
           alpha(1) = alpha(1)*inv_tot
           r_data(i0,i1,0) = (alpha(0)*dq_dx(0)
     &                          +alpha(1)*dq_dx(1))

           dq_dx(0) = scale_y*(q_data(i0,i1)-q_data(i0,i1-1))
           dq_dx(1) = scale_y*(q_data(i0,i1+1)-q_data(i0,i1))

           si(0) = (q_data(i0,i1-1)-q_data(i0,i1))**2
           si(1) = (q_data(i0,i1)-q_data(i0,i1+1))**2

           alpha(0) = 0.5d0/((eps+si(0))**1)
           alpha(1) = 0.5d0/((eps+si(1))**1)
           inv_tot = 1.d0/(alpha(0)+alpha(1))
           alpha(0) = alpha(0)*inv_tot
           alpha(1) = alpha(1)*inv_tot
           r_data(i0,i1,1) =
     &                     (alpha(0)*dq_dx(0)
     &                     +alpha(1)*dq_dx(1))
         endif
        enddo
      enddo
      end subroutine
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       Computes r = grad q
c
c       where r is vector valued cell centered
c       and q is scalar valued cell centered
c       using second order finite differences
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine compute_grad_cent_2d(dx,
     &         q_data, q_gcw, r_data, r_gcw,
     &         ilower0, iupper0,
     &         ilower1, iupper1, k)
      implicit none
      INTEGER ilower0, iupper0
      INTEGER ilower1, iupper1
      INTEGER k

      REAL dx(0:1)
c
c      Q Data
c
      INTEGER q_gcw
      REAL q_data(CELL2d(ilower,iupper,q_gcw))
c
c      Return Data
c
      INTEGER r_gcw
      REAL r_data(CELL2d(ilower,iupper,r_gcw),0:1)

      INTEGER i0, i1, i2
      REAL scale_x, scale_y
      scale_x = 1.d0/(2.d0*dx(0))
      scale_y = 1.d0/(2.d0*dx(1))

      do i1 = ilower1, iupper1
        do i0 = ilower0, iupper0
          r_data(i0,i1,0) = scale_x*
     &      (q_data(i0+1,i1)-q_data(i0-1,i1))
          r_data(i0,i1,1) = scale_y*
     &      (q_data(i0,i1+1)-q_data(i0,i1-1))
        enddo
      enddo
      endsubroutine

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       Reconstructs data on patches using a weno scheme
c       the convex and interpolation weights must be supplied
c
c       q_data is cell centered with depth 1
c       r_data_* are face centered with depth 1
c         and returns the values reconstructed from q_data
c         using upwinding from speeds stored in u_data
c       u_data is cell centered with depth 2
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine reconstruct_data_on_patch_upwind_2d(q_data, q_gcw,
     &            r_data_0, r_data_1, r_gcw,
     &            u_data, u_gcw,
     &            ilower0, ilower1, iupper0, iupper1,
     &            weights_cell_sides,
     &            smooth_weights_sides,
     &            k)
      implicit none
      INTEGER k
      INTEGER ilower0, iupper0
      INTEGER ilower1, iupper1

      INTEGER q_gcw
      REAL q_data(CELL2d(ilower,iupper,q_gcw))

      INTEGER u_gcw
      REAL u_data(CELL2d(ilower,iupper,u_gcw),0:1)

      INTEGER r_gcw
      REAL r_data_0(FACE2d0(ilower,iupper,r_gcw))
      REAL r_data_1(FACE2d1(ilower,iupper,r_gcw))

      REAL weights_cell_sides(0:(k-1),0:(k-1))
      REAL smooth_weights_sides(0:k-1)

      REAL interp_values(0:k-1)
      REAL smooth_id(0:k-1)
      REAL weights(0:k-1)
      REAL s_vals_x(ilower0-k:(iupper0+k),
     &               (ilower1):(iupper1))
      REAL s_vals_y(ilower1-k:(iupper1+k),
     &               (ilower0):(iupper0))

      INTEGER i0, i1
      INTEGER j, r, step

      REAL eps, total, alpha
      REAL u_ij, v_ij
      eps = 1.0d-6
c     FIRST INTERPOLANT
c     X DIRECTION
      do i1 = ilower1, iupper1
       do i0=ilower0,iupper0+1
         u_ij = u_data(i0,i1,0)
         if(i0 == iupper0+1) then
           u_ij = u_data(i0-1,i1,0)
         endif
         do r=0,k-1
           interp_values(r) = 0.d0
           do j=0,k-1
             if(u_ij > 0) then
               step = -1-r+j
             else
               step = r-j
             endif
               interp_values(r) = interp_values(r)
     &           + weights_cell_sides(r,j)*q_data(i0+step,i1)
           enddo
         enddo
         if(u_ij > 0) then
           smooth_id(0) = 13.d0/12.d0*(q_data(i0-1,i1)
     &         -2.d0*q_data(i0,i1)+q_data(i0+1,i1))**2
     &       + 0.25d0*(3.d0*q_data(i0-1,i1)-4.d0*q_data(i0,i1)
     &         +q_data(i0+1,i1))**2
           smooth_id(1) = 13.d0/12.d0*(q_data(i0-2,i1)
     &         -2.d0*q_data(i0-1,i1)+q_data(i0,i1))**2
     &       + 0.25d0*(q_data(i0-2,i1)-q_data(i0,i1))**2
           smooth_id(2) = 13.d0/12.d0*(q_data(i0-3,i1)
     &         -2.d0*q_data(i0-2,i1)+q_data(i0-1,i1))**2
     &       + 0.25d0*(3.d0*q_data(i0-1,i1)
     &         -4.d0*q_data(i0-2,i1)+q_data(i0-3,i1))**2
        else
          smooth_id(0) = 13.d0/12.d0*(q_data(i0,i1)
     &         -2.d0*q_data(i0-1,i1)+q_data(i0-2,i1))**2
     &       + 0.25d0*(3.d0*q_data(i0,i1)-4.d0*q_data(i0-1,i1)
     &         +q_data(i0-2,i1))**2
           smooth_id(1) = 13.d0/12.d0*(q_data(i0+1,i1)
     &         -2.d0*q_data(i0,i1)+q_data(i0-1,i1))**2
     &       + 0.25d0*(q_data(i0+1,i1)-q_data(i0-1,i1))**2
           smooth_id(2) = 13.d0/12.d0*(q_data(i0+2,i1)
     &         -2.d0*q_data(i0+1,i1)+q_data(i0,i1))**2
     &       + 0.25d0*(3.d0*q_data(i0,i1)
     &         -4.d0*q_data(i0+1,i1)+q_data(i0+2,i1))**2
        endif
        total = 0.d0
          do j=0,k-1
            alpha = smooth_weights_sides(j)
     &          /(eps+smooth_id(j))**2
            total = total + alpha
            weights(j) = alpha
          enddo
          do j=0,k-1
            weights(j) = weights(j)/total
          enddo
          r_data_0(i0,i1) = 0.d0
          do r=0,k-1
            r_data_0(i0,i1) = r_data_0(i0,i1)
     &              + weights(r)*interp_values(r)
          enddo
        enddo
      enddo
c     INTERPOLATE IN Y DIRECTION
      do i1 = ilower1,iupper1+1
        do i0 = ilower0,iupper0
          v_ij = u_data(i0,i1,1)
          if(i1 == iupper1+1) then
            v_ij = u_data(i0,i1-1,1)
          endif
          do r=0,k-1
            interp_values(r) = 0.d0
            do j=0,k-1
              if(v_ij > 0) then
                step = -1-r+j
              else
                step = r-j
              endif
              interp_values(r) = interp_values(r)
     &          + weights_cell_sides(r,j)*q_data(i0,i1+step)
            enddo
          enddo
          if(v_ij > 0) then
            smooth_id(0) = 13.d0/12.d0*(q_data(i0,i1-1)
     &          -2.d0*q_data(i0,i1)+q_data(i0,i1+1))**2
     &        + 0.25d0*(3.d0*q_data(i0,i1-1)-4.d0*q_data(i0,i1)
     &          +q_data(i0,i1+1))**2
            smooth_id(1) = 13.d0/12.d0*(q_data(i0,i1-2)
     &          -2.d0*q_data(i0,i1-1)+q_data(i0,i1))**2
     &        + 0.25d0*(q_data(i0,i1-2)-q_data(i0,i1))**2
            smooth_id(2) = 13.d0/12.d0*(q_data(i0,i1-3)
     &          -2.d0*q_data(i0,i1-2)+q_data(i0,i1-1))**2
     &        + 0.25d0*(3.d0*q_data(i0,i1-1)
     &          -4.d0*q_data(i0,i1-2)+q_data(i0,i1-3))**2
          else
            smooth_id(0) = 13.d0/12.d0*(q_data(i0,i1)
     &          -2.d0*q_data(i0,i1-1)+q_data(i0,i1-2))**2
     &        + 0.25d0*(3.d0*q_data(i0,i1)-4.d0*q_data(i0,i1-1)
     &          +q_data(i0,i1-2))**2
            smooth_id(1) = 13.d0/12.d0*(q_data(i0,i1+1)
     &          -2.d0*q_data(i0,i1)+q_data(i0,i1-1))**2
     &        + 0.25d0*(q_data(i0,i1+1)-q_data(i0,i1-1))**2
            smooth_id(2) = 13.d0/12.d0*(q_data(i0,i1+2)
     &          -2.d0*q_data(i0,i1+1)+q_data(i0,i1))**2
     &        + 0.25d0*(3.d0*q_data(i0,i1)
     &          -4.d0*q_data(i0,i1+1)+q_data(i0,i1+2))**2
          endif
          total = 0.d0
          do j=0,k-1
            alpha = smooth_weights_sides(j)
     &          /((eps+smooth_id(j))**2)
            total = total + alpha
            weights(j) = alpha
          enddo
          do j=0,k-1
            weights(j) = weights(j)/total
          enddo
          r_data_1(i1,i0) = 0.d0
          do r=0,k-1
            r_data_1(i1,i0) = r_data_1(i1,i0)
     &        + weights(r)*interp_values(r)
          enddo
        enddo
      enddo
      end subroutine
