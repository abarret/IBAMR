c234567
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes cell centered Oldroyd-B type Convective Operator
c
c     where u is vector valued face centered velocity
c     and tau is the symmetric stress tensor cell centered
c     c_data is convective term u.grad(tau)
c     returns s_data
c     computes grad(u) using centered differencing
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine conform_tens_conv_u_f_oper_2d
     &        (dx, u_data_0, u_data_1,
     &        u_gcw, s_data, s_gcw, rhs_data, rhs_gcw,
     &        c_data, c_gcw, r_data, r_gcw,
     &        ilower0, iupper0, ilower1, iupper1, lambda)
      implicit none
      Integer ilower0, iupper0
      Integer ilower1, iupper1
      real*8 lambda
      real*8 dx(0:1)

c
c     Velocity Data
c
      Integer u_gcw
      real*8 u_data_0((ilower0-u_gcw):(iupper0+u_gcw+1),
     &                (ilower1-u_gcw):(iupper1+u_gcw))
      real*8 u_data_1((ilower1-u_gcw):(iupper1+u_gcw+1),
     &                (ilower0-u_gcw):(iupper0+u_gcw))
c
c     Tensor Data
c
      Integer s_gcw
      real*8 s_data((ilower0-s_gcw):(iupper0+s_gcw),
     &              (ilower1-s_gcw):(iupper1+s_gcw),
     &              0:2)
c
c     RHS Data
c
      Integer rhs_gcw
      real*8 rhs_data((ilower0-rhs_gcw):(iupper0+rhs_gcw),
     &                (ilower1-rhs_gcw):(iupper1+rhs_gcw),
     &                 0:2)
c
c     Return Data
c
      Integer r_gcw
      real*8 r_data((ilower0-r_gcw):(iupper0+r_gcw),
     &              (ilower1-r_gcw):(iupper1+r_gcw),
     &              0:2)
c
c     Convective Data
c
      Integer c_gcw
      real*8 c_data((ilower0-c_gcw):(iupper0+c_gcw),
     &              (ilower1-c_gcw):(iupper1+c_gcw),
     &              0:2)

      real*8 du(ilower0:iupper0,ilower1:iupper1,0:1)
      real*8 dv(ilower0:iupper0,ilower1:iupper1,0:1)

      Integer i0, i1
      real*8 du_dx, dv_dx
      real*8 du_dy, dv_dy
      real*8 scale_ux, scale_uy
      real*8 scale_vx, scale_vy
      real*8 qxx, qyy, qxy
      real*8 l_inv

      l_inv = 1.d0/lambda
      scale_vx = 1.d0/(4.d0*dx(0))
      scale_vy = 1.d0/dx(1)
      scale_ux = 1.d0/dx(0)
      scale_uy = 1.d0/(4.d0*dx(1))
!       scale_ux = 1.d0/(24.d0*dx(0))
!       scale_vy = 1.d0/(24.d0*dx(1))
!       scale_uy = 1.d0/(12.d0*dx(1))
!       scale_vx = 1.d0/(12.d0*dx(0))

      call weno_grad_s_to_c(u_data_0, u_data_1, u_gcw, du, dv, 0, dx,
     &               0.5d0, ilower0, iupper0, ilower1, iupper1,3,3)

      do i1 = ilower1, iupper1
         do i0 = ilower0, iupper0
!           2nd order approximations to derivatives
!            du_dx = scale_ux*(u_data_0(i0+1,i1)-u_data_0(i0,i1))
!            du_dy = scale_uy*(u_data_0(i0+1,i1+1)+u_data_0(i0,i1+1)
!      &              -u_data_0(i0+1,i1-1)-u_data_0(i0,i1-1))
!            dv_dy = scale_vy*(u_data_1(i1+1,i0)-u_data_1(i1,i0))
!            dv_dx = scale_vx*(u_data_1(i1+1,i0+1)+u_data_1(i1,i0+1)
!      &              -u_data_1(i1,i0-1) - u_data_1(i1+1,i0-1))
!           du(i0,i1,0) = scale_ux*(u_data_0(i0+1,i1)-u_data_0(i0,i1))
!           du(i0,i1,1) = scale_uy*(u_data_0(i0+1,i1+1)+u_data_0(i0,i1+1)
!      &              -u_data_0(i0+1,i1-1)-u_data_0(i0,i1-1))
!           dv(i0,i1,1) = scale_vy*(u_data_1(i1+1,i0)-u_data_1(i1,i0))
!           dv(i0,i1,0) = scale_vx*(u_data_1(i1+1,i0+1)+u_data_1(i1,i0+1)
!      &              -u_data_1(i1,i0-1) - u_data_1(i1+1,i0-1))
!           4th order approximations (except for interpolation...)
!            du_dx = scale_ux*(u_data_0(i0-1,i1)-27.d0*u_data_0(i0,i1)
!      &           +27.d0*u_data_0(i0+1,i1)-u_data_0(i0+2,i1))
!            dv_dy = scale_vy*(u_data_1(i1-1,i0)-27.d0*u_data_1(i1,i0)
!      &           +27.d0*u_data_1(i1+1,i0)-u_data_1(i1+2,i0))
!            du_dy = scale_uy*(-u_data_0(i0-1,i1-2)/16.d0+9.d0
!      &           *u_data_0(i0,i1-2)/16.d0+9.d0*u_data_0(i0+1,i1-2)
!      &           /16.d0-u_data_0(i0+2,i1-2)/16.d0
!      &           -8.d0*(-u_data_0(i0-1,i1-1)/16.d0+9.d0
!      &           *u_data_0(i0,i1-1)/16.d0+9.d0*u_data_0(i0+1,i1-1)
!      &           /16.d0-u_data_0(i0+2,i1-1)/16.d0)
!      &           +8.d0*(-u_data_0(i0-1,i1+1)/16.d0+9.d0
!      &           *u_data_0(i0,i1+1)/16.d0+9.d0*u_data_0(i0+1,i1+1)
!      &           /16.d0-u_data_0(i0+2,i1+1)/16.d0)
!      &           -(-u_data_0(i0-1,i1+2)/16.d0+9.d0*u_data_0(i0,i1+2)
!      &           /16.d0+9.d0*u_data_0(i0+1,i1+2)/16.d0
!      &           -u_data_0(i0+2,i1+2)/16.d0))
!            dv_dx = scale_vx*(-u_data_1(i1-1,i0-2)/16.d0+9.d0
!      &           *u_data_1(i1,i0-2)/16.d0+9.d0*u_data_1(i1+1,i0-2)
!      &           /16.d0-u_data_1(i1+2,i0-2)/16.d0
!      &           -8.d0*(-u_data_1(i1-1,i0-1)/16.d0+9.d0
!      &           *u_data_1(i1,i0-1)/16.d0+9.d0*u_data_1(i1+1,i0-1)
!      &           /16.d0-u_data_1(i1+2,i0-1)/16.d0)
!      &           +8.d0*(-u_data_1(i1-1,i0+1)/16.d0+9.d0
!      &           *u_data_1(i1,i0+1)/16.d0+9.d0*u_data_1(i1+1,i0+1)
!      &           /16.d0-u_data_1(i1+2,i0+1)/16.d0)
!      &           -(-u_data_1(i1-1,i0+2)/16.d0+9.d0
!      &           *u_data_1(i1,i0+2)/16.d0+9.d0*u_data_1(i1+1,i0+2)
!      &           /16.d0-u_data_1(i1+2,i0+2)/16.d0))
            du_dx = du(i0,i1,0)
            du_dy = du(i0,i1,1)
            dv_dx = dv(i0,i1,0)
            dv_dy = dv(i0,i1,1)
            qxx = s_data(i0,i1,0)
            qyy = s_data(i0,i1,1)
            qxy = s_data(i0,i1,2)

            r_data(i0,i1,0) = c_data(i0,i1,0)
     &        - 2.d0*du_dx*qxx - 2.d0*du_dy*qxy
     &        - l_inv*rhs_data(i0,i1,0)
            r_data(i0,i1,1) = c_data(i0,i1,1)
     &       - 2.d0*dv_dx*qxy - 2.d0*dv_dy*qyy
     &       - l_inv*rhs_data(i0,i1,1)
            r_data(i0,i1,2) = c_data(i0,i1,2)
     &       - qyy*du_dy - qxx*dv_dx
     &       - l_inv*rhs_data(i0,i1,2)
         enddo
      enddo
      end subroutine
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes cell centered Oldroyd-B type Convective Operator
c
c     where u is vector valued cell centered velocity
c     and tau is symmetric tensor valued cell centered
c     c_data is convective term u.grad(tau)
c     returns s_data
c     computes grad(u) using centered differencing
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine conform_tens_conv_u_c_oper_2d
     &        (dx, u_data,
     &        u_gcw, s_data, s_gcw, rhs_data, rhs_gcw,
     &        c_data, c_gcw, r_data, r_gcw,
     &        ilower0, iupper0, ilower1, iupper1, lambda)
      implicit none
      Integer ilower0, iupper0
      Integer ilower1, iupper1

      real*8 lambda
      real*8 dx(0:1)

c
c     Velocity Data
c
      Integer u_gcw
      real*8 u_data((ilower0-u_gcw):(iupper0+u_gcw),
     &              (ilower1-u_gcw):(iupper1+u_gcw),
     &               0:1)
c
c     Tensor Data
c
      Integer s_gcw
      real*8 s_data((ilower0-s_gcw):(iupper0+s_gcw),
     &              (ilower1-s_gcw):(iupper1+s_gcw),
     &              0:2)
c
c     RHS Data
c
      Integer rhs_gcw
      real*8 rhs_data((ilower0-rhs_gcw):(iupper0+rhs_gcw),
     &                (ilower1-rhs_gcw):(iupper1+rhs_gcw),
     &                 0:2)
c
c     Return Data
c
      Integer r_gcw
      real*8 r_data((ilower0-r_gcw):(iupper0+r_gcw),
     &              (ilower1-r_gcw):(iupper1+r_gcw),
     &              0:2)
c
c     Convective Data
c
      Integer c_gcw
      real*8 c_data((ilower0-c_gcw):(iupper0+c_gcw),
     &              (ilower1-c_gcw):(iupper1+c_gcw),
     &              0:2)

      Integer i0, i1
      real*8 u_ij, v_ij
      real*8 dq_00_dx, dq_11_dx, dq_01_dx
      real*8 dq_00_dy, dq_11_dy, dq_01_dy
      real*8 du(ilower0:iupper0,
     &          ilower1:iupper1,
     &          0:1)
      real*8 dv(ilower0:iupper0,
     &          ilower1:iupper1,
     &          0:1)
      real*8 scale0, scale1
      real*8 l_inv

      call compute_grad_adapt_2d(dx,
     &   u_data(:,:,0), u_gcw, du, 0,
     &   ilower0, iupper0, ilower1, iupper1, 3)
      call compute_grad_adapt_2d(dx,
     &   u_data(:,:,1), u_gcw, dv, 0,
     &   ilower0, iupper0, ilower1, iupper1, 3)

      l_inv = 1.d0/lambda
      scale0 = 1.d0/(2.d0*dx(0))
      scale1 = 1.d0/(2.d0*dx(1))
      do i1 = ilower1, iupper1
         do i0 = ilower0, iupper0
            r_data(i0,i1,0) = c_data(i0,i1,0) -
     &       2.d0*du(i0,i1,0)*s_data(i0,i1,0) -
     &       2.d0*du(i0,i1,1)*s_data(i0,i1,2)
     &       - l_inv*rhs_data(i0,i1,0)
            r_data(i0,i1,1) = c_data(i0,i1,1) -
     &       2.d0*dv(i0,i1,0)*s_data(i0,i1,2) -
     &       2.d0*dv(i0,i1,1)*s_data(i0,i1,1)
     &       - l_inv*rhs_data(i0,i1,1)
            r_data(i0,i1,2) = c_data(i0,i1,2) -
     &      s_data(i0,i1,1)*du(i0,i1,1) -
     &      s_data(i0,i1,0)*dv(i0,i1,0)
     &      - l_inv*rhs_data(i0,i1,2)
         enddo
      enddo
      end subroutine
