define(NDIM,2)dnl
define(REAL,`double precision')dnl
define(INTEGER,`integer')dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim2d.i)dnl
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
      subroutine conform_tens_conv_u_s_oper_2d
     &        (dx, u_data_0, u_data_1,
     &        u_gcw, s_data, s_gcw, rhs_data, rhs_gcw,
     &        c_data, c_gcw, r_data, r_gcw,
     &        ilower0, iupper0, ilower1, iupper1, lambda)
      implicit none
      INTEGER ilower0, iupper0
      INTEGER ilower1, iupper1
      REAL lambda
      REAL dx(0:1)

c
c     Velocity Data
c
      INTEGER u_gcw
      REAL u_data_0(SIDE2d0(ilower,iupper,u_gcw))
      REAL u_data_1(SIDE2d1(ilower,iupper,u_gcw))
c
c     Tensor Data
c
      INTEGER s_gcw
      REAL s_data(CELL2d(ilower,iupper,s_gcw),0:2)
c
c     RHS Data
c
      INTEGER rhs_gcw
      REAL rhs_data(CELL2d(ilower,iupper,rhs_gcw),0:2)
c
c     Return Data
c
      INTEGER r_gcw
      REAL r_data(CELL2d(ilower,iupper,r_gcw),0:2)
c
c     Convective Data
c
      INTEGER c_gcw
      REAL c_data(CELL2d(ilower,iupper,c_gcw),0:2)

      REAL du(CELL2d(ilower,iupper,0),0:1)
      REAL dv(CELL2d(ilower,iupper,0),0:1)

      INTEGER i0, i1
      REAL du_dx, dv_dx
      REAL du_dy, dv_dy
      REAL scale_ux, scale_uy
      REAL scale_vx, scale_vy
      REAL qxx, qyy, qxy
      REAL l_inv

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
           du_dx = scale_ux*(u_data_0(i0+1,i1)-u_data_0(i0,i1))
           du_dy = scale_uy*(u_data_0(i0+1,i1+1)+u_data_0(i0,i1+1)
     &              -u_data_0(i0+1,i1-1)-u_data_0(i0,i1-1))
           dv_dy = scale_vy*(u_data_1(i0,i1+1)-u_data_1(i0,i1))
           dv_dx = scale_vx*(u_data_1(i0+1,i1+1)+u_data_1(i0+1,i1)
     &              -u_data_1(i0-1,i1) - u_data_1(i0-1,i1+1))
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
!             du_dx = du(i0,i1,0)
!             du_dy = du(i0,i1,1)
!             dv_dx = dv(i0,i1,0)
!             dv_dy = dv(i0,i1,1)
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
      INTEGER ilower0, iupper0
      INTEGER ilower1, iupper1

      REAL lambda
      REAL dx(0:1)

c
c     Velocity Data
c
      INTEGER u_gcw
      REAL u_data(CELL2d(ilower,iupper,u_gcw),0:1)
c
c     Tensor Data
c
      INTEGER s_gcw
      REAL s_data(CELL2d(ilower,iupper,s_gcw),0:2)
c
c     RHS Data
c
      INTEGER rhs_gcw
      REAL rhs_data(CELL2d(ilower,iupper,rhs_gcw),0:2)
c
c     Return Data
c
      INTEGER r_gcw
      REAL r_data(CELL2d(ilower,iupper,r_gcw),0:2)
c
c     Convective Data
c
      INTEGER c_gcw
      REAL c_data(CELL2d(ilower,iupper,c_gcw),0:2)

      INTEGER i0, i1
      REAL u_ij, v_ij
      REAL dq_00_dx, dq_11_dx, dq_01_dx
      REAL dq_00_dy, dq_11_dy, dq_01_dy
      REAL du(CELL2d(ilower,iupper,0),0:1)
      REAL dv(CELL2d(ilower,iupper,0),0:1)
      REAL scale0, scale1
      REAL l_inv

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