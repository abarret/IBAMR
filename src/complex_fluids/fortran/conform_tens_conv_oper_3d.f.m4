define(NDIM,3)dnl
define(REAL,`double precision')dnl
define(INTEGER,`integer')dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim3d.i)dnl
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes cell centered Oldroyd-B type Convective Operator
c
c     where u is vector valued face centered velocity
c     and tau is symmetric tensor valued cell centered
c     c_data is convective term u.grad(tau)
c     returns s_data
c     computes grad(u) using centered differencing
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine conform_tens_conv_u_f_oper_3d
     &        (dx, u_data_0, u_data_1, u_data_2,
     &        u_gcw, s_data, s_gcw, rhs_data, rhs_gcw,
     &        c_data, c_gcw, r_data, r_gcw,
     &        ilower0, iupper0, ilower1,
     &        iupper1, ilower2, iupper2, lambda)
      implicit none

      INTEGER ilower0, iupper0
      INTEGER ilower1, iupper1
      INTEGER ilower2, iupper2

      REAL lambda
      REAL dx(0:2)

c
c    Velocity Data
c
      INTEGER u_gcw
      REAL u_data_0(FACE3d0(ilower,iupper,u_gcw))
      REAL u_data_1(FACE3d1(ilower,iupper,u_gcw))
      REAL u_data_2(FACE3d2(ilower,iupper,u_gcw))
c
c    Tensor Data
c
      INTEGER s_gcw
      REAL s_data(CELL3d(ilower,iupper,s_gcw),0:5)
c
c    RHS Data
c
      INTEGER rhs_gcw
      REAL rhs_data(CELL3d(ilower,iupper,rhs_gcw),0:5)
c
c    Convec Data
c
      INTEGER c_gcw
      REAL c_data(CELL3d(ilower,iupper,c_gcw),0:5)
c
c    Return Data
c
      INTEGER r_gcw
      REAL r_data(CELL3d(ilower,iupper,r_gcw),0:5)

      REAL du(CELL3d(ilower,iupper,0),0:2)
      REAL dv(CELL3d(ilower,iupper,0),0:2)
      REAL dw(CELL3d(ilower,iupper,0),0:2)

      INTEGER i0, i1, i2
      REAL du_dx, dv_dx, dw_dx
      REAL du_dy, dv_dy, dw_dy
      REAL du_dz, dv_dz, dw_dz
      REAL scale_ux, scale_uy, scale_uz
      REAL scale_vx, scale_vy, scale_vz
      REAL scale_wx, scale_wy, scale_wz
      REAL l_inv
      REAL qxx_ij, qyy_ij, qzz_ij
      REAL qyz_ij, qxz_ij, qxy_ij

      l_inv = 1.d0/lambda
      scale_ux = 1.d0/dx(0)
      scale_uy = 1.d0/(4.d0*dx(1))
      scale_uz = 1.d0/(4.d0*dx(2))
      scale_vx = 1.d0/(4.d0*dx(0))
      scale_vy = 1.d0/dx(1)
      scale_vz = 1.d0/(4.d0*dx(2))
      scale_wx = 1.d0/(4.d0*dx(0))
      scale_wy = 1.d0/(4.d0*dx(1))
      scale_wz = 1.d0/dx(2)

!       call weno_grad_s_to_c_3d(u_data_0, u_data_1, u_data_2, u_gcw,
!      &                         du, dv, dw, 0, dx, 0.5d0,
!      &                         ilower0, iupper0, ilower1, iupper1,
!      &                         ilower2, iupper2, 3, 3)

      do i2 = ilower2, iupper2
        do i1 = ilower1, iupper1
          do i0 = ilower0, iupper0
            du_dx = 
     &        scale_ux*(u_data_0(i0+1,i1,i2)-u_data_0(i0,i1,i2))
            du_dy = 
     &        scale_uy*(u_data_0(i0+1,i1+1,i2)+u_data_0(i0,i1+1,i2)
     &              -u_data_0(i0+1,i1-1,i2)-u_data_0(i0,i1-1,i2))
            du_dz = 
     &        scale_uz*(u_data_0(i0+1,i1,i2+1)+u_data_0(i0,i1,i2+1)
     &              -u_data_0(i0+1,i1,i2-1)-u_data_0(i0,i1,i2-1))
            dv_dy = 
     &        scale_vy*(u_data_1(i1+1,i2,i0)-u_data_1(i1,i2,i0))
            dv_dx = 
     &        scale_vx*(u_data_1(i1+1,i2,i0+1)+u_data_1(i1,i2,i0+1)
     &              -u_data_1(i1+1,i2,i0-1) - u_data_1(i1,i2,i0-1))
            dv_dz = 
     &        scale_vz*(u_data_1(i1+1,i2+1,i0)+u_data_1(i1,i2+1,i0)
     &              -u_data_1(i1+1,i2-1,i0) - u_data_1(i1,i2-1,i0))
            dw_dx = 
     &        scale_wx*(u_data_2(i2+1,i0+1,i1)+u_data_2(i2,i0+1,i1)
     &              -u_data_2(i2+1,i0-1,i1)-u_data_2(i2,i0-1,i1))
            dw_dy = 
     &        scale_wy*(u_data_2(i2+1,i0,i1+1)+u_data_2(i2,i0,i1+1)
     &              -u_data_2(i2+1,i0,i1-1)-u_data_2(i2,i0,i1-1))
            dw_dz = 
     &        scale_wz*(u_data_2(i2+1,i0,i1)-u_data_2(i2,i0,i1))

            qxx_ij = s_data(i0,i1,i2,0)
            qyy_ij = s_data(i0,i1,i2,1)
            qzz_ij = s_data(i0,i1,i2,2)
            qyz_ij = s_data(i0,i1,i2,3)
            qxz_ij = s_data(i0,i1,i2,4)
            qxy_ij = s_data(i0,i1,i2,5)

            r_data(i0,i1,i2,0) =
     &        c_data(i0,i1,i2,0)
     &        - 2.d0*du_dx*qxx_ij - 2.d0*du_dy*qxy_ij
     &        - 2.d0*du_dz*qxz_ij
     &        - l_inv*rhs_data(i0,i1,i2,0)
            r_data(i0,i1,i2,1) =
     &        c_data(i0,i1,i2,1)
     &        - 2.d0*dv_dx*qxy_ij - 2.d0*dv_dy*qyy_ij
     &        - 2.d0*dv_dz*qyz_ij
     &        - l_inv*rhs_data(i0,i1,i2,1)
            r_data(i0,i1,i2,2) =
     &        c_data(i0,i1,i2,2)
     &        - 2.d0*dw_dx*qxz_ij - 2.d0*dw_dy*qyz_ij
     &        - 2.d0*dw_dz*qzz_ij
     &        - l_inv*rhs_data(i0,i1,i2,2)
            r_data(i0,i1,i2,3) =
     &        c_data(i0,i1,i2,3)
     &        - qyy_ij*dw_dy - qzz_ij*dv_dz
     &        + qyz_ij*du_dx - qxz_ij*dv_dx
     &        - qxy_ij*dw_dx
     &        - l_inv*rhs_data(i0,i1,i2,3)
            r_data(i0,i1,i2,4) =
     &        c_data(i0,i1,i2,4)
     &        - qxx_ij*dw_dx - qzz_ij*du_dz
     &        + qxz_ij*dv_dy - qyz_ij*du_dy
     &        - qxy_ij*dw_dy
     &        - l_inv*rhs_data(i0,i1,i2,4)
            r_data(i0,i1,i2,5) =
     &        c_data(i0,i1,i2,5)
     &        - qxx_ij*dv_dx - qyy_ij*du_dy
     &        + qxy_ij*dw_dz - qxz_ij*dv_dz
     &        - qyz_ij*du_dz
     &        - l_inv*rhs_data(i0,i1,i2,5)
          enddo
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

      subroutine conform_tens_conv_u_c_oper_3d
     &        (dx, u_data,
     &        u_gcw, s_data, s_gcw, rhs_data, rhs_gcw,
     &        c_data, c_gcw, r_data, r_gcw,
     &        ilower0, iupper0, ilower1,
     &        iupper1, ilower2, iupper2, lambda)
      implicit none

      INTEGER ilower0, iupper0
      INTEGER ilower1, iupper1
      INTEGER ilower2, iupper2

      REAL lambda
      REAL dx(0:2)

c
c    Velocity Data
c
      INTEGER u_gcw
      REAL u_data(CELL3d(ilower,iupper,u_gcw),0:2)
c
c    Tensor Data
c
      INTEGER s_gcw
      REAL s_data(CELL3d(ilower,iupper,s_gcw),0:5)
c
c    RHS Data
c
      INTEGER rhs_gcw
      REAL rhs_data(CELL3d(ilower,iupper,rhs_gcw),0:5)
c
c    Convec Data
c
      INTEGER c_gcw
      REAL c_data(CELL3d(ilower,iupper,c_gcw),0:5)
c
c    Return Data
c
      INTEGER r_gcw
      REAL r_data(CELL3d(ilower,iupper,r_gcw),0:5)

      INTEGER i0, i1, i2
      REAL du(CELL3d(ilower,iupper,0),0:2)
      REAL dv(CELL3d(ilower,iupper,0),0:2)
      REAL dw(CELL3d(ilower,iupper,0),0:2)
      REAL scale0, scale1, scale2
      REAL l_inv
      REAL qxx_ij, qyy_ij, qzz_ij
      REAL qyz_ij, qxz_ij, qxy_ij

      l_inv = 1.d0/lambda
      scale0 = 1.0/(2.0*dx(0))
      scale1 = 1.0/(2.0*dx(1))
      scale2 = 1.0/(2.0*dx(2))

      call compute_grad_adapt_3d(dx,
     &      u_data(:,:,:,0), u_gcw,
     &      du, 0, ilower0, iupper0,
     &      ilower1, iupper1,
     &      ilower2, iupper2, 3)
      call compute_grad_adapt_3d(dx,
     &      u_data(:,:,:,1), u_gcw,
     &      dv, 0, ilower0, iupper0,
     &      ilower1, iupper1,
     &      ilower2, iupper2, 3)
      call compute_grad_adapt_3d(dx,
     &      u_data(:,:,:,2), u_gcw,
     &      dw, 0, ilower0, iupper0,
     &      ilower1, iupper1,
     &      ilower2, iupper2, 3)

      do i2 = ilower2, iupper2
        do i1 = ilower1, iupper1
          do i0 = ilower0, iupper0

            qxx_ij = s_data(i0,i1,i2,0)
            qyy_ij = s_data(i0,i1,i2,1)
            qzz_ij = s_data(i0,i1,i2,2)
            qyz_ij = s_data(i0,i1,i2,3)
            qxz_ij = s_data(i0,i1,i2,4)
            qxy_ij = s_data(i0,i1,i2,5)

            r_data(i0,i1,i2,0) =
     &        c_data(i0,i1,i2,0)
     &        - 2.d0*du(i0,i1,i2,0)*qxx_ij - 2.d0*du(i0,i1,i2,1)*qxy_ij
     &        - 2.d0*du(i0,i1,i2,2)*qxz_ij - l_inv*rhs_data(i0,i1,i2,0)
            r_data(i0,i1,i2,1) =
     &        c_data(i0,i1,i2,1)
     &        - 2.d0*dv(i0,i1,i2,0)*qxy_ij - 2.d0*dv(i0,i1,i2,1)*qyy_ij
     &        - 2.d0*dv(i0,i1,i2,2)*qyz_ij - l_inv*rhs_data(i0,i1,i2,1)
            r_data(i0,i1,i2,2) =
     &        c_data(i0,i1,i2,2)
     &        - 2.d0*dw(i0,i1,i2,0)*qxz_ij - 2.d0*dw(i0,i1,i2,1)*qyz_ij
     &        - 2.d0*dw(i0,i1,i2,2)*qzz_ij - l_inv*rhs_data(i0,i1,i2,2)
            r_data(i0,i1,i2,3) =
     &        c_data(i0,i1,i2,3)
     &        - qyy_ij*dw(i0,i1,i2,1) - qzz_ij*dv(i0,i1,i2,2)
     &        + qyz_ij*du(i0,i1,i2,0) - qxz_ij*dv(i0,i1,i2,0)
     &        - qxy_ij*dw(i0,i1,i2,0) - l_inv*rhs_data(i0,i1,i2,3)
            r_data(i0,i1,i2,4) =
     &        c_data(i0,i1,i2,4)
     &        - qxx_ij*dw(i0,i1,i2,0) - qzz_ij*du(i0,i1,i2,2)
     &        + qxz_ij*dv(i0,i1,i2,1) - qyz_ij*du(i0,i1,i2,1)
     &        - qxy_ij*dw(i0,i1,i2,1) - l_inv*rhs_data(i0,i1,i2,4)
            r_data(i0,i1,i2,5) =
     &        c_data(i0,i1,i2,5)
     &        - qxx_ij*dv(i0,i1,i2,0) - qyy_ij*du(i0,i1,i2,1)
     &        + qxy_ij*dw(i0,i1,i2,2) - qxz_ij*dv(i0,i1,i2,2)
     &        - qyz_ij*du(i0,i1,i2,2) - l_inv*rhs_data(i0,i1,i2,5)
          enddo
        enddo
      enddo
      end subroutine