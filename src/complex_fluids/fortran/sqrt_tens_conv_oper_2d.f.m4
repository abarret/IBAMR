define(NDIM,2)dnl
define(REAL,`double precision')dnl
define(INTEGER,`integer')dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim2d.i)dnl
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes cell centered Oldroyd-B type Convective Operator
c     for the square root evolution equation
c
c     where u is vector valued face centered velocity
c     and tau is symmetric sqrt tensor valued cell centered
c     c_data is u.grad(tau)
c     computes grad(u) using centered differences
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine sqrt_tens_conv_u_s_oper_2d
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
c    Velocity Data
c
      INTEGER u_gcw
      REAL u_data_0(SIDE2d0(ilower,iupper,u_gcw))
      REAL u_data_1(SIDE2d1(ilower,iupper,u_gcw))
c
c    Tensor Data
c
      INTEGER s_gcw
      REAL s_data(CELL2d(ilower,iupper,s_gcw),0:2)
c
c    RHS Data
c
      INTEGER rhs_gcw
      REAL rhs_data(CELL2d(ilower,iupper,rhs_gcw),0:2)
c
c    Return Data
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
      REAL du_dx, dv_dx
      REAL du_dy, dv_dy
      REAL scale0_q, scale1_q
      REAL scale_ux, scale_uy
      REAL scale_vx, scale_vy
      REAL l_inv
      REAL det
      REAL g01
      INTEGER i,j
      REAL temp(-2:2)

      REAL du(CELL2d(ilower,iupper,0),0:1)
      REAL dv(CELL2d(ilower,iupper,0),0:1)

      call weno_grad_s_to_c(u_data_0,u_data_1,u_gcw,du,dv,0,dx,
     &               0.5d0,ilower0,iupper0,ilower1,iupper1,3,3)

      l_inv = 1.d0/lambda
      scale_vx = 1.d0/(4.d0*dx(0))
      scale_vy = 1.d0/dx(1)
      scale_ux = 1.d0/dx(0)
      scale_uy = 1.d0/(4.d0*dx(1))
      scale0_q = 1.d0/(2.d0*dx(0))
      scale1_q = 1.d0/(2.d0*dx(1))

      do i1 = ilower1, iupper1
         do i0 = ilower0, iupper0
            det = s_data(i0,i1,0)*s_data(i0,i1,1)-s_data(i0,i1,2)**2

            du_dx = scale_ux*(u_data_0(i0+1,i1)-u_data_0(i0,i1))
            du_dy = scale_uy*(u_data_0(i0+1,i1+1)+u_data_0(i0,i1+1)
     &              -u_data_0(i0+1,i1-1)-u_data_0(i0,i1-1))
            dv_dy = scale_vy*(u_data_1(i0,i1+1)-u_data_1(i0,i1))
            dv_dx = scale_vx*(u_data_1(i0+1,i1+1)+u_data_1(i0+1,i1)
     &              -u_data_1(i0-1,i1)-u_data_1(i0-1,i1+1))
!             du_dx = du(i0,i1,0)
!             du_dy = du(i0,i1,1)
!             dv_dx = dv(i0,i1,0)
!             dv_dy = dv(i0,i1,1)
            g01 =
     &       ((s_data(i0,i1,2)*du_dx-s_data(i0,i1,0)*dv_dx)
     &       +(s_data(i0,i1,1)*du_dy-s_data(i0,i1,2)*dv_dy))
     &       /(s_data(i0,i1,0)+s_data(i0,i1,1))

            r_data(i0,i1,0) = c_data(i0,i1,0)
     &       -du_dx*s_data(i0,i1,0) - du_dy*s_data(i0,i1,2)
     &       -g01*s_data(i0,i1,2) - 0.5d0*l_inv*
     &       (s_data(i0,i1,1)*rhs_data(i0,i1,0)-
     &        s_data(i0,i1,2)*rhs_data(i0,i1,2))/det
            r_data(i0,i1,1) = c_data(i0,i1,1)
     &       -dv_dx*s_data(i0,i1,2) - dv_dy*s_data(i0,i1,1)
     &       +g01*s_data(i0,i1,2) - 0.5d0*l_inv*
     &       (s_data(i0,i1,0)*rhs_data(i0,i1,1)-
     &        s_data(i0,i1,2)*rhs_data(i0,i1,2))/det
            r_data(i0,i1,2) = c_data(i0,i1,2)
     &       -s_data(i0,i1,0)*dv_dx - s_data(i0,i1,2)*dv_dy
     &       -g01*s_data(i0,i1,1) - 0.5d0*l_inv*
     &       (s_data(i0,i1,1)*rhs_data(i0,i1,2)-
     &        s_data(i0,i1,2)*rhs_data(i0,i1,1))/det
         enddo
      enddo
      end subroutine
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes cell centered Oldroyd-B type Convective Operator
c     for the square root evolution
c
c     where u is vector valued cell centered velocity
c     and tau is symmetric sqrt tensor valued cell centered
c     c_data is convective term u.grad(tau)
c     computes grad(u) using an adaptive stencil
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine sqrt_tens_conv_u_c_oper_2d
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
c    Velocity Data
c
      INTEGER u_gcw
      REAL u_data(CELL2d(ilower,iupper,u_gcw),0:1)
c
c    Tensor Data
c
      INTEGER s_gcw
      REAL s_data(CELL2d(ilower,iupper,s_gcw),0:2)
c
c    RHS Data
c
      INTEGER rhs_gcw
      REAL rhs_data(CELL2d(ilower,iupper,rhs_gcw),0:2)
c
c    Return Data
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
      REAL det
      REAL g01
      INTEGER i,j
      REAL temp(-2:2)

      l_inv = 1.d0/lambda
      scale0 = 1.d0/(2.d0*dx(0))
      scale1 = 1.d0/(2.d0*dx(1))

      call compute_grad_adapt_2d(dx,
     &   u_data(:,:,0), u_gcw, du, 0,
     &   ilower0, iupper0, ilower1, iupper1, 3)
      call compute_grad_adapt_2d(dx,
     &   u_data(:,:,1), u_gcw, dv, 0,
     &   ilower0, iupper0, ilower1, iupper1, 3)

      do i1 = ilower1, iupper1
         do i0 = ilower0, iupper0
            det = s_data(i0,i1,0)*s_data(i0,i1,1)-s_data(i0,i1,2)**2

            g01 =
     &       ((s_data(i0,i1,2)*du(i0,i1,0)+s_data(i0,i1,1)*du(i0,i1,1))
     &       -(s_data(i0,i1,0)*dv(i0,i1,0)+s_data(i0,i1,2)*dv(i0,i1,1)))
     &       /(s_data(i0,i1,0)+s_data(i0,i1,1))
            r_data(i0,i1,0) = c_data(i0,i1,0) -
     &       du(i0,i1,0)*s_data(i0,i1,0) - dv(i0,i1,1)*s_data(i0,i1,2)-
     &       g01*s_data(i0,i1,2)-0.5d0*l_inv*
     &       (s_data(i0,i1,1)*rhs_data(i0,i1,0)-
     &        s_data(i0,i1,2)*rhs_data(i0,i1,2))/det
            r_data(i0,i1,1) = c_data(i0,i1,1) -
     &       dv(i0,i1,0)*s_data(i0,i1,2) - dv(i0,i1,1)*s_data(i0,i1,1)+
     &       g01*s_data(i0,i1,2)-0.5d0*l_inv*
     &       (s_data(i0,i1,0)*rhs_data(i0,i1,1)-
     &        s_data(i0,i1,2)*rhs_data(i0,i1,2))/det
            r_data(i0,i1,2) = c_data(i0,i1,2) -
     &       s_data(i0,i1,0)*dv(i0,i1,0) - s_data(i0,i1,2)*dv(i0,i1,1)-
     &       g01*s_data(i0,i1,1)-0.5d0*l_inv*
     &       (s_data(i0,i1,1)*rhs_data(i0,i1,2)-
     &        s_data(i0,i1,2)*rhs_data(i0,i1,1))/det
         enddo
      enddo
      end subroutine