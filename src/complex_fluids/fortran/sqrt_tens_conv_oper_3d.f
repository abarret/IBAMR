c234567
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
      subroutine sqrt_tens_conv_u_f_oper_3d
     &        (dx, u_data_0, u_data_1, u_data_2,
     &        u_gcw, s_data, s_gcw, rhs_data, rhs_gcw,
     &        c_data, c_gcw, r_data, r_gcw,
     &        ilower0, iupper0, ilower1,
     &        iupper1, ilower2, iupper2, lambda)

      implicit none
      Integer ilower0, iupper0
      Integer ilower1, iupper1
      Integer ilower2, iupper2

      real*8 lambda
      real*8 dx(0:2)

c
c    Velocity Data
c
      Integer u_gcw
      real*8 u_data_0((ilower0-u_gcw):(iupper0+u_gcw+1),
     &                (ilower1-u_gcw):(iupper1+u_gcw),
     &                (ilower2-u_gcw):(iupper2+u_gcw))
      real*8 u_data_1((ilower1-u_gcw):(iupper1+u_gcw+1),
     &                (ilower2-u_gcw):(iupper2+u_gcw),
     &                (ilower0-u_gcw):(iupper0+u_gcw))
      real*8 u_data_2((ilower2-u_gcw):(iupper2+u_gcw+1),
     &                (ilower0-u_gcw):(iupper0+u_gcw),
     &                (ilower1-u_gcw):(iupper1+u_gcw))
c
c    Tensor Data
c
      Integer s_gcw
      real*8 s_data((ilower0-s_gcw):(iupper0+s_gcw),
     &             (ilower1-s_gcw):(iupper1+s_gcw),
     &             (ilower2-s_gcw):(iupper2+s_gcw),
     &              0:5)
c
c    RHS Data
c
      Integer rhs_gcw
      real*8 rhs_data((ilower0-rhs_gcw):(iupper0+rhs_gcw),
     &                (ilower1-rhs_gcw):(iupper1+rhs_gcw),
     &                (ilower2-rhs_gcw):(iupper2+rhs_gcw),
     &                 0:5)
c
c    Convec Data
c
      Integer c_gcw
      real*8 c_data((ilower0-c_gcw):(iupper0+c_gcw),
     &             (ilower1-c_gcw):(iupper1+c_gcw),
     &             (ilower2-c_gcw):(iupper2+c_gcw),
     &              0:5)
c
c    Return Data
c
      Integer r_gcw
      real*8 r_data((ilower0-r_gcw):(iupper0+r_gcw),
     &             (ilower1-r_gcw):(iupper1+r_gcw),
     &             (ilower2-r_gcw):(iupper2+r_gcw),
     &              0:5)

      Integer i0, i1, i2
      real*8 du_dx, dv_dx, dw_dx
      real*8 du_dy, dv_dy, dw_dy
      real*8 du_dz, dv_dz, dw_dz
      real*8 scale_ux, scale_uy, scale_uz
      real*8 scale_vx, scale_vy, scale_vz
      real*8 scale_wx, scale_wy, scale_wz
      real*8 l_inv
      real*8 qxx_ij, qyy_ij, qzz_ij
      real*8 qyz_ij, qxz_ij, qxy_ij
      real*8 det
      real*8 g12, g23, g13

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

            det = qxx_ij*qyy_ij*qzz_ij+
     &       2.d0*qxy_ij*qxz_ij*qyz_ij-
     &       qxz_ij**2*qyy_ij-
     &       qyz_ij**2*qxx_ij-
     &       qxy_ij**2*qzz_ij

            call find_G_vals(qxx_ij, qyy_ij, qzz_ij,
     &                       qyz_ij, qxz_ij, qxy_ij,
     &                       du_dx, du_dy, du_dz,
     &                       dv_dx, dv_dy, dv_dz,
     &                       dw_dx, dw_dy, dw_dz,
     &                       g12, g13, g23)

            r_data(i0,i1,i2,0) = 
     &        c_data(i0,i1,i2,0)
     &        - qxx_ij*du_dx - qxy_ij*du_dy - qxz_ij*du_dz
     &        - g12*qxy_ij - g13*qxz_ij
     &        -0.5d0*l_inv/det*(rhs_data(i0,i1,i2,0)*(
     &        qyy_ij*qzz_ij-qyz_ij**2)+rhs_data(i0,i1,i2,4)*(
     &        qxy_ij*qyz_ij-qxz_ij*qyy_ij)+
     &        rhs_data(i0,i1,i2,5)*(qxz_ij*qyz_ij-qxy_ij*qzz_ij))
            r_data(i0,i1,i2,1) = 
     &        c_data(i0,i1,i2,1)
     &        - qxy_ij*dv_dx - qyy_ij*dv_dy - qyz_ij*dv_dz
     &        + g12*qxy_ij - g23*qyz_ij
     &        -0.5d0*l_inv/det*(rhs_data(i0,i1,i2,1)*(
     &        qxx_ij*qzz_ij-qxz_ij**2)+rhs_data(i0,i1,i2,3)*(
     &        qxy_ij*qxz_ij-qxx_ij*qyz_ij)+
     &        rhs_data(i0,i1,i2,5)*(qxz_ij*qyz_ij-qxy_ij*qzz_ij))
            r_data(i0,i1,i2,2) = 
     &        c_data(i0,i1,i2,2)
     &        - qxz_ij*dw_dx - qyz_ij*dw_dy - qzz_ij*dw_dz
     &        + g13*qxz_ij + g23*qyz_ij
     &        -0.5d0*l_inv/det*(rhs_data(i0,i1,i2,2)*(
     &        qxx_ij*qyy_ij-qxy_ij**2)+rhs_data(i0,i1,i2,3)*(
     &        qxy_ij*qxz_ij-qxx_ij*qyz_ij)+
     &        rhs_data(i0,i1,i2,4)*(qxy_ij*qyz_ij-qxz_ij*qyy_ij))
            r_data(i0,i1,i2,3) = 
     &        c_data(i0,i1,i2,3)
     &        - qxy_ij*dw_dx - qyy_ij*dw_dy - qyz_ij*dw_dz
     &        + g12*qxz_ij - g23*qzz_ij
     &        -0.5d0*l_inv/det*(rhs_data(i0,i1,i2,2)*(
     &        qxy_ij*qxz_ij-qxx_ij*qyz_ij)+rhs_data(i0,i1,i2,3)
     &        *(qxx_ij*qzz_ij-qxz_ij**2)+rhs_data(i0,i1,i2,4)*(
     &        qxz_ij*qyz_ij-qxy_ij*qzz_ij))
            r_data(i0,i1,i2,4) = 
     &        c_data(i0,i1,i2,4)
     &        - qxx_ij*dw_dx - qxy_ij*dw_dy - qxz_ij*dw_dz
     &        - g12*qyz_ij - g13*qzz_ij
     &        -0.5d0*l_inv/det*(rhs_data(i0,i1,i2,2)*(
     &        qxy_ij*qyz_ij-qxz_ij*qyy_ij)+rhs_data(i0,i1,i2,3)
     &        *(qxz_ij*qyz_ij-qxy_ij*qzz_ij)+rhs_data(i0,i1,i2,4)
     &        *(qyy_ij*qzz_ij-qyz_ij**2))
            r_data(i0,i1,i2,5) = 
     &        c_data(i0,i1,i2,5)
     &        - qxx_ij*dv_dx - qxy_ij*dv_dy - qxz_ij*dv_dz
     &        - g12*qyy_ij - g13*qyz_ij
     &        -0.5d0*l_inv/det*(rhs_data(i0,i1,i2,1)*(
     &        qxz_ij*qyz_ij-qxy_ij*qzz_ij)+rhs_data(i0,i1,i2,3)
     &        *(qxy_ij*qyz_ij-qxz_ij*qyy_ij)+rhs_data(i0,i1,i2,5)
     &        *(qyy_ij*qzz_ij-qyz_ij**2))
          enddo
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
      subroutine sqrt_tens_conv_u_c_oper_3d
     &        (dx, u_data,
     &        u_gcw, s_data, s_gcw, rhs_data, rhs_gcw,
     &        c_data, c_gcw, r_data, r_gcw,
     &        ilower0, iupper0, ilower1,
     &        iupper1, ilower2, iupper2, lambda)

      implicit none
      Integer ilower0, iupper0
      Integer ilower1, iupper1
      Integer ilower2, iupper2

      real*8 lambda
      real*8 dx(0:2)

c
c    Velocity Data
c
      Integer u_gcw
      real*8 u_data((ilower0-u_gcw):(iupper0+u_gcw),
     &              (ilower1-u_gcw):(iupper1+u_gcw),
     &              (ilower2-u_gcw):(iupper2+u_gcw),
     &               0:2)
c
c    Tensor Data
c
      Integer s_gcw
      real*8 s_data((ilower0-s_gcw):(iupper0+s_gcw),
     &             (ilower1-s_gcw):(iupper1+s_gcw),
     &             (ilower2-s_gcw):(iupper2+s_gcw),
     &              0:5)
c
c    RHS Data
c
      Integer rhs_gcw
      real*8 rhs_data((ilower0-rhs_gcw):(iupper0+rhs_gcw),
     &                (ilower1-rhs_gcw):(iupper1+rhs_gcw),
     &                (ilower2-rhs_gcw):(iupper2+rhs_gcw),
     &                 0:5)
c
c    Convec Data
c
      Integer c_gcw
      real*8 c_data((ilower0-c_gcw):(iupper0+c_gcw),
     &             (ilower1-c_gcw):(iupper1+c_gcw),
     &             (ilower2-c_gcw):(iupper2+c_gcw),
     &              0:5)
c
c    Return Data
c
      Integer r_gcw
      real*8 r_data((ilower0-r_gcw):(iupper0+r_gcw),
     &             (ilower1-r_gcw):(iupper1+r_gcw),
     &             (ilower2-r_gcw):(iupper2+r_gcw),
     &              0:5)

      Integer i0, i1, i2
      real*8 du(ilower0:iupper0,
     &          ilower1:iupper1,
     &          ilower2:iupper2,
     &          0:2)
      real*8 dv(ilower0:iupper0,
     &          ilower1:iupper1,
     &          ilower2:iupper2,
     &          0:2)
      real*8 dw(ilower0:iupper0,
     &          ilower1:iupper1,
     &          ilower2:iupper2,
     &          0:2)
      real*8 scale0, scale1, scale2
      real*8 l_inv
      real*8 qxx_ij, qyy_ij, qzz_ij
      real*8 qyz_ij, qxz_ij, qxy_ij
      real*8 det
      real*8 g12, g23, g13

      l_inv = 1.d0/lambda
      scale0 = 1.d0/(2.d0*dx(0))
      scale1 = 1.d0/(2.d0*dx(1))
      scale2 = 1.d0/(2.d0*dx(2))

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

            det = qxx_ij*qyy_ij*qzz_ij+
     &       2.d0*qxy_ij*qxz_ij*qyz_ij-
     &       qxz_ij**2*qyy_ij-
     &       qyz_ij**2*qxx_ij-
     &       qxy_ij**2*qzz_ij
     
!          print *, dw(i0,i1,i2,0), dw(i0,i1,i2,1), dw(i0,i1,i2,2)

            call find_G_vals(qxx_ij, qyy_ij, qzz_ij,
     &             qyz_ij, qxz_ij, qxy_ij,
     &             du(i0,i1,i2,0), du(i0,i1,i2,1), du(i0,i1,i2,2),
     &             dv(i0,i1,i2,0), dv(i0,i1,i2,1), dv(i0,i1,i2,2),
     &             dw(i0,i1,i2,0), dw(i0,i1,i2,1), dw(i0,i1,i2,2),
     &             g12, g13, g23)

            r_data(i0,i1,i2,0) = 
     &        c_data(i0,i1,i2,0)
     &        - qxx_ij*du(i0,i1,i2,0) - qxy_ij*du(i0,i1,i2,1) - 
     &          qxz_ij*du(i0,i1,i2,2)
     &        - g12*qxy_ij - g13*qxz_ij
     &        -0.5d0*l_inv/det*(rhs_data(i0,i1,i2,0)*(
     &        qyy_ij*qzz_ij-qyz_ij**2)+rhs_data(i0,i1,i2,4)*(
     &        qxy_ij*qyz_ij-qxz_ij*qyy_ij)+
     &        rhs_data(i0,i1,i2,5)*(qxz_ij*qyz_ij-qxy_ij*qzz_ij))
            r_data(i0,i1,i2,1) = 
     &        c_data(i0,i1,i2,1)
     &        - qxy_ij*dv(i0,i1,i2,0) - qyy_ij*dv(i0,i1,i2,1) - 
     &          qyz_ij*dv(i0,i1,i2,2)
     &        + g12*qxy_ij - g23*qyz_ij
     &        -0.5d0*l_inv/det*(rhs_data(i0,i1,i2,1)*(
     &        qxx_ij*qzz_ij-qxz_ij**2)+rhs_data(i0,i1,i2,3)*(
     &        qxy_ij*qxz_ij-qxx_ij*qyz_ij)+
     &        rhs_data(i0,i1,i2,5)*(qxz_ij*qyz_ij-qxy_ij*qzz_ij))
            r_data(i0,i1,i2,2) = 
     &        c_data(i0,i1,i2,2)
     &        - qxz_ij*dw(i0,i1,i2,0) - qyz_ij*dw(i0,i1,i2,1) - 
     &          qzz_ij*dw(i0,i1,i2,2)
     &        + g13*qxz_ij + g23*qyz_ij
     &        -0.5d0*l_inv/det*(rhs_data(i0,i1,i2,2)*(
     &        qxx_ij*qyy_ij-qxy_ij**2)+rhs_data(i0,i1,i2,3)*(
     &        qxy_ij*qxz_ij-qxx_ij*qyz_ij)+
     &        rhs_data(i0,i1,i2,4)*(qxy_ij*qyz_ij-qxz_ij*qyy_ij))
            r_data(i0,i1,i2,3) = 
     &        c_data(i0,i1,i2,3)
     &        - qxy_ij*dw(i0,i1,i2,0) - qyy_ij*dw(i0,i1,i2,1) - 
     &          qyz_ij*dw(i0,i1,i2,2)
     &        + g12*qxz_ij - g23*qzz_ij
     &        -0.5d0*l_inv/det*(rhs_data(i0,i1,i2,2)*(
     &        qxy_ij*qxz_ij-qxx_ij*qyz_ij)+rhs_data(i0,i1,i2,3)
     &        *(qxx_ij*qzz_ij-qxz_ij**2)+rhs_data(i0,i1,i2,4)*(
     &        qxz_ij*qyz_ij-qxy_ij*qzz_ij))
            r_data(i0,i1,i2,4) = 
     &        c_data(i0,i1,i2,4)
     &        - qxx_ij*dw(i0,i1,i2,0) - qxy_ij*dw(i0,i1,i2,1) -
     &          qxz_ij*dw(i0,i1,i2,2)
     &        - g12*qyz_ij - g13*qzz_ij
     &        -0.5d0*l_inv/det*(rhs_data(i0,i1,i2,2)*(
     &        qxy_ij*qyz_ij-qxz_ij*qyy_ij)+rhs_data(i0,i1,i2,3)
     &        *(qxz_ij*qyz_ij-qxy_ij*qzz_ij)+rhs_data(i0,i1,i2,4)
     &        *(qyy_ij*qzz_ij-qyz_ij**2))
            r_data(i0,i1,i2,5) = 
     &        c_data(i0,i1,i2,5)
     &        - qxx_ij*dv(i0,i1,i2,0) - qxy_ij*dv(i0,i1,i2,1) -
     &          qxz_ij*dv(i0,i1,i2,2)
     &        - g12*qyy_ij - g13*qyz_ij
     &        -0.5d0*l_inv/det*(rhs_data(i0,i1,i2,1)*(
     &        qxz_ij*qyz_ij-qxy_ij*qzz_ij)+rhs_data(i0,i1,i2,3)
     &        *(qxy_ij*qyz_ij-qxz_ij*qyy_ij)+rhs_data(i0,i1,i2,5)
     &        *(qyy_ij*qzz_ij-qyz_ij**2))
          enddo
        enddo
      enddo
      end subroutine

      subroutine find_G_vals(q11,q22,q33,
     &           q23, q13, q12,
     &           dudx, dudy, dudz,
     &           dvdx, dvdy, dvdz,
     &           dwdx, dwdy, dwdz,
     &           g12, g13, g23)
      implicit none
c     INPUTS
      real*8 q11, q22, q33
      real*8 q23, q13, q12
      real*8 dudx, dudy, dudz
      real*8 dvdx, dvdy, dvdz
      real*8 dwdx, dwdy, dwdz
c     TO RETURN
      real*8 g12, g13, g23
c     SCRATCH
      real*8 A(3,3), x(3), b(3)
      Integer ipiv(3), info

      A(1,1) = q11 + q22; A(1,2) = q23; A(1,3) = -q13
      A(2,1) = q23; A(2,2) = q11 + q33; A(2,3) = q12
      A(3,1) = -q13; A(3,2) = q12; A(3,3) = q22+q33

      b(1) = q12*dudx-q11*dvdx+q22*dudy-q12*dvdy+q23*dudz-q13*dvdz
      b(2) = q13*dudx-q11*dwdx+q33*dudz-q13*dwdz+q23*dudy-q12*dwdy
      b(3) = q13*dvdx-q12*dwdx+q23*dvdy-q22*dwdy+q33*dvdz-q23*dwdz

      call dgesv(3, 1, A, 3, ipiv, b, 3, info)

      if (info /= 0) then
        print *, 'ERROR IN DGESV'
      endif

      g12 = b(1)
      g13 = b(2)
      g23 = b(3)

      endsubroutine
