c234567
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes cell centered Oldroyd-B type Convective Operator
c     for the stress tensor
c
c     where u is vector valued face centered velocity
c     and tau is symmetric tensor valued cell centered
c     computes grad(u) using centered differences
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine stress_tens_conv_oper_3d
     &        (dx, u_data_0, u_data_1, u_data_2, u_gcw,
     &         t_data, t_gcw, ra_data, ra_gcw,
     &         re_data, re_gcw, r_data, r_gcw,
     &         ilower0, iupper0, ilower1,
     &         iupper1, ilower2, iupper2,
     &         lambda, eta)
      implicit none
      Integer ilower0, iupper0
      Integer ilower1, iupper1
      Integer ilower2, iupper2

      real*8 lambda, eta
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
      Integer t_gcw
      real*8 t_data((ilower0-t_gcw):(iupper0+t_gcw),
     &              (ilower1-t_gcw):(iupper1+t_gcw),
     &              (ilower2-t_gcw):(iupper2+t_gcw),
     &               0:5)
c
c    RHS Data
c
      Integer ra_gcw
      real*8 ra_data((ilower0-ra_gcw):(iupper0+ra_gcw),
     &               (ilower1-ra_gcw):(iupper1+ra_gcw),
     &               (ilower2-ra_gcw):(iupper2+ra_gcw),
     &                0:5)

      Integer re_gcw
      real*8 re_data((ilower0-re_gcw):(iupper0+re_gcw),
     &               (ilower1-re_gcw):(iupper1+re_gcw),
     &               (ilower2-re_gcw):(iupper2+re_gcw),
     &                0:5)
c
c    Return Data
c
      Integer r_gcw
      real*8 r_data((ilower0-r_gcw):(iupper0+r_gcw),
     &              (ilower1-r_gcw):(iupper1+r_gcw),
     &              (ilower2-r_gcw):(iupper2+r_gcw),
     &               0:5)

      Integer i0, i1, i2
      real*8 u_ij, v_ij, w_ij
      real*8 dq_xx_dx, dq_yy_dx, dq_zz_dx
      real*8 dq_yz_dx, dq_xz_dx, dq_xy_dx
      real*8 dq_xx_dy, dq_yy_dy, dq_zz_dy
      real*8 dq_yz_dy, dq_xz_dy, dq_xy_dy
      real*8 dq_xx_dz, dq_yy_dz, dq_zz_dz
      real*8 dq_yz_dz, dq_xz_dz, dq_xy_dz
      real*8 du_dx, dv_dx, dw_dx
      real*8 du_dy, dv_dy, dw_dy
      real*8 du_dz, dv_dz, dw_dz
      real*8 scalex_q, scaley_q, scalez_q
      real*8 scale_ux, scale_uy, scale_uz
      real*8 scale_vx, scale_vy, scale_vz
      real*8 scale_wx, scale_wy, scale_wz
      real*8 wi_inv
      real*8 qxx_ij, qyy_ij, qzz_ij
      real*8 qyz_ij, qxz_ij, qxy_ij

      wi_inv = 1.d0/lambda
      scale_ux = 1.d0/dx(0)
      scale_uy = 1.d0/(4.d0*dx(1))
      scale_uz = 1.d0/(4.d0*dx(2))
      scale_vx = 1.d0/(4.d0*dx(0))
      scale_vy = 1.d0/dx(1)
      scale_vz = 1.d0/(4.d0*dx(2))
      scale_wx = 1.d0/(4.d0*dx(0))
      scale_wy = 1.d0/(4.d0*dx(1))
      scale_wz = 1.d0/dx(2)
      scalex_q = 1.d0/(2.d0*dx(0))
      scaley_q = 1.d0/(2.d0*dx(1))
      scalez_q = 1.d0/(2.d0*dx(2))
      do i2 = ilower2, iupper2
        do i1 = ilower1, iupper1
          do i0 = ilower0, iupper0
            u_ij = (u_data_0(i0,i1,i2)+u_data_0(i0+1,i1,i2))*0.5
            du_dx = 
     &        scale_ux*(u_data_0(i0+1,i1,i2)-u_data_0(i0,i1,i2))
            du_dy = 
     &        scale_uy*(u_data_0(i0+1,i1+1,i2)+u_data_0(i0,i1+1,i2)
     &              -u_data_0(i0+1,i1-1,i2)-u_data_0(i0,i1-1,i2))
            du_dz = 
     &        scale_uz*(u_data_0(i0+1,i1,i2+1)+u_data_0(i0,i1,i2+1)
     &              -u_data_0(i0+1,i1,i2-1)-u_data_0(i0,i1,i2-1))
            v_ij = (u_data_1(i1,i2,i0)+u_data_1(i1+1,i2,i0))*0.5
            dv_dy = 
     &        scale_vy*(u_data_1(i1+1,i2,i0)-u_data_1(i1,i2,i0))
            dv_dx = 
     &        scale_vx*(u_data_1(i1+1,i2,i0+1)+u_data_1(i1,i2,i0+1)
     &              -u_data_1(i1+1,i2,i0-1) - u_data_1(i1,i2,i0-1))
            dv_dz = 
     &        scale_vz*(u_data_1(i1+1,i2+1,i0)+u_data_1(i1,i2+1,i0)
     &              -u_data_1(i1+1,i2-1,i0) - u_data_1(i1,i2-1,i0))
            w_ij = (u_data_2(i2,i0,i1)+u_data_2(i2+1,i0,i1))*0.5
            dw_dx = 
     &        scale_wx*(u_data_2(i2+1,i0+1,i1)+u_data_2(i2,i0+1,i1)
     &              -u_data_2(i2+1,i0-1,i1)-u_data_2(i2,i0-1,i1))
            dw_dy = 
     &        scale_wy*(u_data_2(i2+1,i0,i1+1)+u_data_2(i2,i0,i1+1)
     &              -u_data_2(i2+1,i0,i1-1)-u_data_2(i2,i0,i1-1))
            dw_dz = 
     &        scale_wz*(u_data_2(i2+1,i0,i1)-u_data_2(i2,i0,i1))
            dq_xx_dx = scalex_q*(t_data(i0+1,i1,i2,0)
     &                 -t_data(i0-1,i1,i2,0))
            dq_xx_dy = scaley_q*(t_data(i0,i1+1,i2,0)
     &                 -t_data(i0,i1-1,i2,0))
            dq_xx_dz = scalez_q*(t_data(i0,i1,i2+1,0)
     &                 -t_data(i0,i1,i2-1,0))
            qxx_ij = t_data(i0,i1,i2,0)
            dq_yy_dx = scalex_q*(t_data(i0+1,i1,i2,1)
     &                 -t_data(i0-1,i1,i2,1))
            dq_yy_dy = scaley_q*(t_data(i0,i1+1,i2,1)
     &                 -t_data(i0,i1-1,i2,1))
            dq_yy_dz = scalez_q*(t_data(i0,i1,i2+1,1)
     &                 -t_data(i0,i1,i2-1,1))
            qyy_ij = t_data(i0,i1,i2,1)
            dq_zz_dx = scalex_q*(t_data(i0+1,i1,i2,2)
     &                 -t_data(i0-1,i1,i2,2))
            dq_zz_dy = scaley_q*(t_data(i0,i1+1,i2,2)
     &                 -t_data(i0,i1-1,i2,2))
            dq_zz_dz = scalez_q*(t_data(i0,i1,i2+1,2)
     &                 -t_data(i0,i1,i2-1,2))
            qzz_ij = t_data(i0,i1,i2,2)
            dq_yz_dx = scalex_q*(t_data(i0+1,i1,i2,3)
     &                 -t_data(i0-1,i1,i2,3))
            dq_yz_dy = scaley_q*(t_data(i0,i1+1,i2,3)
     &                 -t_data(i0,i1-1,i2,3))
            dq_yz_dz = scalez_q*(t_data(i0,i1,i2+1,3)
     &                 -t_data(i0,i1,i2-1,3))
            qyz_ij = t_data(i0,i1,i2,3)
            dq_xz_dx = scalex_q*(t_data(i0+1,i1,i2,4)
     &                 -t_data(i0-1,i1,i2,4))
            dq_xz_dy = scaley_q*(t_data(i0,i1+1,i2,4)
     &                 -t_data(i0,i1-1,i2,4))
            dq_xz_dz = scalez_q*(t_data(i0,i1,i2+1,4)
     &                 -t_data(i0,i1,i2-1,4))
            qxz_ij = t_data(i0,i1,i2,4)
            dq_xy_dx = scalex_q*(t_data(i0+1,i1,i2,5)
     &                 -t_data(i0-1,i1,i2,5))
            dq_xy_dy = scaley_q*(t_data(i0,i1+1,i2,5)
     &                 -t_data(i0,i1-1,i2,5))
            dq_xy_dz = scalez_q*(t_data(i0,i1,i2+1,5)
     &                 -t_data(i0,i1,i2-1,5))
            qxy_ij = t_data(i0,i1,i2,5)

            r_data(i0,i1,i2,0) = 
     &        u_ij*dq_xx_dx + v_ij*dq_xx_dy + w_ij*dq_xx_dz
     &        - 2.d0*du_dx*qxx_ij - 2.d0*du_dy*qxy_ij
     &        - 2.d0*du_dz*qxz_ij + ra_data(i0,i1,i2,0)
     &        + 1/lambda*(re_data(i0,i1,i2,0)*qxx_ij
     &        + re_data(i0,i1,i2,5)*qxy_ij
     &        + re_data(i0,i1,i2,4)*qxz_ij)
     &        - eta/lambda*(2.0*du_dx)
            r_data(i0,i1,i2,1) = 
     &        u_ij*dq_yy_dx + v_ij*dq_yy_dy + w_ij*dq_yy_dz 
     &        - 2.d0*dv_dx*qxy_ij - 2.d0*dv_dy*qyy_ij
     &        - 2.d0*dv_dz*qyz_ij + ra_data(i0,i1,i2,1)
     &        + 1/lambda*(re_data(i0,i1,i2,5)*qxy_ij
     &        + re_data(i0,i1,i2,1)*qyy_ij
     &        + re_data(i0,i1,i2,3)*qyz_ij)
     &        - eta/lambda*(2.0*dv_dy)
            r_data(i0,i1,i2,2) = 
     &        u_ij*dq_zz_dx + v_ij*dq_zz_dy + w_ij*dq_zz_dz
     &        - 2.d0*dw_dx*qxz_ij - 2.d0*dw_dy*qyz_ij
     &        - 2.d0*dw_dz*qzz_ij + ra_data(i0,i1,i2,2)
     &        + 1/lambda*(re_data(i0,i1,i2,4)*qxz_ij
     &        + re_data(i0,i1,i2,3)*qyz_ij
     &        + re_data(i0,i1,i2,2)*qzz_ij)
     &        - eta/lambda*(2.0*dw_dz)
            r_data(i0,i1,i2,3) = 
     &        u_ij*dq_yz_dx + v_ij*dq_yz_dy + w_ij*dq_yz_dz
     &        - qyy_ij*dw_dy - qzz_ij*dv_dz
     &        + qyz_ij*du_dx - qxz_ij*dv_dx
     &        - qxy_ij*dw_dx + ra_data(i0,i1,i2,3)
     &        + 1/lambda*(re_data(i0,i1,i2,5)*qxz_ij
     &        + re_data(i0,i1,i2,1)*qyz_ij
     &        + re_data(i0,i1,i2,3)*qzz_ij)
     &        - eta/lambda*(dw_dy+dv_dz)
            r_data(i0,i1,i2,4) = 
     &        u_ij*dq_xz_dx + v_ij*dq_xz_dy + w_ij*dq_xz_dz
     &        - qxx_ij*dw_dx - qzz_ij*du_dz
     &        + qxz_ij*dv_dy - qyz_ij*du_dy
     &        - qxy_ij*dw_dy + ra_data(i0,i1,i2,4)
     &        + 1/lambda*(re_data(i0,i1,i2,0)*qxz_ij
     &        + re_data(i0,i1,i2,5)*qyz_ij
     &        + re_data(i0,i1,i2,4)*qzz_ij)
     &        - eta/lambda*(dw_dx+du_dz)
            r_data(i0,i1,i2,5) = 
     &        u_ij*dq_xy_dx + v_ij*dq_xy_dy + w_ij*dq_xy_dz
     &        - qxx_ij*dv_dx - qyy_ij*du_dy
     &        + qxy_ij*dw_dz - qxz_ij*dv_dz
     &        - qyz_ij*du_dz + ra_data(i0,i1,i2,5)
     &        + 1/lambda*(re_data(i0,i1,i2,0)*qxy_ij
     &        + re_data(i0,i1,i2,5)*qyy_ij
     &        + re_data(i0,i1,i2,4)*qyz_ij)
     &        - eta/lambda*(dv_dx+du_dy)
          enddo
        enddo
      enddo
      end subroutine
