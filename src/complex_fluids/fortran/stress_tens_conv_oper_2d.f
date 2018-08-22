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
      subroutine stress_tens_conv_oper_2d
     &        (dx, u_data_0, u_data_1, u_gcw,
     &        t_data, t_gcw, ra_data, ra_gcw,
     &        re_data, re_gcw, r_data, r_gcw,
     &        ilower0, iupper0, ilower1, iupper1,
     &        lambda, eta)

      parameter (NDIM = 2)
      Integer ilower0, iupper0
      Integer ilower1, iupper1

      real*8 lambda, eta
      real*8 dx(0:NDIM-1)

c
c    Velocity Data
c
      Integer u_gcw
      real*8 u_data_0((ilower0-u_gcw):(iupper0+u_gcw+1),
     &                (ilower1-u_gcw):(iupper1+u_gcw))
      real*8 u_data_1((ilower1-u_gcw):(iupper1+u_gcw+1),
     &                (ilower0-u_gcw):(iupper0+u_gcw))
c
c    Tensor Data
c
      Integer t_gcw
      real*8 t_data((ilower0-t_gcw):(iupper0+t_gcw),
     &              (ilower1-t_gcw):(iupper1+t_gcw),
     &              0:2)
c
c    RATE Data
c
      Integer ra_gcw
      real*8 ra_data((ilower0-ra_gcw):(iupper0+ra_gcw),
     &                (ilower1-ra_gcw):(iupper1+ra_gcw),
     &                 0:2)
c
c    RELAX Data
c
      Integer re_gcw
      real*8 re_data((ilower0-re_gcw):(iupper0+re_gcw),
     &               (ilower1-re_gcw):(iupper1+re_gcw),
     &                0:2)
c
c    Return Data
c
      Integer r_gcw
      real*8 r_data((ilower0-r_gcw):(iupper0+r_gcw),
     &              (ilower1-r_gcw):(iupper1+r_gcw),
     &               0:2)

      Integer i0, i1
      real*8 u_ij, v_ij
      real*8 dq_00_dx, dq_11_dx, dq_01_dx
      real*8 dq_00_dy, dq_11_dy, dq_01_dy
      real*8 du_dx, dv_dx
      real*8 du_dy, dv_dy
      real*8 scale0_q, scale1_q
      real*8 scale_ux, scale_uy
      real*8 scale_vx, scale_vy
      real*8 wi_inv

      wi_inv = 1.d0/wi
      scale_vx = 1.d0/(4.d0*dx(0))
      scale_vy = 1.d0/dx(1)
      scale_ux = 1.d0/dx(0)
      scale_uy = 1.d0/(4.d0*dx(1))
      scale0_q = 1.d0/(2.d0*dx(0))
      scale1_q = 1.d0/(2.d0*dx(1))
      do i1 = ilower1, iupper1
         do i0 = ilower0, iupper0
            u_ij = u_data_0(i0,i1)*0.5+u_data_0(i0+1,i1)*0.5
            v_ij = u_data_1(i1,i0)*0.5+u_data_1(i1+1,i0)*0.5
            dq_00_dx = scale0_q*(t_data(i0+1,i1,0)-t_data(i0-1,i1,0))
            dq_00_dy = scale1_q*(t_data(i0,i1+1,0)-t_data(i0,i1-1,0))
            dq_11_dx = scale0_q*(t_data(i0+1,i1,1)-t_data(i0-1,i1,1))
            dq_11_dy = scale1_q*(t_data(i0,i1+1,1)-t_data(i0,i1-1,1))
            dq_01_dx = scale0_q*(t_data(i0+1,i1,2)-t_data(i0-1,i1,2))
            dq_01_dy = scale1_q*(t_data(i0,i1+1,2)-t_data(i0,i1-1,2))
            du_dx = scale_ux*(u_data_0(i0+1,i1)-u_data_0(i0,i1))
            du_dy = scale_uy*(u_data_0(i0+1,i1+1)+u_data_0(i0,i1+1)
     &              -u_data_0(i0+1,i1-1)-u_data_0(i0,i1-1))
            dv_dy = scale_vy*(u_data_1(i1+1,i0)-u_data_1(i1,i0))
            dv_dx = scale_vx*(u_data_1(i1+1,i0+1)+u_data_1(i1,i0+1)
     &              -u_data_1(i1,i0-1) - u_data_1(i1+1,i0-1))
            r_data(i0,i1,0) = u_ij*dq_00_dx + v_ij*dq_00_dy - 
     &       2.d0*du_dx*t_data(i0,i1,0) - 2.d0*du_dy*t_data(i0,i1,2)
     &       + ra_data(i0,i1,0) + 1/lambda*(re_data(i0,i1,0)*
     &       t_data(i0,i1,0) + re_data(i0,i1,2)*t_data(i0,i1,2)) -
     &       eta/lambda*(2.0*du_dx)
            r_data(i0,i1,1) = u_ij*dq_11_dx + v_ij*dq_11_dy - 
     &       2.d0*dv_dx*t_data(i0,i1,2) - 2.d0*dv_dy*t_data(i0,i1,1)
     &       + ra_data(i0,i1,1) + 1/lambda*(re_data(i0,i1,2)*
     &       t_data(i0,i1,2) + re_data(i0,i1,1)*t_data(i0,i1,1)) - 
     &       eta/lambda*(2.0*dv_dy)
            r_data(i0,i1,2) = u_ij*dq_01_dx + v_ij*dq_01_dy - 
     &       t_data(i0,i1,1)*du_dy - t_data(i0,i1,0)*dv_dx
     &       + ra_data(i0,i1,2) + 1/lambda*(re_data(i0,i1,0)*
     &       t_data(i0,i1,2) + re_data(i0,i1,2)*t_data(i0,i1,1)) - 
     &       eta/lambda*(dv_dx+du_dy)
         enddo
      enddo
      end subroutine