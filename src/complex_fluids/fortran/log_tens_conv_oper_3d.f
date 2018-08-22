c234567
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes cell centered Oldroyd-B type Convective Operator
c     for the log-evolution equation
c
c     where u is vector valued face centered velocity
c     and tau is symmetric log tensor valued cell centered
c     c_data is convective terms u.grad(tau)
c     computes grad(u) using centered differences
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine log_tens_conv_u_f_oper_3d
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
      real*8 qxx, qyy, qzz
      real*8 qyz, qxz, qxy

      real*8 L(3,3), vecs(3,3)
      real*8 vals(3)
      real*8 convec_vals(3,3)
      real*8 sigma(0:5)
      real*8 temp(15)

      Integer info

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

            qxx = s_data(i0,i1,i2,0)
            qyy = s_data(i0,i1,i2,1)
            qzz = s_data(i0,i1,i2,2)
            qyz = s_data(i0,i1,i2,3)
            qxz = s_data(i0,i1,i2,4)
            qxy = s_data(i0,i1,i2,5)

            L(1,1) = du_dx; L(1,2) = du_dy; L(1,3) = du_dz
            L(2,1) = dv_dx; L(2,2) = dv_dy; L(2,3) = dv_dz
            L(3,1) = dw_dx; L(3,2) = dw_dy; L(3,3) = dw_dz

            vecs(1,1) = qxx; vecs(1,2) = qxy; vecs(1,3) = qxz
            vecs(2,1) = qxy; vecs(2,2) = qyy; vecs(2,3) = qyz
            vecs(3,1) = qxz; vecs(3,2) = qyz; vecs(3,3) = qzz

            call dsyev('V','U',3,vecs,3,vals,temp,15,info)
            if (info /= 0) then
              print *, "ERROR IN DSYEV!!"
            endif
            call log_sum(L,vecs,vals,3,convec_vals)

            sigma(0) = vecs(1,1)**2*exp(-vals(1))
     &          +vecs(1,2)**2*exp(-vals(2))
     &          +vecs(1,3)**2*exp(-vals(3))
            sigma(1) = vecs(2,1)**2*exp(-vals(1))
     &          +vecs(2,2)**2*exp(-vals(2))
     &          +vecs(2,3)**2*exp(-vals(3))
            sigma(2) = vecs(3,1)**2*exp(-vals(1))
     &          +vecs(3,2)**2*exp(-vals(2))
     &          +vecs(3,3)**2*exp(-vals(3))
            sigma(3) = vecs(2,1)*vecs(3,1)*exp(-vals(1))
     &          +vecs(2,2)*vecs(3,2)*exp(-vals(2))
     &          +vecs(2,3)*vecs(3,3)*exp(-vals(3))
            sigma(4) = vecs(1,1)*vecs(3,1)*exp(-vals(1))
     &          +vecs(1,2)*vecs(3,2)*exp(-vals(2))
     &          +vecs(1,3)*vecs(3,3)*exp(-vals(3))
            sigma(5) = vecs(1,1)*vecs(2,1)*exp(-vals(1))
     &          +vecs(1,2)*vecs(2,2)*exp(-vals(2))
     &          +vecs(1,3)*vecs(2,3)*exp(-vals(3))

            r_data(i0,i1,i2,0) = 
     &        c_data(i0,i1,i2,0)
     &        - convec_vals(1,1) - l_inv*(
     &        sigma(0)*rhs_data(i0,i1,i2,0) +
     &        sigma(4)*rhs_data(i0,i1,i2,4) +
     &        sigma(5)*rhs_data(i0,i1,i2,5))

            r_data(i0,i1,i2,1) =
     &        c_data(i0,i1,i2,1)
     &        - convec_vals(2,2) - l_inv*(
     &        sigma(1)*rhs_data(i0,i1,i2,1) +
     &        sigma(3)*rhs_data(i0,i1,i2,3) +
     &        sigma(5)*rhs_data(i0,i1,i2,5))

            r_data(i0,i1,i2,2) =
     &        c_data(i0,i1,i2,2)
     &        - convec_vals(3,3) - l_inv*(
     &        sigma(2)*rhs_data(i0,i1,i2,2) +
     &        sigma(3)*rhs_data(i0,i1,i2,3) +
     &        sigma(4)*rhs_data(i0,i1,i2,4))

            r_data(i0,i1,i2,3) =
     &        c_data(i0,i1,i2,3)
     &        - convec_vals(2,3) - l_inv*(
     &        sigma(1)*rhs_data(i0,i1,i2,3) +
     &        sigma(3)*rhs_data(i0,i1,i2,2) +
     &        sigma(5)*rhs_data(i0,i1,i2,4))

            r_data(i0,i1,i2,4) =
     &        c_data(i0,i1,i2,4)
     &        - convec_vals(1,3) - l_inv*(
     &        sigma(0)*rhs_data(i0,i1,i2,4) +
     &        sigma(4)*rhs_data(i0,i1,i2,2) +
     &        sigma(5)*rhs_data(i0,i1,i2,3))

            r_data(i0,i1,i2,5) =
     &        c_data(i0,i1,i2,5)
     &        - convec_vals(1,2) - l_inv*(
     &        sigma(0)*rhs_data(i0,i1,i2,5) +
     &        sigma(4)*rhs_data(i0,i1,i2,3) +
     &        sigma(5)*rhs_data(i0,i1,i2,1))
          enddo
        enddo
      enddo
      end subroutine
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes cell centered Oldroyd-B type Convective Operator
c     for the log-evolution equation
c
c     where u is vector valued cell centered velocity
c     and tau is symmetric log tensor valued cell centered
c     c_data is convective term u.grad(tau)
c     computes grad(u) using an adaptive stencil
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine log_tens_conv_u_c_oper_3d
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
      real*8 qxx, qyy, qzz
      real*8 qyz, qxz, qxy

      real*8 L(3,3), vecs(3,3)
      real*8 vals(3)
      real*8 convec_vals(3,3)
      real*8 sigma(0:5)
      real*8 temp(15)

      Integer info

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
            qxx = s_data(i0,i1,i2,0)
            qyy = s_data(i0,i1,i2,1)
            qzz = s_data(i0,i1,i2,2)
            qyz = s_data(i0,i1,i2,3)
            qxz = s_data(i0,i1,i2,4)
            qxy = s_data(i0,i1,i2,5)

            L(1,1) = du(i0,i1,i2,0); L(1,2) = du(i0,i1,i2,1)
                     L(1,3) = du(i0,i1,i2,2)
            L(2,1) = dv(i0,i1,i2,0); L(2,2) = dv(i0,i1,i2,1)
                     L(2,3) = dv(i0,i1,i2,2)
            L(3,1) = dw(i0,i1,i2,0); L(3,2) = dw(i0,i1,i2,1)
                     L(3,3) = dw(i0,i1,i2,2)

            vecs(1,1) = qxx; vecs(1,2) = qxy; vecs(1,3) = qxz
            vecs(2,1) = qxy; vecs(2,2) = qyy; vecs(2,3) = qyz
            vecs(3,1) = qxz; vecs(3,2) = qyz; vecs(3,3) = qzz

            call dsyev('V','U',3,vecs,3,vals,temp,15,info)
            if (info /= 0) then
              print *, "ERROR IN DSYEV!!"
            endif
            call log_sum(L,vecs,vals,3,convec_vals)

            sigma(0) = vecs(1,1)**2*exp(-vals(1))
     &          +vecs(1,2)**2*exp(-vals(2))
     &          +vecs(1,3)**2*exp(-vals(3))
            sigma(1) = vecs(2,1)**2*exp(-vals(1))
     &          +vecs(2,2)**2*exp(-vals(2))
     &          +vecs(2,3)**2*exp(-vals(3))
            sigma(2) = vecs(3,1)**2*exp(-vals(1))
     &          +vecs(3,2)**2*exp(-vals(2))
     &          +vecs(3,3)**2*exp(-vals(3))
            sigma(3) = vecs(2,1)*vecs(3,1)*exp(-vals(1))
     &          +vecs(2,2)*vecs(3,2)*exp(-vals(2))
     &          +vecs(2,3)*vecs(3,3)*exp(-vals(3))
            sigma(4) = vecs(1,1)*vecs(3,1)*exp(-vals(1))
     &          +vecs(1,2)*vecs(3,2)*exp(-vals(2))
     &          +vecs(1,3)*vecs(3,3)*exp(-vals(3))
            sigma(5) = vecs(1,1)*vecs(2,1)*exp(-vals(1))
     &          +vecs(1,2)*vecs(2,2)*exp(-vals(2))
     &          +vecs(1,3)*vecs(2,3)*exp(-vals(3))

            r_data(i0,i1,i2,0) = 
     &        c_data(i0,i1,i2,0)
     &        - convec_vals(1,1) - l_inv*(
     &        sigma(0)*rhs_data(i0,i1,i2,0) +
     &        sigma(4)*rhs_data(i0,i1,i2,4) +
     &        sigma(5)*rhs_data(i0,i1,i2,5))

            r_data(i0,i1,i2,1) =
     &        c_data(i0,i1,i2,1)
     &        - convec_vals(2,2) - l_inv*(
     &        sigma(1)*rhs_data(i0,i1,i2,1) +
     &        sigma(3)*rhs_data(i0,i1,i2,3) +
     &        sigma(5)*rhs_data(i0,i1,i2,5))

            r_data(i0,i1,i2,2) =
     &        c_data(i0,i1,i2,2)
     &        - convec_vals(3,3) - l_inv*(
     &        sigma(2)*rhs_data(i0,i1,i2,2) +
     &        sigma(3)*rhs_data(i0,i1,i2,3) +
     &        sigma(4)*rhs_data(i0,i1,i2,4))

            r_data(i0,i1,i2,3) =
     &        c_data(i0,i1,i2,3)
     &        - convec_vals(2,3) - l_inv*(
     &        sigma(1)*rhs_data(i0,i1,i2,3) +
     &        sigma(3)*rhs_data(i0,i1,i2,2) +
     &        sigma(5)*rhs_data(i0,i1,i2,4))

            r_data(i0,i1,i2,4) =
     &        c_data(i0,i1,i2,4)
     &        - convec_vals(1,3) - l_inv*(
     &        sigma(0)*rhs_data(i0,i1,i2,4) +
     &        sigma(4)*rhs_data(i0,i1,i2,2) +
     &        sigma(5)*rhs_data(i0,i1,i2,3))

            r_data(i0,i1,i2,5) =
     &        c_data(i0,i1,i2,5)
     &        - convec_vals(1,2) - l_inv*(
     &        sigma(0)*rhs_data(i0,i1,i2,5) +
     &        sigma(4)*rhs_data(i0,i1,i2,3) +
     &        sigma(5)*rhs_data(i0,i1,i2,1))
          enddo
        enddo
      enddo
      end subroutine

      subroutine log_sum(L,vecs,vals,d,to_ret)

      implicit none
      Integer d
      real*8 L(d,d)
      real*8 vecs(d,d)
      real*8 vals(d)
      real*8 to_ret(d,d)
      real*8 Lij, Lji

      real*8 exp_vals(d)
      Integer i, j, k, kk, ii, jj

      do i=1,d
        do j=1,d
          to_ret(i,j) = 0.d0
        enddo
        exp_vals(i) = exp(vals(i))
      enddo

      do k=1,d
        do kk=1,d
          do i=1,d
            Lij = 0
            do ii=1,d
              do jj=1,d
                Lij = Lij + vecs(jj,i)*L(jj,ii)*vecs(ii,i)
              enddo
            enddo
            to_ret(k,kk) = to_ret(k,kk) +
     &        2*Lij*vecs(k,i)*vecs(kk,i)
          enddo
        enddo
      enddo
      do k=1,d
        do kk=1,d
          do i=1,d; do j=1,d
            if (i /= j) then
              Lij = 0.d0; Lji = 0.d0
              do ii=1,d
                do jj=1,d
                  Lij = Lij + vecs(jj,i)*L(jj,ii)*vecs(ii,j)
                  Lji = Lji + vecs(jj,j)*L(jj,ii)*vecs(ii,i)
                enddo
              enddo
              if(abs(vals(k)-vals(kk))<1d-10) then
                  to_ret(k,kk) = to_ret(k,kk) +
     &              (Lji+Lij)*vecs(k,i)*vecs(kk,j)
              else
                  to_ret(k,kk) = to_ret(k,kk) +
     &              (vals(i)-vals(j))/(exp_vals(i)-exp_vals(j))
     &              *(Lij*exp_vals(j)+Lji*exp_vals(i))
     &              *vecs(k,i)*vecs(kk,j)
              endif
            endif
          enddo; enddo
        enddo
      enddo
      endsubroutine
