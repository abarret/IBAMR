c234567
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       Computes cell centered Oldroyd-B Convective Operator
c
c       where u is vector valued face centered velocity
c       and tau is symmetric log tensor valued cell centered
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine log_tens_conv_u_f_oper_2d
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

      Integer i0, i1
      real*8 u, v
      real*8 q00, q01, q10, q11
      real*8 dq_00_dx, dq_11_dx, dq_01_dx, dq_10_dx
      real*8 dq_00_dy, dq_11_dy, dq_01_dy, dq_10_dy
      real*8 du_dx, dv_dx
      real*8 du_dy, dv_dy
      real*8 scale0_q, scale1_q
      real*8 scale_ux, scale_uy
      real*8 scale_vx, scale_vy
      real*8 l_inv

      real*8 L(2,2), vecs(2,2)
      real*8 vals(2)
      real*8 convec_vals(2,2)
      real*8 sigma(0:2)
      real*8 temp(5)
      Integer info

      l_inv = 1.d0/lambda
      scale_vx = 1.d0/(4.d0*dx(0))
      scale_vy = 1.d0/dx(1)
      scale_ux = 1.d0/dx(0)
      scale_uy = 1.d0/(4.d0*dx(1))
      scale0_q = 1.d0/(2.d0*dx(0))
      scale1_q = 1.d0/(2.d0*dx(1))
      do i1 = ilower1, iupper1
        do i0 = ilower0, iupper0
          q00 = s_data(i0,i1,0)
          q01 = s_data(i0,i1,2)
          q11 = s_data(i0,i1,1)
          du_dx = scale_ux*(u_data_0(i0+1,i1)-u_data_0(i0,i1))
          du_dy = scale_uy*(u_data_0(i0+1,i1+1)+u_data_0(i0,i1+1)
     &            -u_data_0(i0+1,i1-1)-u_data_0(i0,i1-1))
          dv_dy = scale_vy*(u_data_1(i1+1,i0)-u_data_1(i1,i0))
          dv_dx = scale_vx*(u_data_1(i1+1,i0+1)+u_data_1(i1,i0+1)
     &            -u_data_1(i1,i0-1) - u_data_1(i1+1,i0-1))

        L(1,1) = du_dx; L(1,2) = du_dy
        L(2,1) = dv_dx; L(2,2) = dv_dy

        vecs(1,1) = q00; vecs(1,2) = q01
        vecs(2,1) = q01; vecs(2,2) = q11

        call dsyev('V','U',2,vecs,2,vals,temp,5,info)
        if (info /= 0) then
          print *, "\nERROR IN DSYEV!!\n\n"
        endif

        call log_sum(L,vecs,vals,2,convec_vals)

        sigma(0) = vecs(1,1)**2*exp(-vals(1))
     &    + vecs(1,2)**2*exp(-vals(2))
        sigma(1) = vecs(2,1)**2*exp(-vals(1))
     &    + vecs(2,2)**2*exp(-vals(2))
        sigma(2) = vecs(1,1)*vecs(2,1)*exp(-vals(1))
     &    + vecs(1,2)*vecs(2,2)*exp(-vals(2))

        r_data(i0,i1,0) =
     &    c_data(i0,i1,0)
     &    - convec_vals(1,1) - l_inv*(
     &    sigma(0)*rhs_data(i0,i1,0) +
     &    sigma(2)*rhs_data(i0,i1,2))
        r_data(i0,i1,1) =
     &    c_data(i0,i1,1)
     &    - convec_vals(2,2) - l_inv*(
     &    sigma(1)*rhs_data(i0,i1,1) +
     &    sigma(2)*rhs_data(i0,i1,2))
        r_data(i0,i1,2) =
     &    c_data(i0,i1,2)
     &    - convec_vals(1,2) - l_inv*(
     &    sigma(0)*rhs_data(i0,i1,2) +
     &    sigma(2)*rhs_data(i0,i1,1))
        enddo
      enddo
      end subroutine

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       Computes cell centered Oldroyd-B Convective Operator
c
c       where u is vector valued cell centered velocity
c       and tau is symmetric log tensor valued cell centered
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine log_tens_conv_u_c_oper_2d
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
      real*8 u, v
      real*8 q00, q01, q10, q11
      real*8 dq_00_dx, dq_11_dx, dq_01_dx, dq_10_dx
      real*8 dq_00_dy, dq_11_dy, dq_01_dy, dq_10_dy
      real*8 du(ilower0:iupper0,
     &          ilower1:iupper1,
     &          0:1)
      real*8 dv(ilower0:iupper0,
     &          ilower1:iupper1,
     &          0:1)
      real*8 scale0, scale1
      real*8 l_inv

      real*8 L(2,2), vecs(2,2)
      real*8 vals(2)
      real*8 convec_vals(2,2)
      real*8 sigma(0:2)
      real*8 temp(5)
      Integer info

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
          q00 = s_data(i0,i1,0)
          q01 = s_data(i0,i1,2)
          q11 = s_data(i0,i1,1)

          L(1,1) = du(i0,i1,0); L(1,2) = du(i0,i1,1)
          L(2,1) = dv(i0,i1,0); L(2,2) = dv(i0,i1,1)

          vecs(1,1) = q00; vecs(1,2) = q01
          vecs(2,1) = q01; vecs(2,2) = q11

          call dsyev('V','U',2,vecs,2,vals,temp,5,info)
          if (info /= 0) then
            print *, "\nERROR IN DSYEV!!\n\n"
          endif

          call log_sum(L,vecs,vals,2,convec_vals)

          sigma(0) = vecs(1,1)**2*exp(-vals(1))
     &      + vecs(1,2)**2*exp(-vals(2))
          sigma(1) = vecs(2,1)**2*exp(-vals(1))
     &      + vecs(2,2)**2*exp(-vals(2))
          sigma(2) = vecs(1,1)*vecs(2,1)*exp(-vals(1))
     &      + vecs(1,2)*vecs(2,2)*exp(-vals(2))

          r_data(i0,i1,0) =
     &      c_data(i0,i1,0)
     &      - convec_vals(1,1) - l_inv*(
     &      sigma(0)*rhs_data(i0,i1,0) +
     &      sigma(2)*rhs_data(i0,i1,2))
          r_data(i0,i1,1) =
     &      c_data(i0,i1,1)
     &      - convec_vals(2,2) - l_inv*(
     &      sigma(1)*rhs_data(i0,i1,1) +
     &      sigma(2)*rhs_data(i0,i1,2))
          r_data(i0,i1,2) =
     &      c_data(i0,i1,2)
     &      - convec_vals(1,2) - l_inv*(
     &      sigma(0)*rhs_data(i0,i1,2) +
     &      sigma(2)*rhs_data(i0,i1,1))
        enddo
      enddo
      end subroutine

      subroutine log_sum(L,vecs,vals,d,to_ret)

      implicit none
      Integer d
      real*8 L(d,d)
      real*8 Lij, Lji
      real*8 vecs(d,d)
      real*8 vals(d)
      real*8 to_ret(d,d)
      real*8 wij

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
     &        2.d0*Lij*vecs(k,i)*vecs(kk,i)
          enddo
        enddo
      enddo

      do k=1,d
        do kk=1,d
          do i=1,d; do j=1,d
            if (i /= j) then
              Lij = 0
              Lji = 0
              do ii=1,d
                do jj=1,d
                  Lij = Lij + vecs(jj,i)*L(jj,ii)*vecs(ii,j)
                  Lji = Lji + vecs(jj,j)*L(jj,ii)*vecs(ii,i)
                enddo
              enddo
              if(abs(exp_vals(i)-exp_vals(j))<1d-12) then
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
