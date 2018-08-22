c234567
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       Computes d = alpha div s
c
c       where d is vector valued side centered
c       and s is tensor valued cell centered
c       using centered differencing
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine div_tensor_s_to_c_2d(dx, d_data_0, d_data_1, d_gcw,
     &        s_data, s_gcw,  ilower0,
     &        iupper0, ilower1,  iupper1, alpha)
      implicit none
c     INPUTS
      INTEGER ilower0,  iupper0
      INTEGER iupper1,  ilower1
      INTEGER s_gcw,  d_gcw

      real*8 alpha
c     RETURNS
      real*8 d_data_0((ilower0-d_gcw):(iupper0+d_gcw+1),
     &                (ilower1-d_gcw):(iupper1+d_gcw))
      real*8 d_data_1((ilower0-d_gcw):(iupper0+d_gcw),
     &                (ilower1-d_gcw):(iupper1+d_gcw+1))
c     TAU DATA
      real*8 s_data((ilower0-s_gcw):(iupper0+s_gcw),
     &                 (ilower1-s_gcw):(iupper1+s_gcw),
     &                  0:2)
      real*8 dx(0:2)

      INTEGER i0, i1
      real*8 scale0_x, scale0_y
      real*8 scale1_x, scale1_y

      scale0_x = alpha/dx(0)
      scale0_y = alpha/(dx(1)*4.d0)
      scale1_y = alpha/dx(1)
      scale1_x = alpha/(dx(0)*4.d0)

      do i1 = ilower1, iupper1
        do i0 = ilower0, iupper0+1
          d_data_0(i0,i1) =
     &      scale0_x*(s_data(i0, i1,0) - s_data(i0-1, i1,0)) +
     &      scale0_y*(s_data(i0-1, i1+1,2)+s_data(i0, i1+1,2)
     &                -s_data(i0-1, i1-1,2)-s_data(i0, i1-1,2))
        enddo
      enddo
      do i1 = ilower1, (iupper1+1)
        do i0 = ilower0, iupper0
          d_data_1(i0,i1) =
     &      scale1_y*(s_data(i0, i1,1) - s_data(i0, i1-1,1)) +
     &      scale1_x*(s_data(i0+1, i1,2) + s_data(i0+1, i1-1,2)
     &                -s_data(i0-1, i1-1,2) - s_data(i0-1, i1,2))
        enddo
      enddo
      end subroutine
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       Computes d = alpha div s
c
c       where d is vector valued cell centered
c       and s is tensor valued cell centered
c       using centered differencing
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine div_tensor_c_to_c_2d(dx, d_data, d_gcw,
     &        s_data, s_gcw,  ilower0,
     &        iupper0, ilower1,  iupper1, alpha)
      implicit none
c     INPUTS
      INTEGER ilower0,  iupper0
      INTEGER iupper1,  ilower1
      INTEGER s_gcw,  d_gcw
      real*8 alpha
c     RETURNS
      real*8 d_data((ilower0-d_gcw):(iupper0+d_gcw),
     &                (ilower1-d_gcw):(iupper1+d_gcw),
     &                 0:1)
c     TAU DATA
      real*8 s_data((ilower0-s_gcw):(iupper0+s_gcw),
     &                 (ilower1-s_gcw):(iupper1+s_gcw),
     &                  0:2)
      real*8 dx(0:2)
      INTEGER i0, i1
      real*8 scale_x, scale_y

      scale_x = alpha/(2.d0*dx(0))
      scale_y = alpha/(2.d0*dx(1))

      do i1 = ilower1, iupper1
        do i0 = ilower0, iupper0
          d_data(i0,i1,0) =
     &     scale_x*(s_data(i0+1, i1,0) - s_data(i0-1, i1,0)) +
     &     scale_y*(s_data(i0, i1+1,2) - s_data(i0, i1-1,2))
          d_data(i0,i1,1) =
     &     scale_y*(s_data(i0, i1+1,1) - s_data(i0, i1-1,1)) +
     &     scale_x*(s_data(i0+1, i1,2) - s_data(i0-1, i1,2))
        enddo
      enddo
      end subroutine
