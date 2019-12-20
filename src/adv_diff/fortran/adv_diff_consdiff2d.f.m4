c ---------------------------------------------------------------------
c
c Copyright (c) 2006 - 2017 by the IBAMR developers
c All rights reserved.
c
c This file is part of IBAMR.
c
c IBAMR is free software and is distributed under the 3-clause BSD
c license. The full text of the license can be found in the file
c COPYRIGHT at the top level directory of IBAMR.
c
c ---------------------------------------------------------------------

define(NDIM,2)dnl
define(REAL,`double precision')dnl
define(INTEGER,`integer')dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim2d.i)dnl
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Update a quantity using flux differencing with a capacity function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine adv_diff_consdiff_cap2d(
     &     dx,dt,
     &     ifirst0,ilast0,ifirst1,ilast1,
     &     nfluxgc0,nfluxgc1,
     &     nqvalgc0,nqvalgc1,
     &     flux0,flux1,
     &     kappa,kappagc,
     &     qval)
c
      implicit none
c
c     Input.
c
      INTEGER ifirst0,ilast0,ifirst1,ilast1

      INTEGER nfluxgc0,nfluxgc1
      INTEGER nqvalgc0,nqvalgc1

      REAL dx(0:NDIM-1),dt

      REAL flux0(FACE2d0VECG(ifirst,ilast,nfluxgc))
      REAL flux1(FACE2d1VECG(ifirst,ilast,nfluxgc))

      INTEGER kappagc
      REAL kappa(CELL2d(ifirst,ilast,kappagc))
c
c     Input/Output.
c
      REAL qval(CELL2dVECG(ifirst,ilast,nqvalgc))
c
c     Local variables.
c
      INTEGER ic0,ic1,d
      REAL dtdx(0:NDIM-1)
      REAL k
c
c     Update a quantity using flux differencing.
c
      do d = 0,NDIM-1
         dtdx(d) = dt*dx(d)
      enddo

      do ic1 = ifirst1,ilast1
         do ic0 = ifirst0,ilast0
c            do i = 0,1
c              call find_l(ls(ic0+i,ic1)+1.d-12,
c     &                    ls(ic0+i,ic1+1)+1.d-12,
c     &                    dx(1), Ly(i))
c              call find_l(ls(ic0,ic1+i)+1.d-12,
c     &                    ls(ic0+1,ic1+i)+1.d-12,
c     &                    dx(0), Lx(i))
c            enddo
c            r = min(dx(0)/(Lx(0)+1.d-12), dx(0)/(Lx(1)+1.d-12),
c     &              dx(1)/(Ly(0)+1.d-12), dx(1)/(Ly(1)+1.d-12))
c            r = max(r, 1.d0)
c            k = max(kappa(ic0,ic1)*r/(dx(0)*dx(1)), 1.d-2*dx(0)*dx(1))
            k = max(kappa(ic0,ic1)/(dx(0)*dx(1)), 1.0d-2*dx(0)*dx(1))
            qval(ic0,ic1) =
     &           -(flux0(ic0+1,ic1)-flux0(ic0,ic1))
     &            /(dtdx(0)*k)
     &           -(flux1(ic1+1,ic0)-flux1(ic1,ic0))
     &            /(dtdx(1)*k)
         enddo
      enddo
c
      return
      endsubroutine
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Update a quantity using flux differencing.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine adv_diff_consdiff2d(
     &     dx,dt,
     &     ifirst0,ilast0,ifirst1,ilast1,
     &     nfluxgc0,nfluxgc1,
     &     nqvalgc0,nqvalgc1,
     &     flux0,flux1,
     &     qval)
c
      implicit none
c
c     Input.
c
      INTEGER ifirst0,ilast0,ifirst1,ilast1

      INTEGER nfluxgc0,nfluxgc1
      INTEGER nqvalgc0,nqvalgc1

      REAL dx(0:NDIM-1),dt

      REAL flux0(FACE2d0VECG(ifirst,ilast,nfluxgc))
      REAL flux1(FACE2d1VECG(ifirst,ilast,nfluxgc))
c
c     Input/Output.
c
      REAL qval(CELL2dVECG(ifirst,ilast,nqvalgc))
c
c     Local variables.
c
      INTEGER ic0,ic1,d
      REAL dtdx(0:NDIM-1)
c
c     Update a quantity using flux differencing.
c
      do d = 0,NDIM-1
         dtdx(d) = dt*dx(d)
      enddo

      do ic1 = ifirst1,ilast1
         do ic0 = ifirst0,ilast0
            qval(ic0,ic1) =
     &           -(flux0(ic0+1,ic1)-flux0(ic0,ic1))
     &            /(dtdx(0))
     &           -(flux1(ic1+1,ic0)-flux1(ic1,ic0))
     &            /(dtdx(1))
         enddo
      enddo
c
      return
      end
      
      subroutine find_l(phi_l, phi_u, dx, L)
      implicit none
      REAL phi_l, phi_u
      REAL dx, L

      L = 0.d0
      if (phi_L .LE. 0.d0 .AND. phi_U .GE. 0.d0) then
          L = phi_L / (phi_L - phi_U)
      else if (phi_L .GE. 0.d0 .AND. phi_U .LE. 0.d0) then
          L = phi_U / (phi_U - phi_L)
      else if (phi_L .LE. 0.d0 .AND. phi_U .LE. 0.d0) then
          L = 1.d0
      endif
      L = L * dx
      end subroutine
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Update a quantity using flux differencing but include the proper
c     source term to account for a non-discretely divergence free
c     advection velocity.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine adv_diff_consdiffwithdivsource2d(
     &     dx,dt,
     &     ifirst0,ilast0,ifirst1,ilast1,
     &     nfluxgc0,nfluxgc1,
     &     nqfluxgc0,nqfluxgc1,
     &     nufluxgc0,nufluxgc1,
     &     nqvalgc0,nqvalgc1,
     &     flux0,flux1,
     &     qflux0,qflux1,
     &     uflux0,uflux1,
     &     qval)
c
      implicit none
      REAL fourth
      parameter (fourth=0.25d0)
c
c     Input.
c
      INTEGER ifirst0,ilast0,ifirst1,ilast1

      INTEGER nfluxgc0,nfluxgc1
      INTEGER nqfluxgc0,nqfluxgc1
      INTEGER nufluxgc0,nufluxgc1
      INTEGER nqvalgc0,nqvalgc1

      REAL dx(0:NDIM-1),dt

      REAL flux0(FACE2d0VECG(ifirst,ilast,nfluxgc))
      REAL flux1(FACE2d1VECG(ifirst,ilast,nfluxgc))

      REAL qflux0(FACE2d0VECG(ifirst,ilast,nqfluxgc))
      REAL qflux1(FACE2d1VECG(ifirst,ilast,nqfluxgc))

      REAL uflux0(FACE2d0VECG(ifirst,ilast,nufluxgc))
      REAL uflux1(FACE2d1VECG(ifirst,ilast,nufluxgc))
c
c     Input/Output.
c
      REAL qval(CELL2dVECG(ifirst,ilast,nqvalgc))
c
c     Local variables.
c
      INTEGER ic0,ic1,d
      REAL dtdx(0:NDIM-1),divsource
c
c     Update a quantity using flux differencing.
c
      do d = 0,NDIM-1
         dtdx(d) = dt*dx(d)
      enddo

      do ic1 = ifirst1,ilast1
         do ic0 = ifirst0,ilast0
            divsource = (fourth/(dt**2.d0))*
     &           ( qflux0(ic0+1,ic1) + qflux0(ic0,ic1)
     &           + qflux1(ic1+1,ic0) + qflux1(ic1,ic0) )*
     &           ( (uflux0(ic0+1,ic1)-uflux0(ic0,ic1))/dx(0)
     &           + (uflux1(ic1+1,ic0)-uflux1(ic1,ic0))/dx(1) )

            qval(ic0,ic1) = divsource
     &           -(flux0(ic0+1,ic1)-flux0(ic0,ic1))/dtdx(0)
     &           -(flux1(ic1+1,ic0)-flux1(ic1,ic0))/dtdx(1)
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
