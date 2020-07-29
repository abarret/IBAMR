c
c  File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/patchdata/fortran/pdat_m4arrdim1d.i $
c  Package:     SAMRAI patchdata
c  Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision: 1917 $
c  Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
c  Description: m4 include file for dimensioning 1d arrays in FORTRAN routines.
c
define(SAMRAICELL1d,`$1$4-$3:$2$4+$3')dnl
define(SAMRAICELL1d0G,`$1$3:$2$3')dnl
define(SAMRAICELL1dVECG,`$1$4-$3$4:$2$4+$3$4')dnl
define(SAMRAIEDGE1d,`$1$4-$3:$2$4+$3')dnl
define(SAMRAIEDGE1d0G,`$1$3:$2$3')dnl
define(SAMRAIEDGE1dVECG,`$1$4-$3$4:$2$4+$3$4')dnl
define(SAMRAIFACE1d,`$1$4-$3:$2$4+1+$3')dnl
define(SAMRAIFACE1d0G,`$1$3:$2$3+1')dnl
define(SAMRAIFACE1dVECG,`$1$4-$3$4:$2$4+1+$3$4')dnl
define(SAMRAINODE1d,`$1$4-$3:$2$4+1+$3')dnl
define(SAMRAINODE1d0G,`$1$3:$2$3+1')dnl
define(SAMRAINODE1dVECG,`$1$4-$3$4:$2$4+1+$3$4')dnl
define(SAMRAIOUTERFACE1d,`1')dnl
define(SAMRAIOUTERSIDE1d,`1')dnl
define(SAMRAIOUTERNODE1d,`1')dnl
define(SAMRAISIDE1d,`$1$4-$3:$2$4+1+$3')dnl
define(SAMRAISIDE1d0G,`$1$3:$2$3+1')dnl
define(SAMRAISIDE1dVECG,`$1$4-$3$4:$2$4+1+$3$4')dnl
define(CELL1d,`ifelse($3,`0',`SAMRAICELL1d0G($1,$2,0)',`SAMRAICELL1d($1,$2,$3,0)')')dnl 
define(EDGE1d,`ifelse($3,`0',`SAMRAIEDGE1d0G($1,$2,0)',`SAMRAIEDGE1d($1,$2,$3,0)')')dnl 
define(FACE1d,`ifelse($3,`0',`SAMRAIFACE1d0G($1,$2,0)',`SAMRAIFACE1d($1,$2,$3,0)')')dnl 
define(NODE1d,`ifelse($3,`0',`SAMRAINODE1d0G($1,$2,0)',`SAMRAINODE1d($1,$2,$3,0)')')dnl 
define(OUTERFACE1d,`SAMRAIOUTERFACE1d')dnl 
define(OUTERSIDE1d,`SAMRAIOUTERSIDE1d')dnl 
define(OUTERNODE1d,`SAMRAIOUTERNODE1d')dnl 
define(SIDE1d,`ifelse($3,`0',`SAMRAISIDE1d0G($1,$2,0)',`SAMRAISIDE1d($1,$2,$3,0)')')dnl
define(CELL1dVECG,`SAMRAICELL1dVECG($1,$2,$3,0)')dnl 
define(EDGE1dVECG,`SAMRAIEDGE1dVECG($1,$2,$3,0)')dnl 
define(FACE1dVECG,`SAMRAIFACE1dVECG($1,$2,$3,0)')dnl 
define(NODE1dVECG,`SAMRAINODE1dVECG($1,$2,$3,0)')dnl
define(SIDE1dVECG,`SAMRAIFACE1dVECG($1,$2,$3,0)')dnl 
