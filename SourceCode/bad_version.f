!**********************************************************************
!  |*|*\    /*|*|  |*|*|*|*|*|  |*|*\   |*|  |*|*|*|*|*| |*|        |*|
!  |*|\*\  /*/|*|  |*|          |*|\*\  |*|  |*|         |*|  /**\  |*|
!  |*| \*\/*/ |*|  |*|          |*| \*\ |*|  |*|*|*|*|*| |*| /*/\*\ |*|  
!  |*|  \**/  |*|  |*|          |*|  \*\|*|  |*|         |*|/*/  \*\|*|     
!  |*|        |*|  |*|*|*|*|*|  |*|   \*|*|  |*|*|*|*|*| |*|*/    \*|*|   
!_______________________________________________________________________
! MCNEW v2022.1
! NEWtonian aerodynamics calculater by Monte Carlo integration
!
! Copyright (c) 2022 Michiko Ahn Furudate
! Released under the MIT license
! https://opensource.org/licenses/mit-license.php
! 
! Permission is hereby granted, free of charge, to any person obtaining
! a copy of this software and associated documentation files (the
! "Software"), to deal in the Software without restriction, including
! without limitation the rights to use, copy, modify, merge, publish,
! distribute, sublicense, and/or sell copies of the Software, and to
! permit persons to whom the Software is furnished to do so, subject to
! the following conditions:
! 
! The above copyright notice and this permission notice shall be
! included in all copies or substantial portions of the Software.
! 
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
! EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
! MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
! IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
! CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
! TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
! SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
! 
!**********************************************************************
! MCNEW v2022.1
! Newtonian aerodynamics coeffs calculater by Monte Carlo Integration
!
!     Michiko Ahn Furudate, Chungnam National University, Korea
!
!======================================================================
! Modules
! << Module list >>
!   1.module funcs 
!       Functions to calculate xyz-coordinates and normal vectors at 
!       sample points on the given geometries
!   2.module geom
!       Subroutine to read inputs and calculating parameters
!       considering connections of geometries
!_______________________________________________________________________
! module funcs 
!----------------------------------------------------------------------
      module funcs
       implicit none
!
       contains
!      ----------------------------------------------------------------
!      Get xyz coordinate of geometries
!      ----------------------------------------------------------------
       function geom_xyz(gtype,var1,var2,param) result(xyz)
         real(8),intent(in)      :: var1,var2,param(5)
         character(8),intent(in) :: gtype
         real(8) :: xyz(3)
         real(8) :: radius,xlngt,sgn,radcen,flngt,fht
!
         radius    = param(1)
         xlngt     = param(2)
         sgn       = param(3)
         flngt     = param(4)
         fht       = param(5)
         if(gtype=="sphere") then
           xyz(1:3) = sph_xyz(var1,var2,radius)
         elseif(gtype=="cone") then
            xyz(1:3) = cone_xyz(var1,var2,radius,xlngt,sgn)
         elseif(gtype=="flap") then
            xyz(1:3) = flap_xyz(var1,var2,radius,xlngt,sgn,flngt,fht)
         elseif(gtype=="cylinder") then
           xyz(1:3) = cylndr_xyz(var1,var2,radius,xlngt,sgn)
         elseif(gtype=="circle") then
           xyz(1:3) = circle_xyz(var1,var2,radius)
         elseif(gtype=="shoulder") then
           radcen   = param(2)
           xyz(1:3) = shldr_xyz(var1,var2,radius,radcen)
         endif
       end function geom_xyz
!
!      ----------------------------------------------------------------
!      Get xyz coordinate and normal vector of geometries
!      ----------------------------------------------------------------
       function geom_vnor(xyz,gtype,var1,var2,param) result(vnor)
         real(8),intent(out)     :: xyz(3)
         real(8),intent(in)      :: var1,var2,param(5)
         character(8),intent(in) :: gtype
         real(8) :: vnor(4)
         real(8) :: radius,xlngt,sgn,radcen,flngt,fht
!
         radius    = param(1)
         xlngt     = param(2)
         sgn       = param(3)
         flngt     = param(4)
         fht       = param(5)
         if(gtype=="sphere") then
           xyz(1:3)  = sph_xyz(var1,var2,radius)
           vnor(1:4) = sph_normal_xyz(xyz(1:3),radius)
         elseif(gtype=="cone") then
           xyz(1:3)  = cone_xyz(var1,var2,radius,xlngt,sgn)
           vnor(1:4) = cone_normal_xyz(xyz(1:3),radius,xlngt)
        elseif(gtype=="flap") then
           xyz(1:3)  = flap_xyz(var1,var2,radius,xlngt,sgn,flngt,fht)
           vnor(1:4)  = flap_normal_xyz(xyz(1:3),var2,radius,
     &          xlngt,sgn,flngt,fht)
         elseif(gtype=="cylinder") then
           xyz(1:3)  = cylndr_xyz(var1,var2,radius,xlngt,sgn)
           vnor(1:4) = cylndr_normal_xyz(xyz(1:3),radius,xlngt,sgn)
         elseif(gtype=="circle") then
           xyz(1:3)  = circle_xyz(var1,var2,radius)
           vnor(1:4) = circle_normal_xyz(xyz(1:3),radius,sgn)
         elseif(gtype=="shoulder") then
           radcen    = param(2)
           xyz(1:3)  = shldr_xyz(var1,var2,radius,radcen)
           vnor(1:4) = shldr_normal_xyz(xyz(1:3),radcen)
         endif
       end function geom_vnor
!
!      ----------------------------------------------------------------
!      Cylinder
!      ----------------------------------------------------------------
       function cylndr_xyz(u,theta,r,h,sgn) result(xyz)
        real(8),intent(in) :: u,theta,r,h,sgn
        real(8) :: xyz(3)
!
        if(sgn .le. 0.d0)then
          xyz(1) =  h*u
          xyz(2) =  r*dcos(theta)
          xyz(3) =  r*dsin(theta)
        else
          xyz(1) =  r*dcos(theta)
          xyz(2) =  r*dsin(theta)
          xyz(3) =  h*u
        endif
!
        return
       end function cylndr_xyz
!
       function cylndr_normal_xyz(xyz,r,h,sgn) result(vnor)
        real(8),intent(in) :: xyz(3),r,h,sgn
        real(8) :: uu,vnor(4)
        if(sgn .le. 0.d0)then
        vnor(1) = 0.d0
        vnor(2) = h*xyz(2)
        vnor(3) = h*xyz(3)
        vnor(4) = dsqrt(vnor(1)**2+vnor(2)**2+vnor(3)**2)
        else
        vnor(1) = h*xyz(1)
        vnor(2) = h*xyz(2)
        vnor(3) = 0.d0
        vnor(4) = dsqrt(vnor(1)**2+vnor(2)**2+vnor(3)**2)
        endif
        return
       end function cylndr_normal_xyz
!
!      ----------------------------------------------------------------
!      Circle
!      ----------------------------------------------------------------
       function circle_xyz(u,theta,r) result(xyz)
        real(8),intent(in) :: u,theta,r
        real(8) :: xyz(3)
!
        xyz(1) =  0.d0
        xyz(2) =  u*r*dcos(theta)
        xyz(3) =  u*r*dsin(theta)
        return
       end function circle_xyz
!
       function circle_normal_xyz(xyz,r,sgn) result(vnor)
        real(8),intent(in) :: xyz(3),r,sgn
        real(8) :: uu,vnor(4)
        vnor(1) = sgn*dsqrt(xyz(2)**2+xyz(3)**2)
        vnor(2) = 0.d0
        vnor(3) = 0.d0
        vnor(4) = dsqrt(vnor(1)**2+vnor(2)**2+vnor(3)**2)
        return
       end function circle_normal_xyz
!
!      ----------------------------------------------------------------
!      Cone
!      ----------------------------------------------------------------
       function cone_xyz(u,theta,r,h,sgn) result(xyz)
        real(8),intent(in) :: u,theta,r,h,sgn
        real(8) :: xyz(3)
!
        xyz(1) =  u*h*sgn
        xyz(2) =  u*r*dcos(theta)
        xyz(3) =  u*r*dsin(theta)
        return
       end function cone_xyz
!
       function cone_normal_xyz(xyz,r,h) result(vnor)
        real(8),intent(in) :: xyz(3),r,h
        real(8) :: uu,vnor(4)
        uu = xyz(1)/h
        vnor(1) =-r*r*uu
        vnor(2) = h*xyz(2)
        vnor(3) = h*xyz(3)
        vnor(4) = dsqrt(vnor(1)**2+vnor(2)**2+vnor(3)**2)
        return
       end function cone_normal_xyz
!
!      ----------------------------------------------------------------
!      Sphere
!      ----------------------------------------------------------------
       function sph_xyz(theta,phi,r) result(xyz)
        real(8),intent(in) :: theta,phi,r
        real(8) :: xyz(3)
!
        xyz(1) = r*dsin(phi)
        xyz(2) = r*dcos(theta)*dcos(phi)
        xyz(3) = r*dsin(theta)*dcos(phi)
        return
       end function sph_xyz
!
       function sph_normal_xyz(xyz,r) result(vnor)
!       use params, only:r => R_nose
        real(8),intent(in) :: xyz(3),r
        real(8) :: rcosphi,vnor(4)
!
        rcosphi = dsqrt(xyz(2)**2+xyz(3)**2)
        vnor(1) = rcosphi*xyz(1)
        vnor(2) = rcosphi*xyz(2)
        vnor(3) = rcosphi*xyz(3)
        vnor(4) = dsqrt(vnor(1)**2+vnor(2)**2+vnor(3)**2)
        return
       end function sph_normal_xyz
!
!      ----------------------------------------------------------------
!      Shoulder
!      ----------------------------------------------------------------
       function shldr_xyz(theta,phi,rs,rc) result(xyz)
        real(8),intent(in) :: theta,phi,rs,rc
        real(8) :: xyz(3)
        xyz(1) = rs*dsin(phi)
        xyz(2) =(rs*dcos(phi)+rc)*dcos(theta)
        xyz(3) =(rs*dcos(phi)+rc)*dsin(theta)
        return
       end function shldr_xyz
!
       function shldr_normal_xyz(xyz,rc) result(vnor)
        real(8),intent(in) :: xyz(3),rc
        real(8) :: rcosphi,vnor(4)
!
        rcosphi = dsqrt(xyz(2)**2+xyz(3)**2)-rc
        vnor(1) =(rc+rcosphi)*xyz(1)
        vnor(2) = rcosphi*xyz(2)
        vnor(3) = rcosphi*xyz(3)
        vnor(4) = dsqrt(vnor(1)**2+vnor(2)**2+vnor(3)**2)
        return
       end function shldr_normal_xyz
!
!      ----------------------------------------------------------------
!      Flap: currently only works with fore setting
!      ----------------------------------------------------------------
      function flap_params(u,theta,r,h,sgn,l,a,b) result(fparams)
      real(8),intent(in) :: u,theta,r,h,sgn,l,a,b
      real(8) :: fparams(6)

!     c, k coeffs so flaps have appropriate length and height
        
        fparams(1) = r / (a*h) + 1 / l ! k
        fparams(2) = ( 1 + r/a ) / fparams(1) - h ! c
        
!     redefine a so it scales with x as determined by c, k
        
        fparams(3) = a * fparams(1) * (u*h*sgn + fparams(2)) ! au
        
!     flaps along y axis (in x-y plane) and along z axis (along x-z plane)
        
        fparams(4) = fparams(3)*b / sqrt((b*dcos(theta))**2 +
     &       (fparams(3)*dsin(theta))**2) ! r_flap_y
        
        fparams(5) = fparams(3)*b / sqrt((fparams(3)*dcos(theta))**2 +
     &       (b*dsin(theta))**2) ! r_flap_z

!     angle subtended by one flap, for calculation of surface area
        
        fparams(6) = 2*b*sqrt(1 / (u*r)**2 - 1 / fparams(3)**2) ! arc

        return
      end function flap_params

      function flap_xyz(u,theta,r,h,sgn,l,a) result(xyz)
        real(8),intent(in) :: u,theta,r,h,sgn,l,a
        real(8) :: xyz(3)
        real(8) :: fparams(6)
        real(8) :: c,k,au,r_flap_y,r_flap_z
        real(8),parameter :: b=0.025 ! flap thickness
!
        xyz(1) =  u*h*sgn
        xyz(2) =  u*r*dcos(theta)
        xyz(3) =  u*r*dsin(theta)

        fparams = flap_params(u,theta,r,h,sgn,l,a,b)
        c = fparams(2)
        k = fparams(1)
        au = fparams(3)
        r_flap_y = fparams(4)
        r_flap_z = fparams(5)

!     check whether on flap or on surface of cone
        
        if ( h .gt. xyz(1) .and. xyz(1) .gt. h-l ) then

           if ( r_flap_z .gt. xyz(1) * r/h ) then
              xyz(2) =  r_flap_z*dcos(theta)
              xyz(3) =  r_flap_z*dsin(theta)
           else if ( r_flap_y .gt. xyz(1) * r/h ) then
              xyz(2) =  r_flap_y*dcos(theta)
              xyz(3) =  r_flap_y*dsin(theta)
           end if
        end if
!        
        return
       end function flap_xyz
!
       function flap_normal_xyz(xyz,theta,r,h,sgn,l,a) result(vnor)
        real(8),intent(in) :: xyz(3),theta,r,h,sgn,l,a
        real(8) :: uu,vnor(4),au,norm,arc
        real(8) :: fparams(6)
        real(8) :: PI
        real(8),parameter :: b=0.025
        PI = 3.141592653
        
!     check whether on flap
        if ( sqrt(xyz(2)**2 + xyz(3)**2) .gt. xyz(1)/h * r ) then
           
           vnor(1) = 0.0        ! approximate normal vector as theta vec
           fparams = flap_params(xyz(1)/h,theta,r,h,sgn,l,a,b)
           au = fparams(3)
           arc = fparams(6)
           norm = abs(h*(au - r)*2/arc)
           
           if ((PI/2 .ge. theta .and. theta .gt. PI/4) .or.
     &          (PI .ge. theta .and. theta .gt. 3*PI/4) .or.
     &          (3*PI/2 .ge. theta .and. theta .gt. 5*PI/4) .or.
     &          (2*PI .ge. theta .and. theta .gt. 7*PI/4)) then
              vnor(2) = dsin(theta) * norm
              vnor(3) = -dcos(theta) * norm
              vnor(4) = norm
           else
              vnor(2) = -dsin(theta) * norm
              vnor(3) = dcos(theta) * norm
              vnor(4) = norm
           end if
        else                    ! if not on flap
           uu = xyz(1)/h
           vnor(1) =-r*r*uu
           vnor(2) = h*xyz(2)
           vnor(3) = h*xyz(3)
           vnor(4) = dsqrt(vnor(1)**2+vnor(2)**2+vnor(3)**2)
        end if

        return
       end function flap_normal_xyz
!
!      ----------------------------------------------------------------
!      Rotation transformation
!      ----------------------------------------------------------------
       function rot_body(unor,alpha,beta) result(unor_rot)
        implicit none
        real(8),intent(in) :: unor(3)
        real(8),intent(in) :: alpha, beta
        real(8) :: rot_pitch(3,3)
     &            ,rot_yaw(3,3)
        real(8),dimension(3) :: temp
        real(8),dimension(3) :: unor_rot
!
        rot_pitch(1,1) = dcos(alpha)
        rot_pitch(1,2) = 0.d0
        rot_pitch(1,3) = dsin(alpha)
        rot_pitch(2,1) = 0.d0
        rot_pitch(2,2) = 1.d0
        rot_pitch(2,3) = 0.d0
        rot_pitch(3,1) =-dsin(alpha)
        rot_pitch(3,2) = 0.d0
        rot_pitch(3,3) = dcos(alpha)
!
        rot_yaw(1,1) = dcos(beta)
        rot_yaw(1,2) =-dsin(beta)
        rot_yaw(1,3) = 0.d0
        rot_yaw(2,1) =-dsin(beta)
        rot_yaw(2,2) = dcos(beta)
        rot_yaw(2,3) = 0.d0
        rot_yaw(3,1) = 0.d0
        rot_yaw(3,2) = 0.d0
        rot_yaw(3,3) = 1.d0
!
        temp=matmul(rot_pitch,unor)
        unor_rot=matmul(rot_yaw,temp)
!
!       temp=matmul(rot_yaw,unor)
!       unor_rot=matmul(rot_pitch,temp)
!
        return
       end function rot_body
!
      end module funcs
!
!======================================================================
! Module geoms
!----------------------------------------------------------------------
      module geoms
       use funcs, only: sph_xyz
     &                 ,cone_xyz
     &                 ,circle_xyz
     &                 ,cylndr_xyz
     &                 ,shldr_xyz
     &                 ,flap_xyz
       implicit none
       real(8),parameter :: pi=dacos(-1.d0)
       integer :: nblk
       character(len=10),allocatable :: gtype(:),ptype(:)
       integer,allocatable :: nshld
       integer,allocatable :: icnctm(:) ! block in the left
     &                       ,icnctp(:) ! block in the right
     &                       ,ishld(:) ! block no. of of shoulder
       real(8),allocatable :: radius(:) ! radius
     &                       ,radius1(:) ! radius
     &                       ,radius2(:) ! radius
     &                       ,radcen(:) ! radius
     &                       ,ahalf(:)  ! half angle
     &                       ,xlngt(:)  ! 
     &                       ,xcntr(:)  ! 
     &                       ,rfac(:)   !
     &                       ,flngt(:)  ! flap length
     &                       ,fht(:)    ! flap height
       real(8),allocatable :: v1min(:)  ! 
     &                       ,v1max(:)  ! 
     &                       ,v2min(:)  ! 
     &                       ,v2max(:)  ! 
       real(8),allocatable :: xmin(:)  ! 
     &                       ,xmax(:)  ! 
     &                       ,ymin(:)  ! 
     &                       ,ymax(:)  ! 
     &                       ,zmin(:)  ! 
     &                       ,zmax(:)  ! 
       real(8) :: x_nose
       real(8) :: area_ref  ! reference area  [m2]
     &           ,len_ref   ! reference radius [m]
!
       contains
!
!     -----------------------------------------------------------------
!     READ GEOMETRY INPUT FILES and setting parameters
!     -----------------------------------------------------------------
        subroutine set_geom_params
          character(len=10) :: cread
          character(len=1) :: cdum
          integer :: i,ib,iblk,icnp,icnm
          real(8) :: h0,angle,angle2,rad,l,a
          real(8) :: xyz(3),sgn,v2,v1mn,v1mx
!
          write(*,*)"input: number of blocks"
          read(*,*) nblk
          write(*,*) nblk
!
          allocate(gtype(nblk))
          allocate(ptype(nblk))
          allocate(radius(nblk))
          allocate(radius1(nblk))
          allocate(radius2(nblk))
          allocate(radcen(nblk))
          allocate(ahalf(nblk))
          allocate(flngt(nblk))
          allocate(fht(nblk))
          allocate(rfac(nblk))
          allocate(xcntr(nblk))
          allocate(xlngt(nblk))
          allocate(v1min(nblk))
          allocate(v1max(nblk))
          allocate(v2min(nblk))
          allocate(v2max(nblk))
          allocate(icnctm(nblk))
          allocate(icnctp(nblk))
          allocate(ishld(nblk))
          allocate(xmin(nblk))
          allocate(xmax(nblk))
          allocate(ymin(nblk))
          allocate(ymax(nblk))
          allocate(zmin(nblk))
          allocate(zmax(nblk))
!
          icnctm(:) = 0
          icnctp(:) = 0
          ishld(:)  = 0
          radius(:) = 0.d0
          radius1(:) = 0.d0
          radius2(:) = 0.d0
          radcen(:) = 0.d0
          ahalf(:)  = 0.d0
          flngt(:)  = 0.d0
          fht(:)  = 0.d0
          xlngt(:)  = 0.d0
          xcntr(:)  = 0.d0
          rfac(:)   = 0.d0
          v1min(:)  = 0.d0
          v1max(:)  = 0.d0
          v2min(:)  = 0.d0
          v2max(:)  = 0.d0
          xmin(:)  = 0.d0
          xmax(:)  = 0.d0
          ymin(:)  = 0.d0
          ymax(:)  = 0.d0
          zmin(:)  = 0.d0
          zmax(:)  = 0.d0
!
          do iblk = 1,nblk
            read(*,*) cdum 
            write(*,*)"input: block numnber"
            read(*,*) i
            write(*,*) i
            write(*,*)"input: shapes of blocks"
            read(*,*) gtype(i) 
            write(*,*) trim(gtype(i)) 
            write(*,*)"input: position of blocks"
            read(*,*) ptype(i) 
            write(*,*) trim(ptype(i)) 
            write(*,*)"input: Neighboring block numbers"
            read(*,*) icnctm(i),icnctp(i)
            write(*,*) icnctm(i),icnctp(i)
            if(gtype(i)=="sphere") then
               write(*,*)"input: radius of spherer, [m]"
               read(*,*) radius(i) 
               write(*,*) radius(i) 
            elseif(gtype(i)=="cone") then
               write(*,*)"input: base radius of cone, [m]"
               read(*,*) radius(i) 
               write(*,*) radius(i) 
               write(*,*)"input: Half angle of cone, [deg]"
               read(*,*) ahalf(i) 
               write(*,*) ahalf(i)
            elseif(gtype(i)=="flap") then
               write(*,*)"input: base radius of cone, [m]"
               read(*,*) radius(i) 
               write(*,*) radius(i) 
               write(*,*)"input: Half angle of cone, [deg]"
               read(*,*) ahalf(i) 
               write(*,*) ahalf(i)
               write(*,*)"input: Longitudinal length of flap, [m]"
               read(*,*) flngt(i)
               write(*,*) flngt(i)
               write(*,*)"input: Max height of flap, [m]"
               read(*,*) fht(i)
               write(*,*) fht(i)
            elseif(gtype(i)=="shoulder") then
               write(*,*)"input: radius of arc, [m]"
               read(*,*) radius(i) 
               write(*,*) radius(i) 
               if(ptype(i)=="sph-cone".or.
     &            ptype(i)=="sph-sph" .or. 
     &            ptype(i)=="cone-sph" ) then
                 write(*,*)"input: radius of circle of intercetion with"
     &                    ," the nose sphere, [m]"
                 read(*,*) radius1(i) 
                 write(*,*) radius1(i) 
               endif
               if(ptype(i)=="sph-sph" )then
                 write(*,*)"input: radius of circle of intercetion with"
     &                    ," the tail sphere, [m]"
                 read(*,*) radius2(i) 
                 write(*,*) radius2(i) 
               endif
            elseif(gtype(i)=="cylinder") then
               write(*,*)"input: base radius, [m]"
               read(*,*) radius(i) 
               write(*,*) radius(i) 
               write(*,*)"input: Height, [m]"
               read(*,*) xlngt(i) 
               write(*,*) xlngt(i) 
            elseif(gtype(i)=="circle") then
               write(*,*)"input: radius, [m]"
               read(*,*) radius(i) 
               write(*,*) radius(i) 
            else
              write(*,*)"!! Error !! : Unknow shape type"
            endif
          enddo
!
! Print input summary 
        write(*,*) " "
        write(*,*) " ========","============","============"
        write(*,*) " Block # "," Shape      "," Position   "
        write(*,*) " --------"," ---------- "," ---------- "
        do i=1,nblk
          write(*,"((2x,i4,2x),2(x,a10,x))") i,trim(gtype(i))
     &                                        ,trim(ptype(i))
        enddo
        write(*,*) " ========","============","============"
!
!----------------------------------------------------------------------
! Shoulder blocks
!----------------------------------------------------------------------
        nshld = 0
        do iblk=1,nblk
          if(gtype(iblk).eq."shoulder") then
             nshld = nshld+1
             ishld(nshld) = iblk
          endif
        enddo
!
        do i = 1,nshld
          iblk = ishld(i)
          icnm = icnctm(iblk)
          icnp = icnctp(iblk)
          v1min(iblk) =  0.0d0
          v1max(iblk) =  2.0d0 * pi
          xcntr(iblk) =  0.d0
          rad = radius(iblk)
          if( ptype(iblk)=="cone-cone" .or. 
     &       (gtype(icnm)=="cone" .and. gtype(icnp)=="cone")) then
            v2min(iblk) =  -ahalf(icnm)/180.d0*pi
            v2max(iblk) =   ahalf(icnp)/180.d0*pi
            angle  =0.5d0*(ahalf(icnm)+ahalf(icnp))/180.d0*pi
            angle2 =0.5d0*(ahalf(icnm)-ahalf(icnp))/180.d0*pi
            xcntr(iblk) = rad/dcos(angle)*dsin(angle2)
            if(radius(icnp).gt.0.d0 .and. radius(icnm).gt.0.d0)then
              rfac(iblk)  =dcos(angle2/180.d0*pi)/dcos(angle/180.d0*pi)
              radcen(iblk)= radius(icnm)-rad*rfac(iblk)
            else
              if(radius(icnm).lt.0.d0)then
                radcen(iblk)= dabs(radius(icnm))-rad
              elseif(radius(icnp).lt.0.d0)then
                radcen(iblk)= dabs(radius(icnp))-rad
              endif
              radius(icnm)= rad/dcos(angle)*dcos(angle2)+radcen(iblk)
              radius(icnp)= rad/dcos(angle)*dcos(angle2)+radcen(iblk)
            endif 
          elseif( ptype(iblk)=="sph-cone" .or. 
     &       (gtype(icnm)=="sphere" .and. gtype(icnp)=="cone")) then
            if(radius1(iblk).lt.0.d0)then
              rfac(iblk)   = -radius1(iblk)/radius(icnm)
              angle        = dacos(rfac(iblk))
              radcen(iblk) = -radius1(iblk)-radius(iblk)*dcos(angle)
            else
              rfac(iblk)   =   (radius1(iblk)-rad)
     &                       / (radius(icnm) -rad)
              angle        = dacos(rfac(iblk))
              radcen(iblk) = radius1(iblk)-rad
              radius1(iblk)= rad*dcos(angle)+radcen(iblk)
            endif
            v2min(iblk)  = -angle
            v2max(iblk)  =  ahalf(icnp)/180.d0*pi
            radius(icnp) =  rad/dcos(v2max(iblk))+radcen(iblk)
          elseif( ptype(iblk)=="cone-sph" .or. 
     &      (gtype(icnm)=="cone" .and. gtype(icnp)=="sphere")) then
            rfac(iblk)   = radius1(iblk)/radius(icnp)
            angle        = dacos(rfac(iblk))
            radcen(iblk) = radius(icnp)-radius(iblk)*rfac(iblk)
            v2min(iblk) = -ahalf(icnm)/180.d0*pi
            v2max(iblk) =  angle
            radcen(iblk) = radius1(iblk)-radius(iblk)*dcos(angle)
          elseif( ptype(iblk)=="sph-sph" .or. 
     &      (gtype(icnm)=="sphere" .and. gtype(icnp)=="sphere")) then
            rfac(iblk)   = radius1(iblk)/radius(icnm)
            angle        = dacos(rfac(iblk))
            angle2      =  dacos(radius2(iblk)/radius(icnp))
            v2min(iblk) = -angle
            v2max(iblk) =  angle2
            radcen(iblk) = radius1(iblk)-radius(iblk)*dcos(angle)
           endif
           xyz(1:3)   = shldr_xyz(0.d0,v2min(iblk),rad,radcen(iblk))
           xmin(iblk) = xyz(1)
           xyz(1:3)   = shldr_xyz(0.d0,v2max(iblk),rad,radcen(iblk))
           xmax(iblk) = xyz(1)
           xyz(1:3)   = shldr_xyz( 1.0d0*pi,0.d0,rad,radcen(iblk))
           ymin(iblk) = xyz(2)
           xyz(1:3)   = shldr_xyz( 0.d0*pi,0.d0,rad,radcen(iblk))
           ymax(iblk) = xyz(2)
           xyz(1:3)   = shldr_xyz( 1.5d0*pi,0.d0,rad,radcen(iblk))
           zmin(iblk) = xyz(3)
           xyz(1:3)   = shldr_xyz( 0.5d0*pi,0.d0,rad,radcen(iblk))
           zmax(iblk) = xyz(3)
         enddo
!
!--------------------------------------------------------------------i-
! Sphere blocks
!----------------------------------------------------------------------
!
        do iblk = 1,nblk
         icnm = icnctm(iblk)
         icnp = icnctp(iblk)
         if(gtype(iblk)=="sphere") then
           v1min(iblk) =  0.0d0
           v1max(iblk) =  2.0d0 * pi
           if(ptype(iblk)=="nose") then
             if(icnp.ne.0) then
               if(gtype(icnp)=="cone".or.gtype(icnp)=="flap") then
                angle       = ahalf(icnp)/180.d0 * pi
                v2min(iblk) = -0.5d0 * pi
                v2max(iblk) = -angle
                xcntr(iblk) = dabs(radius(icnp))/dtan(-angle)
     &                       -radius(iblk)/dsin(-angle)
               elseif(gtype(icnp)=="shoulder") then
                angle      =  dacos(radius1(icnp)/radius(iblk))
                v2min(iblk)= -0.5d0 * pi
                v2max(iblk)= -angle
                xcntr(iblk) = radius(iblk)*dsin(angle)
     &                       -radius(icnp)*dsin(angle)
               elseif(gtype(icnp)=="cylinder") then
                v2min(iblk) = -0.5d0 * pi
                v2max(iblk) =  0.0d0 * pi
                xcntr(iblk) =  0.d0
               else
                call print_err(iblk,icnp)
               endif
             else
               v2min(iblk) = -0.5d0 * pi
               v2max(iblk) =  0.0d0 * pi
               xcntr(iblk) =  0.d0
             endif
           elseif(ptype(iblk).eq."tail") then
             if(icnm.ne.0) then
               if(gtype(icnm)=="cone".or.gtype(icnp)=="flap") then
               angle       = ahalf(icnm)/180.d0 * pi
               v2min(iblk) = angle
               v2max(iblk) = 0.5d0 * pi
               xcntr(iblk) =-dabs(radius(icnm))/dtan(-angle)
     &                      +radius(iblk)/dsin(-angle)
               elseif(gtype(icnm)=="shoulder") then
                angle      =  dacos(radius1(icnm)/radius(iblk))
                v2min(iblk)=  angle
                v2max(iblk)=  0.5d0 * pi
                xcntr(iblk) =-radius(iblk)*dsin(angle)
     &                       +radius(icnm)*dsin(angle)
               elseif(gtype(icnm)=="cylinder") then
                v2min(iblk) =  0.0d0 * pi
                v2max(iblk) =  0.5d0 * pi
                xcntr(iblk) =  0.d0
               else
                call print_err(iblk,icnm)
               endif
             else
               v2min(iblk) =  0.0d0 * pi
               v2max(iblk) =  0.5d0 * pi
               xcntr(iblk) =  0.d0
             endif
           else ! Full sphere
             v2min(iblk) = -0.5d0 * pi
             v2max(iblk) =  0.5d0 * pi
             xcntr(iblk) =  0.d0
           endif
!
           xyz(1:3)   = sph_xyz(0.d0,v2min(iblk),radius(iblk))
           xmin(iblk) = xyz(1)
           xyz(1:3)   = sph_xyz(0.d0,v2max(iblk),radius(iblk))
           xmax(iblk) = xyz(1)
           xyz(1:3)   = sph_xyz( 1.0d0*pi,v2max(iblk),radius(iblk))
           ymin(iblk) = xyz(2)
           xyz(1:3)   = sph_xyz( 0.d0*pi,v2max(iblk),radius(iblk))
           ymax(iblk) = xyz(2)
           xyz(1:3)   = sph_xyz( 1.5d0*pi,v2max(iblk),radius(iblk))
           zmin(iblk) = xyz(3)
           xyz(1:3)   = sph_xyz( 0.5d0*pi,v2max(iblk),radius(iblk))
           zmax(iblk) = xyz(3)
         endif
       enddo
!
!----------------------------------------------------------------------
! Cone blocks
!----------------------------------------------------------------------
        do iblk = 1,nblk
          if(gtype(iblk)=="cone") then
            icnm  = icnctm(iblk)
            icnp  = icnctp(iblk)
            angle = ahalf(iblk)/180.d0 * pi
            rad   = radius(iblk)
            xlngt(iblk) = rad/dtan(angle)
!
            if(ptype(iblk)=="fore") then
              sgn         = 1.d0
              if(icnp.ne.0) then
                if(gtype(icnp)=="shoulder") then
                   v1max(iblk) = 1.d0-(radius(icnp)*dsin(angle)
     &                               -xcntr(icnp))/xlngt(iblk)
                else
                  v1max(iblk) = 1.d0
                endif
              else
                v1max(iblk) = 1.d0
              endif
!
              xcntr(iblk) = xcntr(iblk)-sgn*xlngt(iblk)
!
              if(icnm.ne.0) then
                if(gtype(icnm)=="sphere") then
                  v1min(iblk) = radius(icnm)/rad*dcos(angle)
                elseif(gtype(icnm)=="circle") then
                  v1min(iblk) = radius(icnm)/dtan(angle)/xlngt(iblk)
                elseif(gtype(icnm)=="cone") then
                  if(ptype(icnm)=="rear") then
                    if(radius(icnm)==radius(iblk)) then
                      v1min(iblk) = 0.d0
                    else
                      print *,"The fore cone block ",iblk, 
     &                        " can be connected to the rear cone block"
     &                        ,icnm," only when they have the same base"
     &                       ," rasius."
                      call exit
                    endif
                  else
                    if(ahalf(icnm).ge.ahalf(iblk)) then
                      v1min(iblk) = radius(icnm)/dtan(angle)/xlngt(iblk)
                      xcntr(1:icnm) = xcntr(1:icnm) 
     &                             +xcntr(iblk)*(1.d0-v1min(iblk))
                    else
                      print *,"The cone block ",iblk,"can be connected"
     &                       ," to the fore cone block",icnm,
     &                        " only when the half angle of block",iblk,
     &                        " is smaller then that of ",icnm
                      call exit
                    endif
                  endif
                else
                  call print_err(iblk,icnm)
                endif
              else
                v1min(iblk) = 0.d0
              endif
              v2min(iblk) = 0.0d0
              v2max(iblk) = 2.0d0 * pi
!
            elseif(ptype(iblk)=="rear") then
              sgn         =-1.d0
              if(icnm.ne.0) then
                if(gtype(icnm)=="shoulder") then
                  v1max(iblk) = 1.d0-(radius(icnm)*dsin(angle)
     &                               +xcntr(icnm))/xlngt(iblk)
                else
                  v1max(iblk) = 1.d0
                endif
              endif
!
              xcntr(iblk) = xcntr(iblk)-sgn*xlngt(iblk)
!
              if(icnp.ne.0) then
                if(gtype(icnctp(iblk))=="sphere") then
                  v1min(iblk) = radius(icnp)/rad*dcos(angle)
                elseif(gtype(icnp)=="circle") then
                  v1min(iblk) = radius(icnp)/dtan(angle)/xlngt(iblk)
                elseif(gtype(icnp)=="cone") then
                  if(ptype(icnp)=="fore") then
                    call print_err(iblk,icnp)
                  else
                    if(ahalf(icnp).ge.ahalf(iblk)) then
                      v1min(iblk) = radius(icnp)/dtan(angle)/xlngt(iblk)
                      xcntr(icnp) = xcntr(iblk)*(1.d0-v1min(iblk))
                    else
                      print *,"The cone block ",iblk,"can be connected"
     &                       ," to the fore cone block",icnp,
     &                        " only when the half angle of block",iblk,
     &                        " is larger then that of ",icnp
                      call exit
                    endif
                  endif
                endif
              else
                v1min(iblk) = 0.d0
              endif
              v2min(iblk) =  0.0d0
              v2max(iblk) =  2.0d0 * pi
            endif
!
            v1mn = v1min(iblk)
            v1mx = v1max(iblk)
            rad  = radius(iblk)
            if(sgn.gt.0) then
              xyz(1:3) = cone_xyz(v1mn,0.0d0*pi,rad,xlngt(iblk),sgn)
              xmin(iblk) = xyz(1)
              xyz(1:3) = cone_xyz(v1mx,0.0d0*pi,rad,xlngt(iblk),sgn)
              xmax(iblk) = xyz(1)
              ymax(iblk) = xyz(2)
            else
              xyz(1:3) = cone_xyz(v1mx,0.0d0*pi,rad,xlngt(iblk),sgn)
              xmin(iblk) = xyz(1)
              ymax(iblk) = xyz(2)
              xyz(1:3) = cone_xyz(v1mn,0.0d0*pi,rad,xlngt(iblk),sgn)
              xmax(iblk) = xyz(1)
            endif
            xyz(1:3) = cone_xyz(v1mx,0.5d0*pi,rad,xlngt(iblk),sgn)
            zmax(iblk) = xyz(3)
            xyz(1:3) = cone_xyz(v1mx,1.0d0*pi,rad,xlngt(iblk),sgn)
            ymin(iblk) = xyz(2)
            xyz(1:3) = cone_xyz(v1mx,1.5d0*pi,rad,xlngt(iblk),sgn)
            zmin(iblk) = xyz(3)
          endif
        enddo
!
        
!----------------------------------------------------------------------
! Flap blocks
!----------------------------------------------------------------------
        do iblk = 1,nblk
          if(gtype(iblk)=="flap") then
            icnm  = icnctm(iblk)
            icnp  = icnctp(iblk)
            angle = ahalf(iblk)/180.d0 * pi
            rad   = radius(iblk)
            xlngt(iblk) = rad/dtan(angle)
!
            if(ptype(iblk)=="fore") then
              sgn         = 1.d0
              if(icnp.ne.0) then
                if(gtype(icnp)=="shoulder") then
                   v1max(iblk) = 1.d0-(radius(icnp)*dsin(angle)
     &                               -xcntr(icnp))/xlngt(iblk)
                else
                  v1max(iblk) = 1.d0
                endif
              else
                v1max(iblk) = 1.d0
              endif
!
              xcntr(iblk) = xcntr(iblk)-sgn*xlngt(iblk)
!
              if(icnm.ne.0) then
                if(gtype(icnm)=="sphere") then
                  v1min(iblk) = radius(icnm)/rad*dcos(angle)
                elseif(gtype(icnm)=="circle") then
                  v1min(iblk) = radius(icnm)/dtan(angle)/xlngt(iblk)
                elseif(gtype(icnm)=="cone") then
                  if(ptype(icnm)=="rear") then
                    if(radius(icnm)==radius(iblk)) then
                      v1min(iblk) = 0.d0
                    else
                      print *,"The fore flap block ",iblk, 
     &                        " can be connected to the rear cone block"
     &                        ,icnm," only when they have the same base"
     &                       ," rasius."
                      call exit
                    endif
                  else
                    if(ahalf(icnm).ge.ahalf(iblk)) then
                      v1min(iblk) = radius(icnm)/dtan(angle)/xlngt(iblk)
                      xcntr(1:icnm) = xcntr(1:icnm) 
     &                             +xcntr(iblk)*(1.d0-v1min(iblk))
                    else
                      print *,"The cone block ",iblk,"can be connected"
     &                       ," to the fore cone block",icnm,
     &                        " only when the half angle of block",iblk,
     &                        " is smaller then that of ",icnm
                      call exit
                    endif
                  endif
                else
                  call print_err(iblk,icnm)
                endif
              else
                v1min(iblk) = 0.d0
              endif
              v2min(iblk) = 0.0d0
              v2max(iblk) = 2.0d0 * pi
!
            elseif(ptype(iblk)=="rear") then
              sgn         =-1.d0
              if(icnm.ne.0) then
                if(gtype(icnm)=="shoulder") then
                  v1max(iblk) = 1.d0-(radius(icnm)*dsin(angle)
     &                               +xcntr(icnm))/xlngt(iblk)
                else
                  v1max(iblk) = 1.d0
                endif
              endif
!
              xcntr(iblk) = xcntr(iblk)-sgn*xlngt(iblk)
!
              if(icnp.ne.0) then
                if(gtype(icnctp(iblk))=="sphere") then
                  v1min(iblk) = radius(icnp)/rad*dcos(angle)
                elseif(gtype(icnp)=="circle") then
                  v1min(iblk) = radius(icnp)/dtan(angle)/xlngt(iblk)
                elseif(gtype(icnp)=="cone") then
                  if(ptype(icnp)=="fore") then
                    call print_err(iblk,icnp)
                  else
                    if(ahalf(icnp).ge.ahalf(iblk)) then
                      v1min(iblk) = radius(icnp)/dtan(angle)/xlngt(iblk)
                      xcntr(icnp) = xcntr(iblk)*(1.d0-v1min(iblk))
                    else
                      print *,"The cone block ",iblk,"can be connected"
     &                       ," to the fore cone block",icnp,
     &                        " only when the half angle of block",iblk,
     &                        " is larger then that of ",icnp
                      call exit
                    endif
                  endif
                endif
              else
                v1min(iblk) = 0.d0
              endif
              v2min(iblk) =  0.0d0
              v2max(iblk) =  2.0d0 * pi
            endif
!
            v1mn = v1min(iblk)
            v1mx = v1max(iblk)
            rad  = radius(iblk)
            l    = flngt(iblk)
            a    = fht(iblk)
            if(sgn.gt.0) then
              xyz(1:3)=flap_xyz(v1mn,0.d0,rad,xlngt(iblk),l,a,sgn)
              xmin(iblk) = xyz(1)
              xyz(1:3)=flap_xyz(v1mx,0.d0,rad,xlngt(iblk),l,a,sgn)
              xmax(iblk) = xyz(1)
              ymax(iblk) = xyz(2)
            else
              xyz(1:3) = cone_xyz(v1mx,0.0d0*pi,rad,xlngt(iblk),sgn)
              xmin(iblk) = xyz(1)
              ymax(iblk) = xyz(2)
              xyz(1:3) = cone_xyz(v1mn,0.0d0*pi,rad,xlngt(iblk),sgn)
              xmax(iblk) = xyz(1)
            endif
            xyz(1:3)=flap_xyz(v1mx,.5*pi,rad,xlngt(iblk),l,a,sgn)
            zmax(iblk) = xyz(3)
            xyz(1:3)=flap_xyz(v1mx,pi,rad,xlngt(iblk),l,a,sgn)
            ymin(iblk) = xyz(2)
            xyz(1:3)=flap_xyz(v1mx,1.5*pi,rad,xlngt(iblk),l,a,sgn)
            zmin(iblk) = xyz(3)
          endif
        enddo
!
!----------------------------------------------------------------------
! Cylinder  blocks
!----------------------------------------------------------------------
        do iblk = 1,nblk
          if(gtype(iblk)=="cylinder") then
            if(ptype(iblk)=="horizon") then
              sgn = -1.d0
              v1min(iblk) = -0.5d0
              v1max(iblk) =  0.5d0
            else
              sgn = 1.d0
              v1min(iblk) = 0.d0
              v1max(iblk) = 1.d0
            endif
            v2min(iblk) = 0.0d0
            v2max(iblk) = 2.0d0 * pi
            xcntr(iblk) = 0.d0
!
            rad = radius(iblk)
            v1mn = v1min(iblk)
            v1mx = v1max(iblk)
            xyz(1:3) = cylndr_xyz(v1mn,0.0d0*pi,rad,xlngt(iblk),sgn)
            xmin(iblk) = xyz(1)
            xyz(1:3) = cylndr_xyz(v1mx,0.0d0*pi,rad,xlngt(iblk),sgn)
            xmax(iblk) = xyz(1)
            ymax(iblk) = xyz(2)
            xyz(1:3) = cylndr_xyz(v1mx,0.5d0*pi,rad,xlngt(iblk),sgn)
            zmax(iblk) = xyz(3)
            xyz(1:3) = cylndr_xyz(v1mx,1.0d0*pi,rad,xlngt(iblk),sgn)
            ymin(iblk) = xyz(2)
            xyz(1:3) = cylndr_xyz(v1mx,1.5d0*pi,rad,xlngt(iblk),sgn)
            zmin(iblk) = xyz(3)
          endif
        enddo
!
!----------------------------------------------------------------------
! Circle blocks
!----------------------------------------------------------------------
        do iblk = 1,nblk
          icnm = icnctm(iblk)
          icnp = icnctp(iblk)
          if(gtype(iblk)=="circle") then
            v1min(iblk) = 0.d0
            v1max(iblk) = 1.d0
            v2min(iblk) = 0.0d0
            v2max(iblk) = 2.0d0 * pi
            if (icnm .ne. 0) then
              if(gtype(icnm)=="cylinder")then
                if(ptype(iblk)=="top") then
                  call print_err(iblk,icnm)
                elseif(ptype(iblk)=="bottom") then
                  xcntr(iblk) = xlngt(icnm)
                endif
              elseif(gtype(icnm)=="cone".or.gtype(icnp)=="flap")then
                if(ptype(icnm)=="fore") sgn = 1.d0
                if(ptype(icnm)=="rear") sgn =-1.d0
                angle = ahalf(icnm)/180.d0 * pi
                xcntr(iblk) = xcntr(icnm)+sgn*radius(iblk)/dtan(angle)
              else
                xcntr(iblk) = 0.d0
              endif
            endif
            if (icnp .ne. 0) then
              if(gtype(icnp)=="cylinder")then
                if(ptype(iblk)=="top") then
                  xcntr(iblk) = 0.d0
                elseif(ptype(iblk)=="bottom") then
                   call print_err(iblk,icnp)
                endif
              elseif(gtype(icnp)=="cone".or.gtype(icnp)=="flap")then
                angle = ahalf(icnp)/180.d0 * pi
                xcntr(iblk) = xcntr(icnp)+radius(iblk)/dtan(angle)
              else
               xcntr(iblk) = 0.d0
              endif
            endif
            rad = radius(iblk)
            xyz(1:3) = circle_xyz(0.d0,0.0d0*pi,rad)
            xmin(iblk) = xyz(1)
            xyz(1:3) = circle_xyz(1.d0,0.0d0*pi,rad)
            xmax(iblk) = xyz(1)
            ymax(iblk) = xyz(2)
            xyz(1:3) = circle_xyz(1.d0,0.5d0*pi,rad)
            zmax(iblk) = xyz(3)
            xyz(1:3) = circle_xyz(1.d0,1.0d0*pi,rad)
            ymin(iblk) = xyz(2)
            xyz(1:3) = circle_xyz(1.d0,1.5d0*pi,rad)
            zmin(iblk) = xyz(3)
          endif
        enddo
!
!---------------------------------------------------------------------
! Position of tip of the nose
!---------------------------------------------------------------------
       x_nose=minval(xmin(1:nblk)+xcntr(1:nblk))
       xmin(1:nblk) =(xmin(1:nblk)+xcntr(1:nblk)) - x_nose
       xmax(1:nblk) =(xmax(1:nblk)+xcntr(1:nblk)) - x_nose
!
!----------------------------------------------------------------------
! Read reference length and reference area
!----------------------------------------------------------------------
        read(*,*)
        read(*,*)len_ref
        read(*,*)area_ref
!
        if(len_ref .lt. 0) then
           len_ref = maxval(xmax(1:nblk))-minval(xmin(1:nblk))
        endif
        if(area_ref .lt. 0) then
           area_ref = pi*maxval(zmax(1:nblk))**2
        endif
!
!---------------------------------------------------------------------
! Print summary 
!---------------------------------------------------------------------
        write(*,*) " "
        write(*,115)"========","============","============"
     &             ,"==============","===============","==============="
     &             ,"==============","==============","=============="
        write(*,*) " SUMMARY OF GEOMETRIES"
        write(*,115)"========","============","============"
     &             ,"==============","===============","==============="
     &             ,"==============","==============","=============="
        write(*,112) " Block # "," Shape      "," Position   "
     &            ," R [m]      "," angle[deg]"," h [m]      "
     &            ," x0 [m]     "," flap L [m]"," flap H [m] "
        write(*,115) " --------"," ---------- "," ---------- "
     &              ,"------------","------------","------------"
     &              ,"------------","------------","------------"
        do i=1,nblk
          write(*,116) i,trim(gtype(i)),trim(ptype(i))
     &                 ,radius(i),ahalf(i),xlngt(i)
     &                 ,xcntr(i)-x_nose,flngt(i),fht(i)
        enddo
!
        write(*,111)"========"
     &             ,"==============","===============","==============="
     &             ,"==============","===============","==============="
     &             ,"=============="
        write(*,112) " Block # "
     &            ," Param1 min "," Param1 max "," max - min  "
     &            ," Param2 min "," Param2 max "," max - min  "
     &            ," Area       "
        write(*,112) "--------"
     &              ,"------------","------------","------------"
     &              ,"------------","------------","------------"
     &              ,"------------"
        do i=1,nblk
          write(*,113)i,v1min(i),v1max(i),v1max(i)-v1min(i)
     &                 ,v2min(i),v2max(i),v2max(i)-v2min(i)
     &                 ,(v1max(i)-v1min(i))*(v2max(i)-v2min(i))
        enddo
!
 111    format(2x,a8,9a14)
 112    format(2x,a8,9(x,a12,x))
 113    format(2x,i4,2x,9(2x,f12.8))  
 114    format(2x,a6,9(2x,f12.8))  
 115    format(2x,a8,2a12,6a14)
 116    format(2x,i4,2x,2(x,a10,x),6(2x,f12.6))  
!
        write(*,111)"========"
     &             ,"==============","===============","==============="
     &             ,"==============","===============","==============="
     &             ,"==============","===============","==============="
        write(*,112) " Block # "
     &            ,"   x min    ","   x max   ","   x len   "
     &            ,"   y min    ","   y max   ","   y len   "
     &            ,"   z min    ","   z max   ","   z len   "
        write(*,112) "--------"
     &              ,"------------","------------","------------"
     &              ,"------------","------------","------------"
     &              ,"------------","------------","------------"
        do i=1,nblk
          write(*,113)i,xmin(i),xmax(i),xmax(i)-xmin(i)
     &                 ,ymin(i),ymax(i),ymax(i)-ymin(i)
     &                 ,zmin(i),zmax(i),zmax(i)-zmin(i)
        enddo
        write(*,112) "--------"
     &              ,"------------","------------","------------"
     &              ,"------------","------------","------------"
     &              ,"------------","------------","------------"
        write(*,114) " Total   "
     &               ,minval(xmin(1:nblk)),maxval(xmax(1:nblk))
     &               ,maxval(xmax(1:nblk))-minval(xmin(1:nblk))
     &               ,minval(ymin(1:nblk)),maxval(ymax(1:nblk))
     &               ,maxval(ymax(1:nblk))-minval(ymin(1:nblk))
     &               ,minval(zmin(1:nblk)),maxval(zmax(1:nblk))
     &               ,maxval(zmax(1:nblk))-minval(zmin(1:nblk))
        write(*,111)"========"
     &             ,"==============","===============","==============="
     &             ,"==============","===============","==============="
     &             ,"==============","===============","==============="
        write(*,*) " "
        print *, "  Reference area:", area_ref, "[Deg]"
        print *, "  Reference length:", len_ref, "[Deg]"
        write(*,111)"========"
     &             ,"==============","==============="
!
        return
!         
        end subroutine set_geom_params
!
        subroutine print_err(iblk0,iblk1)
           integer          :: iblk0,iblk1
           character(len=2) :: blk0,blk1
           write (blk0,'(i2)') iblk0
           write (blk1,'(i2)') iblk1
           print *, " "
           print *, "Sorry! "
           print *, "Connection between "
           print *, "  Block #"
     &            ,  trim(blk0)," "
     &            ,  trim(gtype(iblk0))," "
     &            ,  trim(ptype(iblk0))
           print *, "and "
           print *, "  Block #"
     &            ,  trim(blk1)," "
     &            ,  trim(gtype(iblk1))," "
     &            ,  trim(ptype(iblk1))
           print *, "is not supported."
           call exit
        end subroutine
!
      end module geoms
!
!**********************************************************************
! module mcnew
!**********************************************************************
      module mcnew
       use geoms
       use funcs, only: geom_vnor
     &                 ,rot_body
       implicit none
!
       contains
!
       subroutine mcnew_aero(nsmpb,ntry,iout
     &                      ,alp,bet,cg,cp0
     &                      ,aeroc)
       implicit none
! in and out
        integer,intent(in) :: nsmpb(:)! Number of points on blocks
     &                       ,ntry    ! Number of trials
     &                       ,iout    ! Interval for intermediate output
        real(8) :: alp   ! pitch angle in radian
     &            ,bet   ! yaw angle in radian
        real(8),intent(in) :: cg(3)     ! center of gravity
        real(8),intent(in) :: cp0       ! Max pressure for MN
        real(8),intent(out) :: aeroc(13)
                 ! aeroc(1): CD , Drag coef
                 ! aeroc(2): CLy, Lift coef in y-direction
                 ! aeroc(3): CLz, Lift coef in y-direction
                 ! aeroc(4): L/D, Lift-to-drag ratio
                 ! aeroc(5): CA , Axial force coef
                 ! aeroc(6): CNy, Normal force coef in y
                 ! aeroc(7): CNz, Normal force coef in z
                 ! aeroc(8): Cmx, Moment coef aroud x-axis passing c.g.
                 ! aeroc(9): Cmy, Moment coef aroud y-axis passing c.g.
                 ! aeroc(10):Cmz, Moment coef aroud z-axis passing c.g.
                 ! aeroc(11):Cm0x, Moment coef aroud x-axis
                 ! aeroc(12):Cm0y, Moment coef aroud y-axis
                 ! aeroc(13):Cm0z, Moment coef aroud z-axis
!
! variables and parameters related to geometries
       real(8) :: xyz(3),vnor(4),unor(4),unor_rot(3)
       real(8) :: var1,var2
       real(8) :: sgn
       real(8) :: xyz0(3),xyzcg(3)
       real(8) :: param(5)
! variables to calulate aerodynamic coefficients
       real(8),allocatable:: cd_blk(:),cd_blk_sum(:),cd_blk_ave(:)
       real(8),allocatable:: cl_blk(:),cl_blk_sum(:),cl_blk_ave(:)
       real(8),allocatable:: cly_blk(:),cly_blk_sum(:),cly_blk_ave(:)
       real(8),allocatable:: clz_blk(:),clz_blk_sum(:),clz_blk_ave(:)
       real(8),allocatable:: ca_blk(:),ca_blk_sum(:),ca_blk_ave(:)
       real(8),allocatable:: cny_blk(:),cny_blk_sum(:),cny_blk_ave(:)
       real(8),allocatable:: cnz_blk(:),cnz_blk_sum(:),cnz_blk_ave(:)
       real(8),allocatable:: cm0x_blk(:),cm0x_blk_sum(:),cm0x_blk_ave(:)
       real(8),allocatable:: cm0y_blk(:),cm0y_blk_sum(:),cm0y_blk_ave(:)
       real(8),allocatable:: cm0z_blk(:),cm0z_blk_sum(:),cm0z_blk_ave(:)
       real(8),allocatable:: cmx_blk(:),cmx_blk_sum(:),cmx_blk_ave(:)
       real(8),allocatable:: cmy_blk(:),cmy_blk_sum(:),cmy_blk_ave(:)
       real(8),allocatable:: cmz_blk(:),cmz_blk_sum(:),cmz_blk_ave(:)
       real(8),allocatable:: area_blk(:),area_blk_sum(:),area_blk_ave(:)
       real(8),allocatable:: area_yz_blk(:),area_yz_blk_sum(:)
     &                      ,area_yz_blk_ave(:)
       real(8),allocatable:: area_intg(:)
       integer,allocatable:: ivsbl(:),ivsbl_blk_sum(:)
       real(8),allocatable:: cd_vari(:),cl_vari(:),cly_vari(:),
     &      clz_vari(:),ca_vari(:),cny_vari(:),cnz_vari(:),cm0x_vari(:),
     &      cm0y_vari(:),cm0z_vari(:),cmx_vari(:),cmy_vari(:),
     &      cmz_vari(:)
              
       real(8) :: fd,cd,cd_sum,cd_ave,cd_tot
       real(8) :: ca_tot,cny_tot,cnz_tot,cn_tot
       real(8) :: force_xyz(3)
       real(8) :: force_body(3)
       real(8) :: cl,cl_sum,cl_ave,cl_tot
       real(8) :: cly,cly_sum,cly_ave,cly_tot
       real(8) :: clz,clz_sum,clz_ave,clz_tot
       real(8) :: cm0x,cm0x_sum,cm0x_ave,cm0x_tot
       real(8) :: cm0y,cm0y_sum,cm0y_ave,cm0y_tot
       real(8) :: cm0z,cm0z_sum,cm0z_ave,cm0z_tot
       real(8) :: cmx,cmx_sum,cmx_ave,cmx_tot
       real(8) :: cmy,cmy_sum,cmy_ave,cmy_tot
       real(8) :: cmz,cmz_sum,cmz_ave,cmz_tot
       real(8) :: cpl,cp_xyz(3),ca_xyz(3)
       real(8) :: amomt_cg_loc(3),amomt_cg(3)
       real(8) :: amomt_0_loc(3) ,amomt_0(3)
       real(8) :: area_visible,area_tot=0.d0
       real(8) :: area_yz_visible,area_yz_tot=0.d0
       real(8) :: jacob
       real(8) :: cd_var,cl_var,cly_var,clz_var,ca_var,cny_var,cnz_var,
     &      cm0x_var,cm0y_var,cm0z_var,cmx_var,cmy_var,cmz_var
       
! variables and parameters for loop control
       integer :: i,iblk,itry
       integer :: nsmp
!
! variables and parameters for random number generation
       real(8) :: rand1,rand2
       integer :: seedsize=30
       integer, allocatable :: seed(:)
!
!     _________________________________________________________________
!     Initialization of arrays
!     -----------------------------------------------------------------
       allocate(cd_blk(nblk),cd_blk_sum(nblk),cd_blk_ave(nblk))
       allocate(cl_blk(nblk),cl_blk_sum(nblk),cl_blk_ave(nblk))
       allocate(cly_blk(nblk),cly_blk_sum(nblk),cly_blk_ave(nblk))
       allocate(clz_blk(nblk),clz_blk_sum(nblk),clz_blk_ave(nblk))
       allocate(ca_blk(nblk),ca_blk_sum(nblk),ca_blk_ave(nblk))
       allocate(cny_blk(nblk),cny_blk_sum(nblk),cny_blk_ave(nblk))
       allocate(cnz_blk(nblk),cnz_blk_sum(nblk),cnz_blk_ave(nblk))
       allocate(cmx_blk(nblk),cmx_blk_sum(nblk),cmx_blk_ave(nblk))
       allocate(cmy_blk(nblk),cmy_blk_sum(nblk),cmy_blk_ave(nblk))
       allocate(cmz_blk(nblk),cmz_blk_sum(nblk),cmz_blk_ave(nblk))
       allocate(cm0x_blk(nblk),cm0x_blk_sum(nblk),cm0x_blk_ave(nblk))
       allocate(cm0y_blk(nblk),cm0y_blk_sum(nblk),cm0y_blk_ave(nblk))
       allocate(cm0z_blk(nblk),cm0z_blk_sum(nblk),cm0z_blk_ave(nblk))
       allocate(area_blk(nblk),area_blk_sum(nblk),area_blk_ave(nblk))
       allocate(area_yz_blk(nblk),area_yz_blk_sum(nblk)
     &                           ,area_yz_blk_ave(nblk))
       allocate(area_intg(nblk))
       allocate(ivsbl(nblk),ivsbl_blk_sum(nblk))
       allocate(cd_vari(itry),cl_vari(itry),cly_vari(itry),
     &      clz_vari(itry),ca_vari(itry),cny_vari(itry),cnz_vari(itry),
     &      cmx_vari(itry),cmy_vari(itry),cmz_vari(itry),cm0x_vari(itry)
     &      ,cm0y_vari(itry),cm0z_vari(itry))
       
       cd_blk_sum(:) = 0.d0
       cl_blk_sum(:) = 0.d0
       cly_blk_sum(:) = 0.d0
       clz_blk_sum(:) = 0.d0
       ca_blk_sum(:) = 0.d0
       cny_blk_sum(:) = 0.d0
       cnz_blk_sum(:) = 0.d0
       cm0x_blk_sum(:) = 0.d0
       cm0y_blk_sum(:) = 0.d0
       cm0z_blk_sum(:) = 0.d0
       cmx_blk_sum(:) = 0.d0
       cmy_blk_sum(:) = 0.d0
       cmz_blk_sum(:) = 0.d0
       area_blk_sum(:) = 0.d0
       ivsbl_blk_sum(:) = 0
       cd_tot = 0.d0
       cl_tot = 0.d0
       area_tot= 0.d0
       cd_vari(:)  = 0.d0 
       cly_vari(:) = 0.d0 
       clz_vari(:) = 0.d0 
       cl_vari(:)  = 0.d0 
       ca_vari(:)  = 0.d0 
       cny_vari(:) = 0.d0 
       cnz_vari(:) = 0.d0 
       cm0x_vari(:) = 0.d0
       cm0y_vari(:) = 0.d0
       cm0z_vari(:) = 0.d0
       cmx_vari(:) = 0.d0 
       cmy_vari(:) = 0.d0 
       cmz_vari(:) = 0.d0

!
!      _________________________________________________________________
!      Print header line for intermediate reports
!      -----------------------------------------------------------------
       write(*,*)
       write(*,*)
       write(*,*)
       write(*,"(2x,3(a40))")"****************************************"
     &                      ,"********  << Iteration start >>  *******"
     &                      ,"****************************************"
       write(*,117)" Trial  ","Local step values","Cumulative average"
       write(*,"(2x,a10,2(x,a54,x))")"-----------"
     &    ,"-----------------------------------------------------------"
     &    ,"-----------------------------------------------------------"
       write(*,118)" number ","ivisible","CD","CLy","CLz"
     &                       ,"ivisible","CD","CLy","CLz"
       write(*,"(2x,a10,8(x,a12,x))")"-----------"
     &                 ,"------------","------------"
     &                 ,"------------","------------"
     &                 ,"------------","------------"
     &                 ,"------------","------------"
!
       open(21,file="convergence.dat")
       write(21,117)"iTrial  ","Local step values","Cumulative average"
       write(21,118)"iTrial  "
     &     ,"ivisible","CD","CLy","CLz"
     &     ,"ivisible","CD","CLy","CLz"
 117   format('#',x,a7,7x,2(a28,28x))
 118   format('#',x,a7,7x,2(a12,2x,3(x,a6,7x)))
!
!      _________________________________________________________________
!      Open files to save geometry point data
!      -----------------------------------------------------------------
       open(22,file="visible_points.dat")
       open(23,file="hidden_points.dat")
!
!      ________________________________________________________________
!      Preparaion of random seed (Activate when necessary)
!      ----------------------------------------------------------------
!      call random_seed(size=seedsize)
!      allocate(seed(seedsize))
!      do i = 1, seedsize
!         call system_clock(count=seed(i))
!      end do
!      call random_seed(put=seed(:))
       call random_seed
!
!_______________________________________________________________________
! main loop
!-----------------------------------------------------------------------
!
! trial loop start
       do itry=1,ntry
!
! reset arrys
       cd_blk(:)   = 0.d0
       cl_blk(:)   = 0.d0
       cly_blk(:)  = 0.d0
       clz_blk(:)  = 0.d0
       ca_blk(:)   = 0.d0
       cny_blk(:)  = 0.d0
       cnz_blk(:)  = 0.d0
       cm0x_blk(:)  = 0.d0
       cm0y_blk(:)  = 0.d0
       cm0z_blk(:)  = 0.d0
       cmx_blk(:)  = 0.d0
       cmy_blk(:)  = 0.d0
       cmz_blk(:)  = 0.d0
       area_blk(:) = 0.d0
       area_yz_blk(:) = 0.d0
       ivsbl(:)    = 0

       area_intg(1:nblk) = (v1max(1:nblk)-v1min(1:nblk))
     &                    *(v2max(1:nblk)-v2min(1:nblk))
!
!      ________________________________________________________________
!      Loop over blocks to calculate forces on each blocks
!      ----------------------------------------------------------------
       do iblk = 1,nblk
!
         force_xyz(1:3)  = 0.d0
         force_body(1:3) = 0.d0
         amomt_cg(1:3)   = 0.d0
         amomt_0(1:3)    = 0.d0
         param(1:5)      = 0.d0
!        
!        ------------------------------------------------------------
!        set parameters
!        ------------------------------------------------------------
           if(gtype(iblk)=="sphere") then
             param(1) = radius(iblk)
             param(2) = 0.d0
             param(3) = 0.d0
             param(4) = 0.d0
             param(5) = 0.d0
           elseif(gtype(iblk)=="cone") then
             if(ptype(iblk)=="fore") sgn = 1.d0
             if(ptype(iblk)=="rear") sgn =-1.d0
             param(1) = radius(iblk)
             param(2) = xlngt(iblk)
             param(3) = sgn
             param(4) = 0.d0
             param(5) = 0.d0
          elseif(gtype(iblk)=="flap") then
             if(ptype(iblk)=="fore") sgn = 1.d0
             if(ptype(iblk)=="rear") sgn =-1.d0
             param(1) = radius(iblk)
             param(2) = xlngt(iblk)
             param(3) = sgn
             param(4) = flngt(iblk)
             param(5) = fht(iblk)
           elseif(gtype(iblk)=="cylinder") then
             if(ptype(iblk)=="vertical") sgn =  1.d0
             if(ptype(iblk)=="horizon")  sgn = -1.d0
             param(1) = radius(iblk)
             param(2) = xlngt(iblk)
             param(3) = sgn
             param(4) = 0.d0
             param(5) = 0.d0
           elseif(gtype(iblk)=="circle") then
             if(ptype(iblk)=="top") sgn  = -1.d0
             if(ptype(iblk)=="bottom") sgn =  1.d0
             param(1) = radius(iblk)
             param(2) = xlngt(iblk)
             param(3) = sgn
             param(4) = 0.d0
             param(5) = 0.d0
           elseif(gtype(iblk)=="shoulder") then
             param(1) = radius(iblk)
             param(2) = radcen(iblk)
             param(3) = 0.d0
             param(4) = 0.d0
             param(5) = 0.d0
           endif
!        ______________________________________________________________
!        Loop over sample points
!        --------------------------------------------------------------
         nsmp = nsmpb(iblk)
         do i=1,nsmp
!          ------------------------------------------------------------
!          Random sampling
!          ------------------------------------------------------------
           call random_number(rand1)
           call random_number(rand2)
           var1 = (1.d0-rand1)*v1min(iblk) + rand1*v1max(iblk)
           var2 = (1.d0-rand2)*v2min(iblk) + rand2*v2max(iblk)
!
!          ------------------------------------------------------------
!          xyz coordinate and normal vector
!          ------------------------------------------------------------
           vnor(1:4) = geom_vnor(xyz,gtype(iblk),var1,var2,param)
!
!          ------------------------------------------------------------
!          Unit normal vector
!          ------------------------------------------------------------
           unor(1) = vnor(1)/vnor(4)
           unor(2) = vnor(2)/vnor(4)
           unor(3) = vnor(3)/vnor(4)
           unor(4) = vnor(4)
!
!          ------------------------------------------------------------
!          Rotational transformation of the normal vector
!          ------------------------------------------------------------
           unor_rot(1:3) = rot_body(unor(1:3),alp,bet)
!
!          ------------------------------------------------------------
!          Forces on visible points 
!          ------------------------------------------------------------
           if(unor_rot(1).lt.0.d0)then
!           -----------------------------------------------------------
!           Local force coef components at the current sample point
!           -----------------------------------------------------------
            cpl  =  2.d0 * unor_rot(1)**2.0
            cp_xyz(1:3)  = -cpl*unor_rot(1:3)
            ca_xyz(1:3)  = -cpl*unor(1:3)
!           -----------------------------------------------------------
!           Moment components around origin at the current sample point
!           -----------------------------------------------------------
            xyz0(1) = xyz(1)+xcntr(iblk)-x_nose
            xyz0(2) = xyz(2)
            xyz0(3) = xyz(3)
            amomt_0_loc(1)  = xyz0(2) *ca_xyz(3)-xyz0(3) *ca_xyz(2)
            amomt_0_loc(2)  = xyz0(3) *ca_xyz(1)-xyz0(1) *ca_xyz(3)
            amomt_0_loc(3)  = xyz0(1) *ca_xyz(2)-xyz0(2) *ca_xyz(1)
!           -----------------------------------------------------------
!           Moment components around c.g. at the current sample point
!           -----------------------------------------------------------
            xyzcg(1) =(xyz(1)+xcntr(iblk)-x_nose)-cg(1)
            xyzcg(2) = xyz(2)-cg(2)
            xyzcg(3) = xyz(3)-cg(3)
            amomt_cg_loc(1) = xyzcg(2)*ca_xyz(3)-xyzcg(3)*ca_xyz(2)
            amomt_cg_loc(2) = xyzcg(3)*ca_xyz(1)-xyzcg(1)*ca_xyz(3)
            amomt_cg_loc(3) = xyzcg(1)*ca_xyz(2)-xyzcg(2)*ca_xyz(1)
!           ------------------------------------------------------------
!           Summing up the forces over the sample points
!           ------------------------------------------------------------
            force_xyz(1:3)  = force_xyz(1:3)  + cp_xyz(1:3)*vnor(4)
            force_body(1:3) = force_body(1:3) + ca_xyz(1:3)*vnor(4)
            amomt_cg(1:3)   = amomt_cg(1:3) 
     &                      + amomt_cg_loc(1:3)*vnor(4)
            amomt_0(1:3)    = amomt_0(1:3)
     &                      + amomt_0_loc(1:3)*vnor(4)
!           ------------------------------------------------------------
!           Summing up the aera
!           ------------------------------------------------------------
            area_blk(iblk)    = area_blk(iblk) + vnor(4)
            area_yz_blk(iblk) = area_yz_blk(iblk) - unor_rot(1)* vnor(4)
!           ------------------------------------------------------------
!           Counting visible points
!           ------------------------------------------------------------
            ivsbl(iblk) = ivsbl(iblk) + 1
!           ------------------------------------------------------------
!           Write out the visible points (Only once)
!           ------------------------------------------------------------
            if(itry==1)write(22,119)
     &               xyz(1)+xcntr(iblk)-x_nose,xyz(2:3),cpl,cp_xyz(1:3)
!
!          ------------------------------------------------------------
!          Invisible points 
!          ------------------------------------------------------------
           else
            cpl  = 0.d0
            cp_xyz(1:3)  = 0.d0
            ca_xyz(1:3)  = 0.d0
!           ------------------------------------------------------------
!           Write out the invisible points (Only once)
!           ------------------------------------------------------------
            if(itry==1)write(23,119)
     &               xyz(1)+xcntr(iblk)-x_nose,xyz(2:3),cpl,cp_xyz(1:3)
           endif
         enddo 
 119    format(7e14.6)
!        ___End of loop over sample points____________________________
!
!        --------------------------------------------------------------
!        Visible area
!        --------------------------------------------------------------
         area_visible = area_blk(iblk)/dble(nsmp)*area_intg(iblk)
!
!        --------------------------------------------------------------
!        Visible area projected on yz-plane
!        --------------------------------------------------------------
         area_yz_visible = area_yz_blk(iblk)/dble(nsmp)*area_intg(iblk)
!
!        --------------------------------------------------------------
!        Aerodynamic coeffs of the current block at the current trial
!        --------------------------------------------------------------
         jacob = area_intg(iblk)/dble(nsmp)/area_ref
         cd_blk(iblk)  = force_xyz(1)*jacob
         cly_blk(iblk) = force_xyz(2)*jacob
         clz_blk(iblk) = force_xyz(3)*jacob
         cl_blk(iblk)  = dsqrt(force_xyz(2)**2+force_xyz(3)**2)*jacob
         ca_blk(iblk)  = force_body(1)*jacob
         cny_blk(iblk) = force_body(2)*jacob
         cnz_blk(iblk) = force_body(3)*jacob
         cm0x_blk(iblk) = amomt_0(1)*jacob/len_ref
         cm0y_blk(iblk) = amomt_0(2)*jacob/len_ref
         cm0z_blk(iblk) = amomt_0(3)*jacob/len_ref
         cmx_blk(iblk)  = amomt_cg(1)*jacob/len_ref
         cmy_blk(iblk)  = amomt_cg(2)*jacob/len_ref
         cmz_blk(iblk)  = amomt_cg(3)*jacob/len_ref
!
!        ---------------------------------------------------------------
!        Summing up aerodynamic coeffs and areas for averaging
!        --------------------------------------------------------------
         cd_blk_sum(iblk)   = cd_blk_sum(iblk)   + cd_blk(iblk)
         cly_blk_sum(iblk)  = cly_blk_sum(iblk)  + cly_blk(iblk)
         clz_blk_sum(iblk)  = clz_blk_sum(iblk)  + clz_blk(iblk)
         cl_blk_sum(iblk)   = cl_blk_sum(iblk)   + cl_blk(iblk)
         ca_blk_sum(iblk)   = ca_blk_sum(iblk)   + ca_blk(iblk)
         cny_blk_sum(iblk)  = cny_blk_sum(iblk)  + cny_blk(iblk)
         cnz_blk_sum(iblk)  = cnz_blk_sum(iblk)  + cnz_blk(iblk)
         cm0x_blk_sum(iblk) = cm0x_blk_sum(iblk) + cm0x_blk(iblk)
         cm0y_blk_sum(iblk) = cm0y_blk_sum(iblk) + cm0y_blk(iblk)
         cm0z_blk_sum(iblk) = cm0z_blk_sum(iblk) + cm0z_blk(iblk)
         cmx_blk_sum(iblk)  = cmx_blk_sum(iblk)  + cmx_blk(iblk)
         cmy_blk_sum(iblk)  = cmy_blk_sum(iblk)  + cmy_blk(iblk)
         cmz_blk_sum(iblk)  = cmz_blk_sum(iblk)  + cmz_blk(iblk)
         area_blk_sum(iblk)    = area_blk_sum(iblk)    + area_visible
         area_yz_blk_sum(iblk) = area_yz_blk_sum(iblk) + area_yz_visible
         ivsbl_blk_sum(iblk)  = ivsbl_blk_sum(iblk)  + ivsbl(iblk)

!        ---------------------------------------------------------------
!        Prepare to get variance in aerodynamic coeffs of the whole body
!     --------------------------------------------------------------

         print *, "Preparing to calculate variances"
         cd_vari(itry)    = cd_vari(itry)   + cd_blk(iblk)
         cly_vari(itry)   = cly_vari(itry)  + cly_blk(iblk) 
         clz_vari(itry)   = clz_vari(itry)  + clz_blk(iblk) 
         cl_vari(itry)    = cl_vari(itry)   + cl_blk(iblk)  
         ca_vari(itry)    = ca_vari(itry)   + ca_blk(iblk)  
         cny_vari(itry)   = cny_vari(itry)  + cny_blk(iblk) 
         cnz_vari(itry)   = cnz_vari(itry)  + cnz_blk(iblk) 
         cm0x_vari(itry)  = cm0x_vari(itry) + cm0x_blk(iblk)
         cm0y_vari(itry)  = cm0y_vari(itry) + cm0y_blk(iblk)
         cm0z_vari(itry)  = cm0z_vari(itry) + cm0z_blk(iblk)
         cmx_vari(itry)   = cmx_vari(itry)  + cmx_blk(iblk) 
         cmy_vari(itry)   = cmy_vari(itry)  + cmy_blk(iblk) 
         cmz_vari(itry)   = cmz_vari(itry)  + cmz_blk(iblk)
         print *, "Calculated all variances"

       enddo
!      ___End of loop over blocks______________________________________
!
!      ---------------------------------------------------------------
!      Aerodynamic coeffs of the whole body
!      --------------------------------------------------------------
       cd_tot    = sum(cd_blk(1:nblk))
       cly_tot   = sum(cly_blk(1:nblk))
       clz_tot   = sum(clz_blk(1:nblk))
       cl_tot    = sum(cl_blk(1:nblk))
       ca_tot    = sum(ca_blk(1:nblk))
       cny_tot   = sum(cny_blk(1:nblk))
       cnz_tot   = sum(cnz_blk(1:nblk))
       cm0x_tot  = sum(cm0x_blk(1:nblk))
       cm0y_tot  = sum(cm0y_blk(1:nblk))
       cm0z_tot  = sum(cm0z_blk(1:nblk))
       cmx_tot   = sum(cmx_blk(1:nblk))
       cmy_tot   = sum(cmy_blk(1:nblk))
       cmz_tot   = sum(cmz_blk(1:nblk))
!
!      ---------------------------------------------------------------
!      Aerodynamic coeffs of the whole body
!      --------------------------------------------------------------
        if(mod(itry,iout).eq.0) then
          write(*,120)itry
     &          ,sum(ivsbl(1:nblk)),cd_tot,cly_tot,clz_tot
     &          ,sum(ivsbl_blk_sum(1:nblk))
     &          ,sum(cd_blk_sum(1:nblk))/dble(itry)
     &          ,sum(cly_blk_sum(1:nblk))/dble(itry)
     &          ,sum(clz_blk_sum(1:nblk))/dble(itry)
          write(21,120)itry
     &           ,sum(ivsbl(1:nblk)),cd_tot,cly_tot,clz_tot
     &          ,sum(ivsbl_blk_sum(1:nblk))
     &          ,sum(cd_blk_sum(1:nblk))/dble(itry)
     &          ,sum(cly_blk_sum(1:nblk))/dble(itry)
     &          ,sum(clz_blk_sum(1:nblk))/dble(itry)
 120      format(i12,2x,2(i12,2x,3(f12.8,2x)))
        endif
!
       enddo
       write(*,"(2x,3(a40))")"****************************************"
     &                      ,"*******  << End of iteration >>  *******"
     &                      ,"****************************************"
       write(*,*)
       write(*,*)
       write(*,*)
!
       close(21)
       close(22)
       close(23)
!
!      ___End of loop over trials______________________________________
!
! End of main loop
!======================================================================
! Post processing
!======================================================================
!
!     --------------------------------------------------------------
!     Block aerodynamic coeffs and areas averaged over trials
!     -------------------------------------------------------------
      cd_blk_ave(1:nblk)   = cd_blk_sum(1:nblk)/dble(ntry)
      cl_blk_ave(1:nblk)   = cl_blk_sum(1:nblk)/dble(ntry)
      cly_blk_ave(1:nblk)  = cly_blk_sum(1:nblk)/dble(ntry)
      clz_blk_ave(1:nblk)  = clz_blk_sum(1:nblk)/dble(ntry)
      ca_blk_ave(1:nblk)   = ca_blk_sum(1:nblk)/dble(ntry)
      cny_blk_ave(1:nblk)  = cny_blk_sum(1:nblk)/dble(ntry)
      cnz_blk_ave(1:nblk)  = cnz_blk_sum(1:nblk)/dble(ntry)
      cm0x_blk_ave(1:nblk) = cm0x_blk_sum(1:nblk)/dble(ntry)
      cm0y_blk_ave(1:nblk) = cm0y_blk_sum(1:nblk)/dble(ntry)
      cm0z_blk_ave(1:nblk) = cm0z_blk_sum(1:nblk)/dble(ntry)
      cmx_blk_ave(1:nblk)  = cmx_blk_sum(1:nblk)/dble(ntry)
      cmy_blk_ave(1:nblk)  = cmy_blk_sum(1:nblk)/dble(ntry)
      cmz_blk_ave(1:nblk)  = cmz_blk_sum(1:nblk)/dble(ntry)
      area_blk_ave(1:nblk) = area_blk_sum(1:nblk)/dble(ntry)
      area_yz_blk_ave(1:nblk) = area_yz_blk_sum(1:nblk)/dble(ntry)
!
!     --------------------------------------------------------------
!     Aerodynamic coeffs for the whole body
!     -------------------------------------------------------------
       cd_tot   = sum(cd_blk_ave(1:nblk))
       cly_tot  = sum(cly_blk_ave(1:nblk))
       clz_tot  = sum(clz_blk_ave(1:nblk))
       cl_tot   = dsqrt(cly_tot**2+clz_tot**2)
       cm0x_tot = sum(cm0x_blk_ave(1:nblk))
       cm0y_tot = sum(cm0y_blk_ave(1:nblk))
       cm0z_tot = sum(cm0z_blk_ave(1:nblk))
       cmx_tot  = sum(cmx_blk_ave(1:nblk))
       cmy_tot  = sum(cmy_blk_ave(1:nblk))
       cmz_tot  = sum(cmz_blk_ave(1:nblk))
       ca_tot   = sum(ca_blk_ave(1:nblk))
       cny_tot  = sum(cny_blk_ave(1:nblk))
       cnz_tot  = sum(cnz_blk_ave(1:nblk))
       cn_tot   = dsqrt(cny_tot**2+cnz_tot**2)
       area_tot = sum(area_blk_ave(1:nblk))
       area_yz_tot = sum(area_yz_blk_ave(1:nblk))
!
!     ---------------------------------------------------------------
!     Variance in aerodynamic coeffs of the whole body
!     --------------------------------------------------------------

       do itry=1, ntry
          cd_vari(itry)    = (cd_vari(itry)  -cd_tot   )**2
          cly_vari(itry)   = (cly_vari(itry) -cly_tot  )**2
          clz_vari(itry)   = (clz_vari(itry) -clz_tot  )**2
          cl_vari(itry)    = (cl_vari(itry)  -cl_tot   )**2
          ca_vari(itry)    = (ca_vari(itry)  -ca_tot   )**2
          cny_vari(itry)   = (cny_vari(itry) -cny_tot  )**2
          cnz_vari(itry)   = (cnz_vari(itry) -cnz_tot  )**2
          cm0x_vari(itry)  = (cm0x_vari(itry)-cm0x_tot )**2
          cm0y_vari(itry)  = (cm0y_vari(itry)-cm0y_tot )**2
          cm0z_vari(itry)  = (cm0z_vari(itry)-cm0z_tot )**2
          cmx_vari(itry)   = (cmx_vari(itry) -cmx_tot  )**2
          cmy_vari(itry)   = (cmy_vari(itry) -cmy_tot  )**2
          cmz_vari(itry)   = (cmz_vari(itry) -cmz_tot  )**2
       enddo
       
       cd_var = sum(cd_vari(1:ntry)) / (ntry - 1)
       cly_var = sum(cly_vari(1:ntry)) / (ntry - 1)
       clz_var = sum(clz_vari(1:ntry)) / (ntry - 1)
       cl_var = sum(cl_vari(1:ntry)) / (ntry - 1)
       ca_var = sum(ca_vari(1:ntry)) / (ntry - 1)
       cny_var = sum(cny_vari(1:ntry)) / (ntry - 1)
       cnz_var = sum(cnz_vari(1:ntry)) / (ntry - 1)
       cm0x_var = sum(cm0x_vari(1:ntry)) / (ntry - 1)
       cm0y_var = sum(cm0y_vari(1:ntry)) / (ntry - 1)
       cm0z_var = sum(cm0z_vari(1:ntry)) / (ntry - 1)
       cmx_var = sum(cmx_vari(1:ntry)) / (ntry - 1)
       cmy_var = sum(cmy_vari(1:ntry)) / (ntry - 1)
       cmz_var = sum(cmz_vari(1:ntry)) / (ntry - 1)

       print *, cd_var
       print *, cly_var
       print *, clz_var 
       print *, cl_var
       print *, ca_var
       print *, cny_var 
       print *, cnz_var 
       print *, cm0x_var
       print *, cm0y_var
       print *, cm0z_var
       print *, cmx_var 
       print *, cmy_var 
       print *, cmz_var 
       
!     --------------------------------------------------------------
!     Store aerodynamic coefficients in output array
!     -------------------------------------------------------------
       aeroc(1)  = cd_tot
       aeroc(2)  = cly_tot
       aeroc(3)  = clz_tot
       aeroc(4)  = clz_tot/cd_tot
       aeroc(5)  = ca_tot
       aeroc(6)  = cny_tot
       aeroc(7)  = cnz_tot
       aeroc(8)  = cmx_tot
       aeroc(9)  = cmy_tot
       aeroc(10) = cmz_tot
       aeroc(11) = cm0x_tot
       aeroc(12) = cm0y_tot
       aeroc(13) = cm0z_tot
!
!     -----------------------------------------------------------------
!     Disply summary of aerodynamic coefficients
!     -----------------------------------------------------------------
 121    format(2x,3a38) 
 122    format(2x,3a34) 
       write(*,*)
       write(*,121)"=========================================="
     &            
     &           ,"=========================================="
     &           ,"=========================================="
       write(*,*) " SUMMARY "
       write(*,121)"=========================================="
     &            ,"=========================================="
     &            ,"=========================================="
       write(*,*) " Aerodynamic coefficients of each blocks           "
       write(*,121)
     &            "------------------------------------------"
     &           ,"------------------------------------------"
     &           ,"------------------------------------------"
       write(*,"(x,2a9,7(2x,a8,4x))") " Block # ","Shape "
     &                       ,"CD_av ","CLy_av","CLz_av"
     &                       ,"L/D   "
     &                       ,"CA_av ","CNy_av","CNz_av"
       write(*,"(x,2a9,x,7(x,a12,x))") " ------- ","----------"
     &                 ,"------------","------------","------------"
     &                 ,"------------","------------","------------"
     &                 ,"------------"
       do iblk=1,nblk
         write(*,"(2x,i4,4x,a9,7(f14.8))") iblk,trim(gtype(iblk))
     &             ,cd_blk_ave(iblk),cly_blk_ave(iblk),clz_blk_ave(iblk)
     &             ,cl_blk_ave(iblk)/(cd_blk_ave(iblk)+1.d-60)
     &             ,ca_blk_ave(iblk),cny_blk_ave(iblk),cnz_blk_ave(iblk)
       enddo
       write(*,"(x,2a9,x,7(x,a12,x))") " --------","---------"
     &                 ,"------------","------------","------------"
     &                 ,"------------","------------","------------"
     &                 ,"------------"
       write(*,"(2x,a5,3x,a9,7(f14.8))")"Total",""
     &                        ,cd_tot,cly_tot,clz_tot
     &                        ,cl_tot/cd_tot
     &                        ,ca_tot,cny_tot,cnz_tot
       write(*,"(2x,a17,7(f14.8))")"Modified newtonian"
     &                        ,cd_tot*cp0,cly_tot*cp0,clz_tot*cp0
     &                        ,cl_tot/cd_tot
     &      ,ca_tot*cp0,cny_tot*cp0,cnz_tot*cp0
       write(*,"(2x,a17,7(f14.8))")"Standard deviation"
     &      ,sqrt(cd_var)*cp0,sqrt(cly_var)*cp0,sqrt(clz_var)*cp0,0
     &      ,sqrt(ca_var)*cp0,sqrt(cny_var)*cp0,sqrt(cnz_var)*cp0
       write(*,121)
     &            "=========================================="
     &           ,"=========================================="
     &           ,"=========================================="
       write(*,*) " Moment coefficients of each blocks             "
       write(*,122)
     &            "------------------------------------------"
     &           ,"------------------------------------------"
     &           ,"------------------------------------------"
       write(*,"(x,2a9,7(2x,a8,4x))") " Block # ","Shape "
     &                       ,"Cm,0,x ","Cm,0,y ","Cm,0,z "
     &                       ,"Cm,cg,x ","Cm,cg,y ","Cm,cg,z "
       write(*,"(x,2a9,x,7(x,a12,x))") " ------- ","----------"
     &                 ,"------------","------------","------------"
     &                 ,"------------","------------","------------"
       do iblk=1,nblk
         write(*,"(2x,i4,4x,a9,7(f14.8))") iblk,trim(gtype(iblk))
     &        ,cm0x_blk_ave(iblk),cm0y_blk_ave(iblk),cm0z_blk_ave(iblk)
     &        ,cmx_blk_ave(iblk),cmy_blk_ave(iblk),cmz_blk_ave(iblk)
       enddo
       write(*,"(x,2a9,x,7(x,a12,x))") " --------","---------"
     &                 ,"------------","------------","------------"
     &                 ,"------------","------------","------------"
       write(*,"(2x,a5,3x,a9,7(f14.8))")"Total",""
     &                        ,cm0x_tot,cm0y_tot,cm0z_tot
     &                        ,cmx_tot,cmy_tot,cmz_tot
       write(*,"(2x,a17,6(f14.8))")"Modified newtonian"
     &                     ,cm0x_tot*cp0,cm0y_tot*cp0,cm0z_tot*cp0
     &      ,cmx_tot*cp0,cmy_tot*cp0,cmz_tot*cp0
       write(*,"(2x,a17,6(f14.8))")"Standard deviation"
     &      ,sqrt(cm0x_var)*cp0,sqrt(cm0y_var)*cp0,sqrt(cm0z_var)*cp0
     &      ,sqrt(cmx_var)*cp0,sqrt(cmy_var)*cp0,sqrt(cmz_var)*cp0
       write(*,"(2x,a17,6(f14.8))")"Cm,cg-Cm,0        "
     &                     ,cmx_tot-cm0x_tot
     &                     ,cmy_tot-cm0y_tot
     &                     ,cmz_tot-cm0z_tot
       write(*,"(2x,a17,6(f14.8))")"Cm,cg-Cm,0 (MN)   "
     &                     ,(cmx_tot-cm0x_tot)*cp0
     &                     ,(cmy_tot-cm0y_tot)*cp0
     &                     ,(cmz_tot-cm0z_tot)*cp0
       write(*,"(2x,3a34)") 
     &            "=========================================="
     &           ,"=========================================="
     &           ,"=========================================="
       write(*,*) " Visible area                                       "
       write(*,"(2x,3a25)") 
     &            "-----------------------------------------"
     &           ,"-----------------------------------------"
     &           ,"-----------------------------------------"
       write(*,"(x,a9,5a14)") " Block # ","Visible","Projected"
     &                     ," ivisible"," nsample"
       write(*,"(x,a9,5a14)") " --------"
     &                     ,"----------","----------"
     &                     ,"----------","----------"
       do iblk=1,nblk
         write(*,"(2x,i4,4x,2(f14.8),i12,x,a1,i14)") iblk
     &             ,area_blk_ave(iblk),area_yz_blk_ave(iblk)
     &          ,ivsbl_blk_sum(iblk),"/",nsmpb(iblk)*ntry
       enddo
       write(*,"(x,a9,5a14)") " --------","","----------","----------"
     &                     ,"----------"
       write(*,"(2x,a5,3x,2(f14.8),i12,x,a1,i14)")"Total",area_tot
     &                                          ,area_yz_tot
     &           ,sum(ivsbl_blk_sum(1:nblk)),"/"
     &           ,sum(nsmpb(1:nblk))*ntry
       write(*,"(x,a9,3a14)") " --------","","----------",""
       write(*,"(2x,a5,17x,f14.8)")"Ref. ",area_ref
       write(*,"(2x,3a25)") 
     &            "=========================================="
     &           ,"=========================================="
     &           ,"=========================================="
      end subroutine mcnew_aero
!
      end module mcnew
!
!**********************************************************************
!  Main program
!**********************************************************************
      program mcnew_main
        use mcnew, only : mcnew_aero
        use geoms, only : nblk,set_geom_params,xmin,xmax,zmax
        implicit none
!
        real(8),parameter :: pi = dacos(-1.d0)
!
! input parameters for calculation contral
        integer, allocatable :: nsmpb(:)! Number of sample points/trial
        integer :: ntry      ! Number of trials
     &            ,iout      ! Interval for intermediate output
!
! input parameters for geometry
        real(8) :: cg(3)     ! xyz coordinates of center of gravity
!
! input parameters for freesteam 
        real(8) :: alp_deg   ! pitch angle in degree
     &            ,alp_rad   ! pitch angle in radian
     &            ,bet_deg   ! yaw angle in degree
     &            ,bet_rad   ! yaw angle in radian
        real(8) :: gam       ! specific heat ratio
     &            ,amach     ! Mach number
!
! Pressure correction factor in Modified Newtonian (MN)
        real(8) :: cp0       ! Maximum pressure * 0.5
        real(8) :: amach2,gamach2,gamp1,gamm1 ! Temporal storage
!
! Aerodynamic coefficients array
        real(8) :: aeroc(13)
                 ! aeroc(1): CD , Drag coef
                 ! aeroc(2): CLy, Lift coef in y-direction
                 ! aeroc(3): CLz, Lift coef in y-direction
                 ! aeroc(4): L/D, Lift-to-drag ratio
                 ! aeroc(5): CA , Axial force coef
                 ! aeroc(6): CNy, Normal force coef in y
                 ! aeroc(7): CNz, Normal force coef in z
                 ! aeroc(8): Cmx, Moment coef aroud x-axis passing c.g.
                 ! aeroc(9): Cmy, Moment coef aroud y-axis passing c.g.
                 ! aeroc(10):Cmz, Moment coef aroud z-axis passing c.g.
                 ! aeroc(11):Cm0x, Moment coef aroud x-axis
                 ! aeroc(12):Cm0y, Moment coef aroud y-axis
                 ! aeroc(13):Cm0z, Moment coef aroud z-axis
!
! variables for printing results
       character(50) :: fname_result1 = "NEW_aero_coefs.dat"
     &                 ,fname_result2 = "MN_aero_coefs.dat"
       character(50) :: 
     &   ftitle1 = "  Aerodynamic coefs by Newtoninan theory"
     &  ,ftitle2 = "  Aerodynamic coefs by Modified Newtoninan theory"
       character(7) :: vname(17) = (
     &                       / "Mach   ","Gamma  ","Alpha  ","Beta   "
     &                        ,"CD     ","CL,y   ","CL,z   ","L/D    "
     &                        ,"CA     ","CN,y   ","CN,z   "
     &                        ,"Cm,x cg","Cm,y cg","Cm,z cg"
     &                        ,"Cm,x O ","Cm,y O ","Cm,z O "  /)
       real(8)      :: vars(17)
!
! variables for system clock
       integer :: time_begin_c,time_end_c,time_real,CountPerSec
       real(8) :: time_begin_s,time_end_s,time_cpu
       integer :: nsmptot
!
! tempral variables for calculation
       integer :: i
!
!      ----------------------------------------------------------------
!      Read input file and set parameters to define geometries
!      ----------------------------------------------------------------
        call set_geom_params
!
!      ----------------------------------------------------------------
!       Read centor of gravity coordinate for moment calculation
!      ----------------------------------------------------------------
        read(*,*)
        read(*,*) (cg(i),i=1,3)
!
!      ----------------------------------------------------------------
!      Read sample number and trial number
!      ----------------------------------------------------------------
        allocate(nsmpb(nblk))
        read(*,*)
        do i=1,nblk
          read(*,*)nsmpb(i)
        enddo
!
!      ----------------------------------------------------------------
!      Read number of trials and output interval
!      ----------------------------------------------------------------
        read(*,*)
        read(*,*)ntry
        read(*,*)iout
        nsmptot = sum(nsmpb(1:nblk)) * ntry ! Total sample number
!
!      ----------------------------------------------------------------
!      Read angle-of-attack in deg and covert to rad
!      ----------------------------------------------------------------
        read(*,*)
        read(*,*)alp_deg
        read(*,*)bet_deg
        alp_rad = alp_deg / 180.d0 * pi
        bet_rad = bet_deg / 180.d0 * pi
!
!      ----------------------------------------------------------------
!      Read mach number and specific heat ratio
!      ----------------------------------------------------------------
        read(*,*)
        read(*,*)amach
        read(*,*)gam
!
!      ----------------------------------------------------------------
!      Calculate max pressire coef for modifeid newtonian theory
!      ----------------------------------------------------------------
        gamp1  = gam + 1.d0
        gamm1  = gam - 1.d0
        amach2 = amach**2
        gamach2 = gam*amach2
        cp0 = 1.0d0 /(gam*amach2)
     &      *((gamp1**2 *amach2/(4.d0*gamach2-2.d0*gamm1))**(gam/gamm1)
     &      * (-gamm1 + 2.d0*gamach2)/gamp1-1.d0)
!
!      ----------------------------------------------------------------
!      Print summary of condition
!      ----------------------------------------------------------------
        print *, ""
        print *, " ================================================"
        print *, "  Free stream conditions                         "
        print *, " ------------------------------------------------"
        print *, "  Mach number, M:", amach
        print *, "  Specific heat ratio:", gam
        print *, "  Maxumum pressure coef, Cp0:", cp0*2.d0
        print *, "  Pitch angle:", alp_deg, "[Deg]"
        print *, "  Yaw angle:", bet_deg, "[Deg]"
        print *, " ================================================"
        print *, "  Computation conditions                         "
        print *, " ------------------------------------------------"
        print *, "  Number of sample points:"
        do i = 1,nblk
        print *, "     Block #",i,":",nsmpb(i),"/trial " 
        enddo
        print *, "  Number of trial test: ",ntry
        print *, "  Interval for output:" ,iout
        print *, " ================================================"
!
!      ----------------------------------------------------------------
!      Open files for outputting results
!      ----------------------------------------------------------------
        open(80,file=fname_result1)
        open(81,file=fname_result2)
        call print_results(80,0,vname,vars,ftitle1)
        call print_results(81,0,vname,vars,ftitle2)
!
!      ----------------------------------------------------------------
!      Check system clock at the begining
!      ----------------------------------------------------------------
        call system_clock(time_begin_c,CountPerSec)
        call cpu_time(time_begin_s)
!
!      ----------------------------------------------------------------
!      Aerodynamic calculation
!      ----------------------------------------------------------------
        call mcnew_aero(nsmpb,ntry,iout  ! in : Sample numbers
     &                 ,alp_rad,bet_rad  ! in : Angle of attack
     &                 ,cg(1:3)          ! in : Center of gravity
     &                 ,cp0              ! in : Modification factor
     &                 ,aeroc(1:13))     ! out: Aerodynamic coeffs
!
!      ----------------------------------------------------------------
!      Output results
!      ----------------------------------------------------------------
        vars(1)    = amach
        vars(2)    = gam
        vars(3)    = alp_deg
        vars(4)    = bet_deg
        vars(5:17) = aeroc(1:13)
!
        call print_results(80,1,vname,vars)
!
        vars(5:17) = aeroc(1:13)*cp0
        vars(8)    = aeroc(4)
        call print_results(81,1,vname(1:17),vars(1:17))
!
        close(80)
        close(81)
!
!      ----------------------------------------------------------------
!      Check system clock at the end
!      ----------------------------------------------------------------
       call system_clock(time_end_c,CountPerSec)
       call cpu_time(time_end_s)
!
!      ----------------------------------------------------------------
!      Calculate computation time
!      ----------------------------------------------------------------
       time_real = real(time_end_c-time_begin_c)/CountPerSec
       time_cpu = real(time_end_s-time_begin_s)
!
!      ----------------------------------------------------------------
!      Disply computaion info
!      ----------------------------------------------------------------
       print *, ""
       print *, " ================================================"
       print *, "  Computation info                               "
       print *, " ------------------------------------------------"
       print *, "  CPU time :", time_cpu,"sec"
       print *, "  Real time:", time_real,"sec"
       print *, "  Number of trials", ntry,"times"
       print *, "  Total number of sample points", nsmptot,"points"
       print *, " ================================================"
!
!      End of main
!      ================================================================
       contains
!      ----------------------------------------------------------------
!      Subroutine to print results on file
!      ----------------------------------------------------------------
         subroutine print_results(iunit,id,vname,vars,title)
!
         implicit none
         character(7),intent(in)  :: vname(:)
         real(8),intent(in)       :: vars(:)
         integer,intent(in)       :: iunit,id
         character(50),intent(in),optional   :: title
         integer                  :: nvar
         integer*2 tabc / 2313 /
!
          nvar = size(vars)
          if(id == 0) then
            if(present(title)) then
              write(iunit,123) "========================="
     &                         ,"========================="
     &                         ,"========================="
              write(iunit,124) title
              write(iunit,123) "========================="
     &                         ,"========================="
     &                         ,"========================="
            endif
!           write(iunit,125)((vname(i),tabc),i=1,nvar)
            write(iunit,127)(vname(i),i=1,nvar)
        elseif (id==1) then
!           write(iunit,126)((vars(i),tabc),i=1,nvar)
            write(iunit,128)(vars(i),i=1,nvar)
        endif
!
 123    format('#',3a25)
 124    format('#',a50)
!125    format('#',4(x,a7,a1),20(x,a7,4x,a1))
!126    format(" ",4(f8.2,a1),20(f12.8, a1))
 127    format('#',4(x,a7),20(x,a7,5x))
 128    format(" ",4(x,f8.2),20(x,f12.8))
!
        end subroutine print_results
!
      end program mcnew_main
! 
!_______________________________________________________________________
! END OF PROGRAM 
! ************* ************* ************* ************* *************
! **   **    ** ****      *** *  *       ** **        *** ***  ***  ***
! *  **   **  * **   ******** *   *****   * *   ***   *** **   ***   **
! *   ** **   * *   ******* * *   *****   * *       *** * *   ** **   *
! **   ***   ** **  *****  ** *  ******  ** **   *****  * *  **   **  *
! ***  ***  *** ***       *** * ******* *** ***        ** **    **   **
! ************* ************* ************* ************* *************
