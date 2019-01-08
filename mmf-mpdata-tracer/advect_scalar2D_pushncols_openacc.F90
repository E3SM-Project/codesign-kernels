
!Positively definite monotonic advection with non-oscillatory option
program test_advect
  implicit none

  !User changeable
  integer, parameter :: nslices = 48
  integer, parameter :: nz      = 58
  integer, parameter :: nx      = 32
  
  !From grid.F90
 !integer, parameter :: rp      = selected_real_kind(7)    !For single precision, use this line
  integer, parameter :: rp      = selected_real_kind(13)   !For double precision, use this line
  integer, parameter :: nzm     = nz-1
  integer, parameter :: nxp3    = nx+3
  integer, parameter :: nxp2    = nx+2
  integer, parameter :: nxp1    = nx+1
  integer, parameter :: dimx1_u = -1
  integer, parameter :: dimx2_u = nxp3
  integer, parameter :: dimy1_u = 1
  integer, parameter :: dimy2_u = 1
  integer, parameter :: dimx1_w = -1
  integer, parameter :: dimx2_w = nxp2
  integer, parameter :: dimy1_w = 1
  integer, parameter :: dimy2_w = 1
  integer, parameter :: dimx1_s = -2
  integer, parameter :: dimx2_s = nxp3
  integer, parameter :: dimy1_s = 1
  integer, parameter :: dimy2_s = 1
  real(rp) :: adz(nslices,nzm)   ! ratio of the thickness of scalar levels to dz 

  !Globals
  real(rp) :: f    (nslices,dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
  real(rp) :: u    (nslices,dimx1_u:dimx2_u, dimy1_u:dimy2_u, nzm)
  real(rp) :: w    (nslices,dimx1_w:dimx2_w, dimy1_w:dimy2_w, nz )
  real(rp) :: rho  (nslices,nzm)
  real(rp) :: rhow (nslices,nz )
  real(rp) :: flux (nslices,nz )

  !Save for diffs
  real(rp) :: f_save    (nslices,dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
  real(rp) :: u_save    (nslices,dimx1_u:dimx2_u, dimy1_u:dimy2_u, nzm)
  real(rp) :: w_save    (nslices,dimx1_w:dimx2_w, dimy1_w:dimy2_w, nz )
  real(rp) :: rho_save  (nslices,nzm)
  real(rp) :: rhow_save (nslices,nz )
  real(rp) :: flux_save (nslices,nz )

  call init()
  call advect_scalar2D_cpu(f,u,w,rho,rhow,flux)
  call save()
!
!  call init()
!  call advect_scalar2D_openacc_1(f,u,w,rho,rhow,flux)
!  call compare()
!
!  call init()
!  call advect_scalar2D_openacc_2(f,u,w,rho,rhow,flux)
!  call compare()
!
!  call init()
!  call advect_scalar2D_openacc_1(f,u,w,rho,rhow,flux)
!  call compare()
!
!  call init()
!  call advect_scalar2D_openacc_2(f,u,w,rho,rhow,flux)
!  call compare()


contains


  subroutine advect_scalar2D_openacc_2(f, u, w, rho, rhow, flux)
    implicit none
    real(rp), intent(inout) :: f    (nslices,dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
    real(rp), intent(in   ) :: u    (nslices,dimx1_u:dimx2_u, dimy1_u:dimy2_u, nzm)
    real(rp), intent(in   ) :: w    (nslices,dimx1_w:dimx2_w, dimy1_w:dimy2_w, nz )
    real(rp), intent(in   ) :: rho  (nslices,nzm)
    real(rp), intent(in   ) :: rhow (nslices,nz )
    real(rp), intent(  out) :: flux (nslices,nz )
    real(rp) :: mx   (nslices,0:nxp1,1,nzm)
    real(rp) :: mn   (nslices,0:nxp1,1,nzm)
    real(rp) :: uuu  (nslices,-1:nxp3,1,nzm)
    real(rp) :: www  (nslices,-1:nxp2,1,nz)
    real(rp) :: iadz (nslices,nzm)
    real(rp) :: irho (nslices,nzm)
    real(rp) :: irhow(nslices,nzm)
    real(rp) :: eps, dd
    integer  :: i,j,k,ic,ib,kc,kb, sl
    logical  :: nonos
    real(rp) :: x1, x2, a, b, a1, a2, y
    real(rp) :: andiff,across,pp,pn
    integer(8) :: t1, t2, tr
    
    !Statement functions
    andiff(x1,x2,a,b) = (abs(a)-a*a*b)*0.5*(x2-x1)
    across(x1,a1,a2)  = 0.03125*a1*a2*x1
    pp(y)             =  max(0._rp,y)
    pn(y)             = -min(0._rp,y)
    
    !Initialization
    nonos = .true.
    eps = 1.e-10
    j=1

    !$acc enter data pcreate(f,u,w,rho,rhow,flux,mx,mn,uuu,www,iadz,irho,irhow,rhow,adz)

    !$acc update device(f,u,w,rho,adz,rhow) async(1)

    !$acc wait
    call system_clock(t1)
    
    !$acc parallel loop gang vector collapse(3) present(mn,mx,f,www,uuu,u,w) private(kb,kc,ib,ic) async(1)
    do k=1,nzm
      do i=-1,nxp3
        do sl = 1 , nslices
          kc=min(nzm,k+1)
          kb=max(1,k-1)
          ib=i-1
          ic=i+1
          if(nonos) then
            if (i >= 0 .and. i <= nxp1) then
              mx(sl,i,j,k)=max(f(sl,ib,j,k),f(sl,ic,j,k),f(sl,i,j,kb),f(sl,i,j,kc),f(sl,i,j,k))
              mn(sl,i,j,k)=min(f(sl,ib,j,k),f(sl,ic,j,k),f(sl,i,j,kb),f(sl,i,j,kc),f(sl,i,j,k))
            endif
          endif  ! nonos
          if (k == 1) www(sl,i,j,nz) = 0.
          if (i >= -1 .and. i <= nxp3) then
            uuu(sl,i,j,k)=max(0._rp,u(sl,i,j,k))*f(sl,i-1,j,k )+min(0._rp,u(sl,i,j,k))*f(sl,i,j,k)
          endif
          if (i >= -1 .and. i <= nxp2) then
            www(sl,i,j,k)=max(0._rp,w(sl,i,j,k))*f(sl,i  ,j,kb)+min(0._rp,w(sl,i,j,k))*f(sl,i,j,k)
          endif
        enddo
      enddo
    enddo

    !$acc parallel loop gang vector collapse(2) present(flux,www,rho,adz,irho,iadz,irhow,rhow) async(1)
    do k=1,nzm
      do sl = 1 , nslices
        flux(sl,k) = 0.
        do i=1,nx
          flux(sl,k) = flux(sl,k) + www(sl,i,j,k)  
        enddo
        irho(sl,k) = 1./rho(sl,k)
        iadz(sl,k) = 1./adz(sl,k)
        irhow(sl,k)=1./(rhow(sl,k)*adz(sl,k))
      enddo
    enddo

    !$acc parallel loop gang vector collapse(3) present(f,uuu,www,iadz,irho) async(1)
    do k=1,nzm
      do i=-1,nxp2
        do sl = 1 , nslices
          f(sl,i,j,k) = f(sl,i,j,k) - (uuu(sl,i+1,j,k)-uuu(sl,i,j,k) + (www(sl,i,j,k+1)-www(sl,i,j,k))*iadz(sl,k))*irho(sl,k)            
        enddo
      enddo
    enddo 
    
    !$acc parallel loop gang vector collapse(3) present(uuu,f,u,irho,w,www,irhow,mn,mx,adz) private(kc,kb,dd,ib,ic) async(1)
    do k=1,nzm
      do i=0,nxp2
        do sl = 1 , nslices
          kc = min(nzm,k+1)
          kb = max(1,k-1)
          ib = i-1
          ic = i+1
          if (i >= 0 .and. i <= nxp2) then
            dd = 2./(kc-kb)/adz(sl,k)
            uuu(sl,i,j,k)=andiff(f(sl,ib,j,k ),f(sl,i,j,k),u(sl,i,j,k),irho(sl,k)) - &
                          across( dd*(f(sl,ib,j,kc)+f(sl,i,j,kc)-f(sl,ib,j,kb)-f(sl,i,j,kb))   , &
                                  u(sl,i,j,k), w(sl,ib,j,k)+w(sl,ib,j,kc)+w(sl,i,j,k)+w(sl,i,j,kc) ) * irho(sl,k)
          endif
          if (i >= 0 .and. i <= nxp1) then
            www(sl,i,j,k)=andiff(f(sl,i,j,kb),f(sl,i,j,k),w(sl,i,j,k),irhow(sl,k)) - &
                          across( f(sl,ic,j,kb)+f(sl,ic,j,k)-f(sl,ib,j,kb)-f(sl,ib,j,k)        , &
                                  w(sl,i,j,k), u(sl,i,j,kb)+u(sl,i,j,k)+u(sl,ic,j,k)+u(sl,ic,j,kb) ) * irho(sl,k)
          endif
          if (k == 1) www(sl,i,j,1)=0.
          if (i >= 0 .and. i <= nxp1) then
            if(nonos) then
              mx(sl,i,j,k)=max(f(sl,ib,j,k),f(sl,ic,j,k),f(sl,i,j,kb),f(sl,i,j,kc),f(sl,i,j,k),mx(sl,i,j,k))
              mn(sl,i,j,k)=min(f(sl,ib,j,k),f(sl,ic,j,k),f(sl,i,j,kb),f(sl,i,j,kc),f(sl,i,j,k),mn(sl,i,j,k))
            endif
          endif
        enddo
      enddo
    enddo

    if(nonos) then
      !$acc parallel loop gang vector collapse(3) present(mx,mn,f,rho,uuu,www,iadz) private(kc,kb,ib,ic) async(1)
      do k=1,nzm
        do i=0,nxp1
          do sl = 1 , nslices
            kc=min(nzm,k+1)
            kb=max(1  ,k-1)
            ic=i+1
            mx(sl,i,j,k)=rho(sl,k)*(mx(sl,i,j,k)- f(sl,i,j,k))/(pn(uuu(sl,ic,j,k)) + pp(uuu(sl,i,j,k))+&
                         iadz(sl,k)*(pn(www(sl,i,j,kc)) + pp(www(sl,i,j,k)))+eps)  
            mn(sl,i,j,k)=rho(sl,k)*( f(sl,i,j,k)-mn(sl,i,j,k))/(pp(uuu(sl,ic,j,k)) + pn(uuu(sl,i,j,k))+&
                         iadz(sl,k)*(pp(www(sl,i,j,kc)) + pn(www(sl,i,j,k)))+eps)  
          enddo
        enddo
      enddo
      !$acc parallel loop gang vector collapse(3) present(uuu,mx,mn,www) private(kb,ib) async(1)
      do k=1,nzm
        do i=1,nxp1
          do sl = 1 , nslices
            kb=max(1,k-1)
            ib=i-1
            uuu(sl,i,j,k) = pp(uuu(sl,i,j,k))*min(1._rp,mx(sl,i,j,k), mn(sl,ib,j,k)) - pn(uuu(sl,i,j,k))*min(1._rp,mx(sl,ib,j,k),mn(sl,i,j,k))
            if (i >= 0 .and. i <= nx) then
              www(sl,i,j,k) = pp(www(sl,i,j,k))*min(1._rp,mx(sl,i,j,k), mn(sl,i,j,kb)) - pn(www(sl,i,j,k))*min(1._rp,mx(sl,i,j,kb),mn(sl,i,j,k))
            endif
          enddo
        enddo
      enddo
      !$acc parallel loop gang vector collapse(2) present(flux,www) async(1)
      do k=1,nzm
        do sl = 1 , nslices
          do i=1,nx
            flux(sl,k) = flux(sl,k) + www(sl,i,j,k)  
          enddo
        enddo
      enddo
    endif ! nonos
    
    !$acc parallel loop gang collapse(3) present(f,uuu,www,iadz,irho) private(kc) async(1)
    do k=1,nzm
      do i=1,nx
        do sl = 1 , nslices
          kc=k+1
          f(sl,i,j,k)= max( 0._rp , f(sl,i,j,k) - (uuu(sl,i+1,j,k)-uuu(sl,i,j,k) + (www(sl,i,j,k+1)-www(sl,i,j,k))*iadz(sl,k))*irho(sl,k) )
        enddo
      enddo
    enddo 

    !$acc wait(1)
    call system_clock(t2,tr)
    write(*,*) 'OpenACC-2 Timing: ', dble(t2-t1)/dble(tr)

    !$acc update host(flux,f) async(1)
    !$acc wait(1)
  
  end subroutine advect_scalar2D_openacc_2


  subroutine advect_scalar2D_openacc_1(f, u, w, rho, rhow, flux)
    implicit none
    real(rp), intent(inout) :: f    (nslices,dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
    real(rp), intent(in   ) :: u    (nslices,dimx1_u:dimx2_u, dimy1_u:dimy2_u, nzm)
    real(rp), intent(in   ) :: w    (nslices,dimx1_w:dimx2_w, dimy1_w:dimy2_w, nz )
    real(rp), intent(in   ) :: rho  (nslices,nzm)
    real(rp), intent(in   ) :: rhow (nslices,nz )
    real(rp), intent(  out) :: flux (nslices,nz )
    real(rp) :: mx   (nslices,0:nxp1,1,nzm)
    real(rp) :: mn   (nslices,0:nxp1,1,nzm)
    real(rp) :: uuu  (nslices,-1:nxp3,1,nzm)
    real(rp) :: www  (nslices,-1:nxp2,1,nz)
    real(rp) :: iadz (nslices,nzm)
    real(rp) :: irho (nslices,nzm)
    real(rp) :: irhow(nslices,nzm)
    real(rp) :: eps, dd
    integer  :: i,j,k,ic,ib,kc,kb, sl
    logical  :: nonos
    real(rp) :: x1, x2, a, b, a1, a2, y
    real(rp) :: andiff,across,pp,pn
    integer(8) :: t1, t2, tr
    
    !Statement functions
    andiff(x1,x2,a,b) = (abs(a)-a*a*b)*0.5*(x2-x1)
    across(x1,a1,a2)  = 0.03125*a1*a2*x1
    pp(y)             =  max(0._rp,y)
    pn(y)             = -min(0._rp,y)
    
    !Initialization
    nonos = .true.
    eps = 1.e-10
    j=1

    !$acc enter data pcreate(f,u,w,rho,rhow,flux,mx,mn,uuu,www,iadz,irho,irhow,rhow,adz)

    !$acc update device(f,u,w,rho,adz,rhow) async(1)

    !$acc wait
    call system_clock(t1)

    !$acc parallel loop gang vector collapse(2) present(www) async(1)
    do i = dimx1_w , dimx2_w
      do sl = 1 , nslices
        www(sl,i,j,nz)=0.
      enddo
    enddo
    
    if(nonos) then
      !$acc parallel loop gang vector collapse(3) present(mn,mx,f) private(kb,kc,ib,ic) async(1)
      do k=1,nzm
        do i=0,nxp1
          do sl = 1 , nslices
            kc=min(nzm,k+1)
            kb=max(1,k-1)
            ib=i-1
            ic=i+1
            mx(sl,i,j,k)=max(f(sl,ib,j,k),f(sl,ic,j,k),f(sl,i,j,kb),f(sl,i,j,kc),f(sl,i,j,k))
            mn(sl,i,j,k)=min(f(sl,ib,j,k),f(sl,ic,j,k),f(sl,i,j,kb),f(sl,i,j,kc),f(sl,i,j,k))
          enddo
        enddo
      enddo
    endif  ! nonos
    
    !$acc parallel loop gang vector collapse(3) present(uuu,u,f) private(kb) async(1)
    do k=1,nzm
      do i=-1,nxp3
        do sl = 1 , nslices
          kb=max(1,k-1)
          uuu(sl,i,j,k)=max(0._rp,u(sl,i,j,k))*f(sl,i-1,j,k )+min(0._rp,u(sl,i,j,k))*f(sl,i,j,k)
        enddo
      enddo
    enddo
    !$acc parallel loop gang vector collapse(3) present(www,w,f) private(kb) async(1)
    do k=1,nzm
      do i=-1,nxp2
        do sl = 1 , nslices
          kb=max(1,k-1)
          www(sl,i,j,k)=max(0._rp,w(sl,i,j,k))*f(sl,i  ,j,kb)+min(0._rp,w(sl,i,j,k))*f(sl,i,j,k)
        enddo
      enddo
    enddo
    !$acc parallel loop gang vector collapse(2) present(flux,www) async(1)
    do k=1,nzm
      do sl = 1 , nslices
        flux(sl,k) = 0.
        do i=1,nx
          flux(sl,k) = flux(sl,k) + www(sl,i,j,k)  
        enddo
      enddo
    enddo

    !$acc parallel loop gang vector collapse(2) present(rho,adz,irho,iadz) async(1)
    do k=1,nzm
      do sl = 1 , nslices
        irho(sl,k) = 1./rho(sl,k)
        iadz(sl,k) = 1./adz(sl,k)
      enddo
    enddo

    !$acc parallel loop gang vector collapse(3) present(f,uuu,www,iadz,irho) async(1)
    do k=1,nzm
      do i=-1,nxp2
        do sl = 1 , nslices
          f(sl,i,j,k) = f(sl,i,j,k) - (uuu(sl,i+1,j,k)-uuu(sl,i,j,k) + (www(sl,i,j,k+1)-www(sl,i,j,k))*iadz(sl,k))*irho(sl,k)            
        enddo
      enddo
    enddo 
    
    !$acc parallel loop gang vector collapse(2) present(irhow,rhow,adz) async(1)
    do k=1,nzm
      do sl = 1 , nslices
        irhow(sl,k)=1./(rhow(sl,k)*adz(sl,k))
      enddo
    enddo

    !$acc parallel loop gang vector collapse(3) present(uuu,f,u,irho,w,adz) private(kc,kb,dd,ib) async(1)
    do k=1,nzm
      do i=0,nxp2
        do sl = 1 , nslices
          kc = min(nzm,k+1)
          kb = max(1,k-1)
          dd = 2./(kc-kb)/adz(sl,k)
          ib = i-1
          uuu(sl,i,j,k)=andiff(f(sl,ib,j,k ),f(sl,i,j,k),u(sl,i,j,k),irho(sl,k)) - &
                        across( dd*(f(sl,ib,j,kc)+f(sl,i,j,kc)-f(sl,ib,j,kb)-f(sl,i,j,kb))   , &
                                u(sl,i,j,k), w(sl,ib,j,k)+w(sl,ib,j,kc)+w(sl,i,j,k)+w(sl,i,j,kc) ) * irho(sl,k)
        enddo
      enddo
    enddo

    !$acc parallel loop gang vector collapse(3) present(www,f,w,irhow,irho,u) private(kc,kb,ic,ib) async(1)
    do k=1,nzm
      do i=0,nxp1
        do sl = 1 , nslices
          kc = min(nzm,k+1)
          kb = max(1,k-1)
          ib=i-1
          ic=i+1
          www(sl,i,j,k)=andiff(f(sl,i,j,kb),f(sl,i,j,k),w(sl,i,j,k),irhow(sl,k)) - &
                        across( f(sl,ic,j,kb)+f(sl,ic,j,k)-f(sl,ib,j,kb)-f(sl,ib,j,k)        , &
                                w(sl,i,j,k), u(sl,i,j,kb)+u(sl,i,j,k)+u(sl,ic,j,k)+u(sl,ic,j,kb) ) * irho(sl,k)
        enddo
      enddo
    enddo

    !$acc parallel loop gang vector collapse(2) present(www) async(1)
    do i = dimx1_w , dimx2_w
      do sl = 1 , nslices
        www(sl,i,j,1)=0.
      enddo
    enddo
  
    if(nonos) then
      !$acc parallel loop gang vector collapse(3) present(mx,mn,f) private(kc,kb,ib,ic) async(1)
      do k=1,nzm
        do i=0,nxp1
          do sl = 1 , nslices
            kc=min(nzm,k+1)
            kb=max(1  ,k-1)
            ib=i-1
            ic=i+1
            mx(sl,i,j,k)=max(f(sl,ib,j,k),f(sl,ic,j,k),f(sl,i,j,kb),f(sl,i,j,kc),f(sl,i,j,k),mx(sl,i,j,k))
            mn(sl,i,j,k)=min(f(sl,ib,j,k),f(sl,ic,j,k),f(sl,i,j,kb),f(sl,i,j,kc),f(sl,i,j,k),mn(sl,i,j,k))
          enddo
        enddo
      enddo
      !$acc parallel loop gang vector collapse(3) present(mx,mn,f,rho,uuu,www,iadz) private(kc,kb,ib,ic) async(1)
      do k=1,nzm
        do i=0,nxp1
          do sl = 1 , nslices
            kc=min(nzm,k+1)
            kb=max(1  ,k-1)
            ic=i+1
            mx(sl,i,j,k)=rho(sl,k)*(mx(sl,i,j,k)- f(sl,i,j,k))/(pn(uuu(sl,ic,j,k)) + pp(uuu(sl,i,j,k))+&
                         iadz(sl,k)*(pn(www(sl,i,j,kc)) + pp(www(sl,i,j,k)))+eps)  
            mn(sl,i,j,k)=rho(sl,k)*( f(sl,i,j,k)-mn(sl,i,j,k))/(pp(uuu(sl,ic,j,k)) + pn(uuu(sl,i,j,k))+&
                         iadz(sl,k)*(pp(www(sl,i,j,kc)) + pn(www(sl,i,j,k)))+eps)  
          enddo
        enddo
      enddo
      !$acc parallel loop gang vector collapse(3) present(uuu,mx,mn) private(kb,ib) async(1)
      do k=1,nzm
        do i=1,nxp1
          do sl = 1 , nslices
            kb=max(1,k-1)
            ib=i-1
            uuu(sl,i,j,k) = pp(uuu(sl,i,j,k))*min(1._rp,mx(sl,i,j,k), mn(sl,ib,j,k)) - pn(uuu(sl,i,j,k))*min(1._rp,mx(sl,ib,j,k),mn(sl,i,j,k))
          enddo
        enddo
      enddo
      !$acc parallel loop gang vector collapse(3) present(www,mx,mn) private(kb,ib) async(1)
      do k=1,nzm
        do i=1,nx
          do sl = 1 , nslices
            kb=max(1,k-1)
            ib=i-1
            www(sl,i,j,k) = pp(www(sl,i,j,k))*min(1._rp,mx(sl,i,j,k), mn(sl,i,j,kb)) - pn(www(sl,i,j,k))*min(1._rp,mx(sl,i,j,kb),mn(sl,i,j,k))
          enddo
        enddo
      enddo
      !$acc parallel loop gang vector collapse(2) present(flux,www) async(1)
      do k=1,nzm
        do sl = 1 , nslices
          do i=1,nx
            flux(sl,k) = flux(sl,k) + www(sl,i,j,k)  
          enddo
        enddo
      enddo
    endif ! nonos
    
    !$acc parallel loop gang collapse(3) present(f,uuu,www,iadz,irho) private(kc) async(1)
    do k=1,nzm
      do i=1,nx
        do sl = 1 , nslices
          kc=k+1
          f(sl,i,j,k)= max( 0._rp , f(sl,i,j,k) - (uuu(sl,i+1,j,k)-uuu(sl,i,j,k) + (www(sl,i,j,k+1)-www(sl,i,j,k))*iadz(sl,k))*irho(sl,k) )
        enddo
      enddo
    enddo 

    !$acc wait(1)
    call system_clock(t2,tr)
    write(*,*) 'OpenACC-1 Timing: ', dble(t2-t1)/dble(tr)

    !$acc update host(flux,f) async(1)
    !$acc wait(1)
  
  end subroutine advect_scalar2D_openacc_1


  subroutine advect_scalar2D_cpu(f, u, w, rho, rhow, flux)
    implicit none
    real(rp), intent(inout) :: f    (nslices,dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
    real(rp), intent(in   ) :: u    (nslices,dimx1_u:dimx2_u, dimy1_u:dimy2_u, nzm)
    real(rp), intent(in   ) :: w    (nslices,dimx1_w:dimx2_w, dimy1_w:dimy2_w, nz )
    real(rp), intent(in   ) :: rho  (nslices,nzm)
    real(rp), intent(in   ) :: rhow (nslices,nz )
    real(rp), intent(  out) :: flux (nslices,nz )
    real(rp) :: mx   (nslices,0:nxp1,1,nzm)
    real(rp) :: mn   (nslices,0:nxp1,1,nzm)
    real(rp) :: uuu  (nslices,-1:nxp3,1,nzm)
    real(rp) :: www  (nslices,-1:nxp2,1,nz)
    real(rp) :: iadz (nslices,nzm)
    real(rp) :: irho (nslices,nzm)
    real(rp) :: irhow(nslices,nzm)
    real(rp) :: eps, dd
    integer  :: i,j,k,ic,ib,kc,kb, sl
    logical  :: nonos
    real(rp) :: x1, x2, a, b, a1, a2, y
    real(rp) :: andiff,across,pp,pn
    integer(8) :: t1, t2, tr
    
    !Statement functions
    andiff(x1,x2,a,b) = (abs(a)-a*a*b)*0.5*(x2-x1)
    across(x1,a1,a2)  = 0.03125*a1*a2*x1
    pp(y)             =  max(0._rp,y)
    pn(y)             = -min(0._rp,y)

    call system_clock(t1)
    
    !Initialization
    nonos = .true.
    eps = 1.e-10
    j=1
    www(:,:,:,nz)=0.
    
    if(nonos) then
      do k=1,nzm
        kc=min(nzm,k+1)
        kb=max(1,k-1)
        do i=0,nxp1
          do sl = 1 , nslices
            ib=i-1
            ic=i+1
            mx(sl,i,j,k)=max(f(sl,ib,j,k),f(sl,ic,j,k),f(sl,i,j,kb),f(sl,i,j,kc),f(sl,i,j,k))
            mn(sl,i,j,k)=min(f(sl,ib,j,k),f(sl,ic,j,k),f(sl,i,j,kb),f(sl,i,j,kc),f(sl,i,j,k))
          enddo
        enddo
      enddo
    endif  ! nonos
    
    do k=1,nzm
      kb=max(1,k-1)
      do i=-1,nxp3
        do sl = 1 , nslices
          uuu(sl,i,j,k)=max(0._rp,u(sl,i,j,k))*f(sl,i-1,j,k )+min(0._rp,u(sl,i,j,k))*f(sl,i,j,k)
        enddo
      enddo
      do i=-1,nxp2
        do sl = 1 , nslices
          www(sl,i,j,k)=max(0._rp,w(sl,i,j,k))*f(sl,i  ,j,kb)+min(0._rp,w(sl,i,j,k))*f(sl,i,j,k)
        enddo
      enddo
      do sl = 1 , nslices
        flux(sl,k) = 0.
      enddo
      do i=1,nx
        do sl = 1 , nslices
          flux(sl,k) = flux(sl,k) + www(sl,i,j,k)  
        enddo
      enddo
    enddo
    
    do k=1,nzm
      do sl = 1 , nslices
        irho(sl,k) = 1./rho(sl,k)
        iadz(sl,k) = 1./adz(sl,k)
      enddo
      do i=-1,nxp2
        do sl = 1 , nslices
          f(sl,i,j,k) = f(sl,i,j,k) - (uuu(sl,i+1,j,k)-uuu(sl,i,j,k) + (www(sl,i,j,k+1)-www(sl,i,j,k))*iadz(sl,k))*irho(sl,k)            
        enddo
      enddo
    enddo 
    do k=1,nzm
      kc = min(nzm,k+1)
      kb = max(1,k-1)
      do sl = 1 , nslices
        irhow(sl,k)=1./(rhow(sl,k)*adz(sl,k))
      enddo
      do i=0,nxp2
        do sl = 1 , nslices
          dd = 2./(kc-kb)/adz(sl,k)
          ib = i-1
          uuu(sl,i,j,k)=andiff(f(sl,ib,j,k ),f(sl,i,j,k),u(sl,i,j,k),irho(sl,k)) - &
                        across( dd*(f(sl,ib,j,kc)+f(sl,i,j,kc)-f(sl,ib,j,kb)-f(sl,i,j,kb))   , &
                                u(sl,i,j,k), w(sl,ib,j,k)+w(sl,ib,j,kc)+w(sl,i,j,k)+w(sl,i,j,kc) ) * irho(sl,k)
        enddo
      enddo
      do i=0,nxp1
        do sl = 1 , nslices
          ib=i-1
          ic=i+1
          www(sl,i,j,k)=andiff(f(sl,i,j,kb),f(sl,i,j,k),w(sl,i,j,k),irhow(sl,k)) - &
                        across( f(sl,ic,j,kb)+f(sl,ic,j,k)-f(sl,ib,j,kb)-f(sl,ib,j,k)        , &
                                w(sl,i,j,k), u(sl,i,j,kb)+u(sl,i,j,k)+u(sl,ic,j,k)+u(sl,ic,j,kb) ) * irho(sl,k)
        enddo
      enddo
    enddo
    www(:,:,:,1) = 0.
  
    if(nonos) then
      do k=1,nzm
        kc=min(nzm,k+1)
        kb=max(1  ,k-1)
        do i=0,nxp1
          do sl = 1 , nslices
            ib=i-1
            ic=i+1
            mx(sl,i,j,k)=max(f(sl,ib,j,k),f(sl,ic,j,k),f(sl,i,j,kb),f(sl,i,j,kc),f(sl,i,j,k),mx(sl,i,j,k))
            mn(sl,i,j,k)=min(f(sl,ib,j,k),f(sl,ic,j,k),f(sl,i,j,kb),f(sl,i,j,kc),f(sl,i,j,k),mn(sl,i,j,k))
          enddo
        enddo
      enddo
      do k=1,nzm
        kc=min(nzm,k+1)
        do i=0,nxp1
          do sl = 1 , nslices
            ic=i+1
            mx(sl,i,j,k)=rho(sl,k)*(mx(sl,i,j,k)- f(sl,i,j,k))/(pn(uuu(sl,ic,j,k)) + pp(uuu(sl,i,j,k))+&
                      iadz(sl,k)*(pn(www(sl,i,j,kc)) + pp(www(sl,i,j,k)))+eps)  
            mn(sl,i,j,k)=rho(sl,k)*( f(sl,i,j,k)-mn(sl,i,j,k))/(pp(uuu(sl,ic,j,k)) + pn(uuu(sl,i,j,k))+&
                      iadz(sl,k)*(pp(www(sl,i,j,kc)) + pn(www(sl,i,j,k)))+eps)  
          enddo
        enddo
      enddo
      do k=1,nzm
        kb=max(1,k-1)
        do i=1,nxp1
          do sl = 1 , nslices
            ib=i-1
            uuu(sl,i,j,k) = pp(uuu(sl,i,j,k))*min(1._rp,mx(sl,i,j,k), mn(sl,ib,j,k)) - pn(uuu(sl,i,j,k))*min(1._rp,mx(sl,ib,j,k),mn(sl,i,j,k))
          enddo
        enddo
        do i=1,nx
          do sl = 1 , nslices
            www(sl,i,j,k) = pp(www(sl,i,j,k))*min(1._rp,mx(sl,i,j,k), mn(sl,i,j,kb)) - pn(www(sl,i,j,k))*min(1._rp,mx(sl,i,j,kb),mn(sl,i,j,k))
            flux(sl,k) = flux(sl,k) + www(sl,i,j,k)  
          enddo
        enddo
      enddo
    endif ! nonos
    
    do k=1,nzm
      kc=k+1
      do i=1,nx
        do sl = 1 , nslices
          f(sl,i,j,k)= max( 0._rp , f(sl,i,j,k) - (uuu(sl,i+1,j,k)-uuu(sl,i,j,k) + (www(sl,i,j,k+1)-www(sl,i,j,k))*iadz(sl,k))*irho(sl,k) )
        enddo
      enddo
    enddo 

    call system_clock(t2,tr)
    write(*,*) 'CPU Timing: ', dble(t2-t1)/dble(tr)
  
  end subroutine advect_scalar2D_cpu


  subroutine init()
    implicit none
    integer :: seed_size
    integer, allocatable :: seed(:)
    call random_seed(size = seed_size)
    allocate(seed(seed_size))
    seed = 100
    call random_seed(put=seed)

    call random_number(adz )
    call random_number(f   )
    call random_number(u   )
    call random_number(w   )
    call random_number(rho )
    call random_number(rhow)
    call random_number(flux)

    !$acc enter data pcreate(f,u,w,rho,rhow,flux,adz)
    !$acc update device(f,u,w,rho,rhow,flux,adz)

  end subroutine init


  subroutine save()
    implicit none
    f_save    = f   
    u_save    = u   
    w_save    = w   
    rho_save  = rho 
    rhow_save = rhow
    flux_save = flux
  end subroutine save


  subroutine compare()
    implicit none
    write(*,*) 'Relative L1 Error - f    : ' , sum(abs( f    - f_save    )) / sum(abs( f_save    )) 
    write(*,*) 'Relative L1 Error - flux : ' , sum(abs( flux - flux_save )) / sum(abs( flux_save )) 
  end subroutine compare


end program test_advect


