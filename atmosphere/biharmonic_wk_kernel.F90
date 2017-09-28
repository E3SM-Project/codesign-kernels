


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! HELPER ROUTINES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module kernel_utils
  implicit none
  public
  integer             , parameter :: real_kind = selected_real_kind(13)
  integer             , parameter :: np = 4
  integer             , parameter :: nlev = 72
  integer             , parameter :: qsize = 40
  real(kind=real_kind), parameter :: rrearth = 0.00000016666666666666
  integer :: nelemd = 16
  integer :: nets = 1
  integer :: nete = 16

  type :: derivative_t
    real(kind=real_kind) :: Dvv(np,np)
  end type derivative_t

  type :: element_t
     real (kind=real_kind) :: Dinv(np,np,2,2)
     real (kind=real_kind) :: spheremp(np,np)
     real (kind=real_kind) :: tensorVisc(np,np,2,2)
  end type element_t

  type(derivative_t) :: deriv
  type(element_t), allocatable :: elem(:)
  real(kind=real_kind), allocatable :: qtens(:,:,:,:,:)
  real(kind=real_kind), allocatable :: grads(:,:,:,:,:,:)
  real(kind=real_kind), allocatable :: qtens_save(:,:,:,:,:)

contains

  subroutine allocate_data()
    implicit none
    !Allocate and initialize data
    allocate(elem(nelemd))
    allocate(qtens(np,np,nlev,qsize,nelemd))
    allocate(grads(np,np,2,nlev,qsize,nelemd))
    allocate(qtens_save(np,np,nlev,qsize,nelemd))
  end subroutine allocate_data



  subroutine initialize_data()
    implicit none
    integer :: ie
    call myrandom( product(shape(deriv%dvv)) , deriv%dvv , reset = .true. )
    do ie = nets , nete
      call myrandom( product(shape(elem(ie)%Dinv      )) , elem(ie)%Dinv       )
      call myrandom( product(shape(elem(ie)%spheremp  )) , elem(ie)%spheremp   )
      call myrandom( product(shape(elem(ie)%tensorVisc)) , elem(ie)%tensorVisc )
    enddo
    call myrandom( product(shape(qtens)) , qtens )
  end subroutine initialize_data



  subroutine save_qtens()
    implicit none
    qtens_save = qtens
  end subroutine save_qtens



  function compute_l2norm()  result(l2)
    implicit none
    real(kind=real_kind) :: l2
    l2 = sqrt( sum( (qtens - qtens_save)**2 ) / sum(qtens_save**2) )
  end function compute_l2norm



  subroutine myrandom(n,a,reset)
    implicit none
    integer             , intent(in   ) :: n
    real(kind=real_kind), intent(  out) :: a(n)
    logical, optional   , intent(in   ) :: reset
    integer :: i
    integer, save :: old = 11
    if (present(reset)) then
      if (reset) old = 11
    endif
    do i = 1 , n
      old = MOD((1301*old+97), 1024*128)
      a(i) = old / dble(1024*128)
    enddo
  end subroutine myrandom

end module kernel_utils



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! CPU KERNEL ROUTINES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module biharmonic_wk_scalar_cpu
  use kernel_utils
  implicit none
  private

  public :: biharmonic_wk_scalar

contains

  function gradient_sphere(s,deriv,Dinv) result(ds)
    real(kind=real_kind), intent(in) :: s   (np,np)
    type (derivative_t) , intent(in) :: deriv
    real(kind=real_kind), intent(in) :: Dinv(np,np,2,2)
    real(kind=real_kind)             :: ds  (np,np,2)
    integer :: i, j, l
    real(kind=real_kind) :: dsdx00, dsdy00, v1(np,np) , v2(np,np)
    do j = 1 , np
      do l = 1 , np
        dsdx00 = 0.0d0
        dsdy00 = 0.0d0
        do i=1,np
          dsdx00 = dsdx00 + deriv%Dvv(i,l)*s(i,j)
          dsdy00 = dsdy00 + deriv%Dvv(i,l)*s(j,i)
        end do
        v1(l,j) = dsdx00*rrearth
        v2(j,l) = dsdy00*rrearth
      enddo
    enddo
    do j=1,np
      do i=1,np
        ds(i,j,1) = Dinv(i,j,1,1)*v1(i,j) + Dinv(i,j,2,1)*v2(i,j)
        ds(i,j,2) = Dinv(i,j,1,2)*v1(i,j) + Dinv(i,j,2,2)*v2(i,j)
      enddo
    enddo
  end function gradient_sphere



  function divergence_sphere_wk(v,deriv,elem) result(div)
    real(kind=real_kind), intent(in) :: v  (np,np,2)
    type (derivative_t) , intent(in) :: deriv
    type (element_t)    , intent(in) :: elem
    real(kind=real_kind)             :: div(np,np)
    real(kind=real_kind) :: vtemp(np,np,2)
    integer i,j,m,n
    do j = 1 , np
      do i = 1 , np
        vtemp(i,j,1) = (elem%Dinv(i,j,1,1)*v(i,j,1) + elem%Dinv(i,j,1,2)*v(i,j,2))
        vtemp(i,j,2) = (elem%Dinv(i,j,2,1)*v(i,j,1) + elem%Dinv(i,j,2,2)*v(i,j,2))
      enddo
    enddo
    do n = 1 , np
      do m = 1 , np
        div(m,n) = 0
        do j = 1 , np
          div(m,n) = div(m,n) - ( elem%spheremp(j,n)*vtemp(j,n,1)*deriv%Dvv(m,j) + &
                                  elem%spheremp(m,j)*vtemp(m,j,2)*deriv%Dvv(n,j) ) * rrearth
        enddo
      enddo
    enddo
  end function divergence_sphere_wk



  function laplace_sphere_wk(s,deriv,elem) result(laplace)
    real(kind=real_kind), intent(in   ) :: s(np,np) 
    type (derivative_t) , intent(in   ) :: deriv
    type (element_t)    , intent(in   ) :: elem
    real(kind=real_kind) :: laplace(np,np)
    real(kind=real_kind) :: grads(np,np,2), oldgrads(np,np,2)
    integer :: i,j
    grads = gradient_sphere(s,deriv,elem%Dinv)
    oldgrads = grads
    do j = 1 , np
      do i = 1 , np
        grads(i,j,1) = oldgrads(i,j,1)*elem%tensorVisc(i,j,1,1) + &
                       oldgrads(i,j,2)*elem%tensorVisc(i,j,1,2)
        grads(i,j,2) = oldgrads(i,j,1)*elem%tensorVisc(i,j,2,1) + &
                       oldgrads(i,j,2)*elem%tensorVisc(i,j,2,2)
      enddo
    enddo
    laplace = divergence_sphere_wk(grads,deriv,elem)
  end function laplace_sphere_wk



  subroutine biharmonic_wk_scalar(elem,qtens,deriv,nets,nete)
    type (element_t)     , intent(in   ) :: elem(:)
    real (kind=real_kind), intent(inout) :: qtens(np,np,nlev,qsize,nets:nete)
    type (derivative_t)  , intent(in   ) :: deriv
    integer              , intent(in   ) :: nets
    integer              , intent(in   ) :: nete
    integer :: k,i,j,ie,q
    do ie = nets , nete
      do q = 1 , qsize
        do k = 1 , nlev
          qtens(:,:,k,q,ie) = laplace_sphere_wk(qtens(:,:,k,q,ie),deriv,elem(ie))
        enddo
      enddo
    enddo
  end subroutine biharmonic_wk_scalar

end module biharmonic_wk_scalar_cpu



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! GPU KERNEL ROUTINES USING COMPILER INLINING
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module biharmonic_wk_scalar_gpu_compiler_inline
  use kernel_utils
  implicit none
  private
  integer, parameter :: kchunk = 8

  public :: biharmonic_wk_scalar

contains

  subroutine gradient_sphere(s,deriv,Dinv,ds)
    implicit none
    !$acc routine vector
    real(kind=real_kind), intent(in   ) :: s   (np,np,kchunk)
    type (derivative_t) , intent(in   ) :: deriv
    real(kind=real_kind), intent(in   ) :: Dinv(np,np,2,2)
    real(kind=real_kind), intent(  out) :: ds  (np,np,kchunk,2)
    integer :: i, j, l, kk
    real(kind=real_kind) :: dsdx00, dsdy00, v1(np,np) , v2(np,np)
    !$acc loop vector collapse(3) private(dsdx00,dsdy00)
    do kk = 1 , kchunk
      do j = 1 , np
        do i = 1 , np
          dsdx00 = 0.0d0
          dsdy00 = 0.0d0
          do l = 1 , np
            dsdx00 = dsdx00 + deriv%Dvv(l,i)*s(l,j,kk)
            dsdy00 = dsdy00 + deriv%Dvv(l,j)*s(i,l,kk)
          end do
          ds(i,j,kk,1) = ( Dinv(i,j,1,1)*dsdx00 + Dinv(i,j,2,1)*dsdy00 ) * rrearth
          ds(i,j,kk,2) = ( Dinv(i,j,1,2)*dsdx00 + Dinv(i,j,2,2)*dsdy00 ) * rrearth
        enddo
      enddo
    enddo
  end subroutine gradient_sphere



  subroutine divergence_sphere_wk(v,deriv,elem,div)
    implicit none
    !$acc routine vector
    real(kind=real_kind), intent(in   ) :: v  (np,np,kchunk,2)
    type (derivative_t) , intent(in   ) :: deriv
    type (element_t)    , intent(in   ) :: elem
    real(kind=real_kind), intent(  out) :: div(np,np,kchunk)
    real(kind=real_kind) :: vtemp(np,np,kchunk,2), tmp
    integer i,j,m,kk
    !$acc loop vector collapse(3)
    do kk = 1 , kchunk
      do j = 1 , np
        do i = 1 , np
          vtemp(i,j,kk,1) = (elem%Dinv(i,j,1,1)*v(i,j,kk,1) + elem%Dinv(i,j,1,2)*v(i,j,kk,2)) * elem%spheremp(i,j) * rrearth
          vtemp(i,j,kk,2) = (elem%Dinv(i,j,2,1)*v(i,j,kk,1) + elem%Dinv(i,j,2,2)*v(i,j,kk,2)) * elem%spheremp(i,j) * rrearth
        enddo
      enddo
    enddo
    !$acc loop vector collapse(3) private(tmp)
    do kk = 1 , kchunk
      do j = 1 , np
        do i = 1 , np
          tmp = 0
          do m = 1 , np
            tmp = tmp - vtemp(m,j,kk,1)*deriv%Dvv(i,m) - &
                        vtemp(i,m,kk,2)*deriv%Dvv(j,m)
          enddo
          div(i,j,kk) = tmp
        enddo
      enddo
    enddo
  end subroutine divergence_sphere_wk



  subroutine laplace_sphere_wk(s,deriv,elem,laplace)
    implicit none
    !$acc routine vector
    real(kind=real_kind), intent(in   ) :: s(np,np,kchunk) 
    type (derivative_t) , intent(in   ) :: deriv
    type (element_t)    , intent(in   ) :: elem
    real(kind=real_kind), intent(inout) :: laplace(np,np,kchunk)
    real(kind=real_kind) :: grads(np,np,kchunk,2), oldgrads(np,np,kchunk,2)
    integer :: i,j,kk
    call gradient_sphere(s,deriv,elem%Dinv,grads)
    !$acc loop vector collapse(3)
    do kk = 1 , kchunk
      do j = 1 , np
        do i = 1 , np
          oldgrads(i,j,kk,1) = grads(i,j,kk,1)
          oldgrads(i,j,kk,2) = grads(i,j,kk,2)
        enddo
      enddo
    enddo
    !$acc loop vector collapse(3)
    do kk = 1 , kchunk
      do j = 1 , np
        do i = 1 , np
          grads(i,j,kk,1) = oldgrads(i,j,kk,1)*elem%tensorVisc(i,j,1,1) + &
                            oldgrads(i,j,kk,2)*elem%tensorVisc(i,j,1,2)
          grads(i,j,kk,2) = oldgrads(i,j,kk,1)*elem%tensorVisc(i,j,2,1) + &
                            oldgrads(i,j,kk,2)*elem%tensorVisc(i,j,2,2)
        enddo
      enddo
    enddo
    call divergence_sphere_wk(grads,deriv,elem,laplace)
  end subroutine laplace_sphere_wk



  subroutine biharmonic_wk_scalar(elem,qtens,deriv,nets,nete)
    implicit none
    type (element_t)     , intent(in   ) :: elem(:)
    real (kind=real_kind), intent(inout) :: qtens(np,np,nlev,qsize,nets:nete)
    type (derivative_t)  , intent(in   ) :: deriv
    integer              , intent(in   ) :: nets
    integer              , intent(in   ) :: nete
    integer :: k,i,j,ie,q,kc,kk
    real(kind=real_kind) :: qtmp(np,np,kchunk)
    !$acc parallel loop gang collapse(3) private(qtmp) present(qtens,elem,deriv)
    do ie = nets , nete
      do q = 1 , qsize
        !Manually fissioning the klev loop & pushing a block of it down the callstack
        do kc = 1 , nlev/kchunk
          !$acc cache(qtmp)

          !First, load data to the temp variable (for shared memory)
          !$acc loop vector collapse(3) private(k)
          do kk = 1 , kchunk
            do j = 1 , np
              do i = 1 , np
                k = (kc-1)*kchunk + kk
                if (k <= nlev) qtmp(i,j,kk) = qtens(i,j,k,q,ie)
              enddo
            enddo
          enddo

          !Call laplace_sphere_wk only over a block of vertical levels
          call laplace_sphere_wk(qtmp,deriv,elem(ie),qtmp)

          !Finally, load the data back
          !$acc loop vector collapse(3) private(k)
          do kk = 1 , kchunk
            do j = 1 , np
              do i = 1 , np
                k = (kc-1)*kchunk + kk
                if (k <= nlev) qtens(i,j,k,q,ie) = qtmp(i,j,kk)
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
  end subroutine biharmonic_wk_scalar

end module biharmonic_wk_scalar_gpu_compiler_inline



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! GPU KERNEL ROUTINES PUSHING LOOPING DOWN THE CALLSTACK
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module biharmonic_wk_scalar_gpu_push_loop
  use kernel_utils
  implicit none
  private
  integer, parameter :: kchunk = 8

  public :: biharmonic_wk_scalar

contains



  subroutine divergence_sphere_wk_openacc(v,deriv,elem,div,len,nets,nete,ntl,tl)
    implicit none
!   input:  v = velocity in lat-lon coordinates
!   ouput:  div(v)  spherical divergence of v, integrated by parts
!   Computes  -< grad(psi) dot v > 
!   (the integrated by parts version of < psi div(v) > )
!   note: after DSS, divergence_sphere () and divergence_sphere_wk() 
!   are identical to roundoff, as theory predicts.
    real(kind=real_kind), intent(in) :: v(np,np,2,len,ntl,nelemd)  ! in lat-lon coordinates
    type (derivative_t) , intent(in) :: deriv
    type (element_t)    , intent(in) :: elem(:)
    real(kind=real_kind), intent(out):: div(np,np,len,ntl,nelemd)
    integer             , intent(in) :: len
    integer             , intent(in) :: nets , nete , ntl , tl
    ! Local
    integer, parameter :: kchunk = 8
    integer :: i,j,l,k,ie,kc,kk
    real(kind=real_kind) :: vtemp(np,np,2,kchunk), tmp, deriv_tmp(np,np)
    ! latlon- > contra
    !$acc parallel loop gang collapse(2) present(v,elem(:),div,deriv) private(vtemp,deriv_tmp)
    do ie = nets , nete
      do kc = 1 , len/kchunk+1
        !$acc cache(vtemp,deriv_tmp)
        !$acc loop vector collapse(3)
        do kk = 1 , kchunk
          do j = 1 , np
            do i = 1 , np
              k = (kc-1)*kchunk+kk
              if (k <= len) then
                vtemp(i,j,1,kk)=elem(ie)%spheremp(i,j)*(elem(ie)%Dinv(i,j,1,1)*v(i,j,1,k,tl,ie) + elem(ie)%Dinv(i,j,1,2)*v(i,j,2,k,tl,ie))
                vtemp(i,j,2,kk)=elem(ie)%spheremp(i,j)*(elem(ie)%Dinv(i,j,2,1)*v(i,j,1,k,tl,ie) + elem(ie)%Dinv(i,j,2,2)*v(i,j,2,k,tl,ie))
              endif
              if (kk == 1) deriv_tmp(i,j) = deriv%Dvv(i,j)
            enddo
          enddo
        enddo
        !$acc loop vector collapse(3) private(tmp)
        do kk = 1 , kchunk
          do j = 1 , np
            do i = 1 , np
              k = (kc-1)*kchunk+kk
              if (k <= len) then
                tmp = 0.
                do l = 1 , np
                  tmp = tmp - ( vtemp(l,j,1,kk)*deriv_tmp(i,l) + vtemp(i,l,2,kk)*deriv_tmp(j,l) )
                enddo
                div(i,j,k,tl,ie) = tmp * rrearth
              endif
            enddo
          enddo
        enddo
      enddo
    enddo
  end subroutine divergence_sphere_wk_openacc



  subroutine gradient_sphere_openacc(s,deriv,elem,ds,len,nets,nete,ntl,tl)
    implicit none
    !   input s:  scalar
    !   output  ds: spherical gradient of s, lat-lon coordinates
    real(kind=real_kind), intent(in) :: s(np,np,len,ntl,nelemd)
    type(derivative_t)  , intent(in) :: deriv
    type(element_t)     , intent(in) :: elem(:)
    real(kind=real_kind), intent(out):: ds(np,np,2,len,ntl,nelemd)
    integer             , intent(in) :: len
    integer             , intent(in) :: nets,nete,ntl,tl
    integer, parameter :: kchunk = 8
    integer :: i, j, l, k, ie, kc, kk
    real(kind=real_kind) :: dsdx00, dsdy00
    real(kind=real_kind) :: stmp(np,np,kchunk), deriv_tmp(np,np)
    !$acc parallel loop gang collapse(2) present(ds,elem(:),s,deriv%Dvv) private(stmp,deriv_tmp)
    do ie = nets , nete
      do kc = 1 , len/kchunk+1
        !$acc cache(stmp,deriv_tmp)
        !$acc loop vector collapse(3)
        do kk = 1 , kchunk
          do j = 1 , np
            do i = 1 , np
              k = (kc-1)*kchunk+kk
              if (k > len) k = len
              stmp(i,j,kk) = s(i,j,k,tl,ie)
              if (kk == 1) deriv_tmp(i,j) = deriv%Dvv(i,j)
            enddo
          enddo
        enddo
        !$acc loop vector collapse(3)
        do kk = 1 , kchunk
          do j = 1 , np
            do i = 1 , np
              k = (kc-1)*kchunk+kk
              if (k <= len) then
                dsdx00=0.0d0
                dsdy00=0.0d0
                do l = 1 , np
                  dsdx00 = dsdx00 + deriv_tmp(l,i)*stmp(l,j,kk)
                  dsdy00 = dsdy00 + deriv_tmp(l,j)*stmp(i,l,kk)
                enddo
                ds(i,j,1,k,tl,ie) = ( elem(ie)%Dinv(i,j,1,1)*dsdx00 + elem(ie)%Dinv(i,j,2,1)*dsdy00 ) * rrearth
                ds(i,j,2,k,tl,ie) = ( elem(ie)%Dinv(i,j,1,2)*dsdx00 + elem(ie)%Dinv(i,j,2,2)*dsdy00 ) * rrearth
              endif
            enddo
          enddo
        enddo
      enddo
    enddo
  end subroutine gradient_sphere_openacc



  subroutine laplace_sphere_wk_openacc(s,grads,deriv,elem,laplace,len,nets,nete,ntl,tl)
    implicit none
    !input:  s = scalar
    !ouput:  -< grad(PHI), grad(s) >   = weak divergence of grad(s)
    !note: for this form of the operator, grad(s) does not need to be made C0
    real(kind=real_kind) , intent(in   ) :: s(np,np,len,ntl,nelemd)
    real(kind=real_kind) , intent(inout) :: grads(np,np,2,len,nelemd)
    type (derivative_t)  , intent(in   ) :: deriv
    type (element_t)     , intent(in   ) :: elem(:)
    real(kind=real_kind) , intent(  out) :: laplace(np,np,len,ntl,nelemd)
    integer              , intent(in   ) :: len,nets,nete,ntl,tl
    integer :: i,j,k,ie
    ! Local
    real(kind=real_kind) :: oldgrads(2)
    call gradient_sphere_openacc(s,deriv,elem(:),grads,len,nets,nete,ntl,tl)
    !$acc parallel loop gang vector collapse(4) present(grads,elem(:)) private(oldgrads)
    do ie = nets , nete
      do k = 1 , len
        do j = 1 , np
          do i = 1 , np
            oldgrads = grads(i,j,:,k,ie)
            grads(i,j,1,k,ie) = sum(oldgrads(:)*elem(ie)%tensorVisc(i,j,1,:))
            grads(i,j,2,k,ie) = sum(oldgrads(:)*elem(ie)%tensorVisc(i,j,2,:))
          enddo
        enddo
      enddo
    enddo
    ! note: divergnece_sphere and divergence_sphere_wk are identical *after* bndry_exchange
    ! if input is C_0.  Here input is not C_0, so we should use divergence_sphere_wk().  
    call divergence_sphere_wk_openacc(grads,deriv,elem(:),laplace,len,nets,nete,ntl,tl)
  end subroutine laplace_sphere_wk_openacc



  subroutine biharmonic_wk_scalar(elem,qtens,grads,deriv,nets,nete)
    implicit none
    type (element_t)     , intent(inout) :: elem(:)
    real (kind=real_kind), intent(inout) :: qtens(np,np,nlev,qsize,nelemd)
    real(kind=real_kind) , intent(inout) :: grads(np,np,2,nlev,qsize,nelemd)
    type (derivative_t)  , intent(in   ) :: deriv
    integer              , intent(in   ) :: nets,nete
    integer :: k,kptr,i,j,ie,ic,q
    call laplace_sphere_wk_openacc(qtens,grads,deriv,elem,qtens,nlev*qsize,nets,nete,1,1)
  end subroutine biharmonic_wk_scalar

end module biharmonic_wk_scalar_gpu_push_loop





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! MAIN KERNEL DRIVER
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program biharmonic_wk_scalar_kernel
  use kernel_utils
  use biharmonic_wk_scalar_cpu                , only: bws_cpu  => biharmonic_wk_scalar
  use biharmonic_wk_scalar_gpu_compiler_inline, only: bws_gpu1 => biharmonic_wk_scalar
  use biharmonic_wk_scalar_gpu_push_loop      , only: bws_gpu2 => biharmonic_wk_scalar
  implicit none
  integer*8 :: t1, t2, rt
  call allocate_data()
  !$acc enter data pcreate(qtens,elem,deriv,grads)

  call initialize_data()
  call system_clock(t1)
  call bws_cpu(elem,qtens,deriv,nets,nete)
  call system_clock(t2,rt)
  call save_qtens()
  write(*,*) 'CPU  time: ', dble(t2-t1)/dble(rt)

  call initialize_data()
  call system_clock(t1)
  !$acc update device(qtens,elem,deriv)
  call bws_gpu1(elem,qtens,deriv,nets,nete)
  !$acc update host(qtens)
  call system_clock(t2,rt)
  write(*,*) 'GPU1 time: ', dble(t2-t1)/dble(rt)
  write(*,*) 'GPU1 L2 norm: ', compute_l2norm()

  call initialize_data()
  call system_clock(t1)
  !$acc update device(qtens,elem,deriv)
  call bws_gpu2(elem,qtens,grads,deriv,nets,nete)
  !$acc update host(qtens)
  call system_clock(t2,rt)
  write(*,*) 'GPU2 time: ', dble(t2-t1)/dble(rt)
  write(*,*) 'GPU2 L2 norm: ', compute_l2norm()

end program biharmonic_wk_scalar_kernel




