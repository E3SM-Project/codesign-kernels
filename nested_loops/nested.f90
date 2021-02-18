program nested

 implicit none

  type Timer
    integer :: start, rate=-1
  end type Timer
  

 integer :: i,j,k, n, iEdge, iCell
 integer, parameter :: nEdges = 16*1600, nVertLevels = 100, nadv = 10, ncells = 4*700
 real*8, parameter :: coef3rdOrder = 2.14
 integer :: nAdvCellsForEdge(nEdges), advCellsForEdge(nadv, nEdges)
 real*8 :: edgeFlx(nadv), tracerCur(nVertLevels,ncells), advCoefs(nadv,nEdges), &
           advCoefs3rd(nadv,nEdges), &
         wgtTmp(nVertLevels),            &! vertical temporaries for
         flxTmp(nVertLevels), sgnTmp(nVertLevels)      !   high-order flux computation

 real*8 :: coef1,coef2,coef3, csgn, normalThicknessFlux(nVertLevels,nEdges), &
            highOrderAdvectionMask(nVertLevels,nEdges), highOrderFlx(nVertLevels,nEdges), &
            tmpTracer(nadv,nVertLevels,nEdges), maxLevelCell(ncells), tmpidx
 integer :: nranks, rank
 
 logical :: new
 
type(Timer) :: crono1

 
nAdvCellsForEdge = nadv
highOrderAdvectionMask = 1
highOrderFlx = 0.0
maxLevelCell = nVertLevels

do j = 1,nEdges
do i = 1,nadv
    call random_number(tmpidx)
    advCellsForEdge(i,j) = int(ncells * tmpidx) + 1
end do
enddo

do j = 1,ncells
do i = 1,nVertLevels
    call random_number(tmpidx)
    tracerCur(i,j) = 15 * tmpidx
end do
enddo

do j = 1,nEdges
do i = 1,nVertLevels
    call random_number(tmpidx)
    normalThicknessFlux(i,j) = 15 * tmpidx
end do
enddo

do j = 1,nEdges
do i = 1,nVertLevels
do n = 1,nadv
    iCell = advCellsForEdge(n,j)
    tmpTracer(n,i,j) = tracerCur(i,iCell)
end do
end do
enddo

do j = 1,nEdges
do i = 1,nadv
    call random_number(tmpidx)
    advCoefs(i,j) = 20 * tmpidx
    call random_number(tmpidx)
    advCoefs3rd(i,j) = 20 * tmpidx
end do
enddo

    !$acc enter data copyin(normalThicknessFlux,highOrderAdvectionMask,nAdvCellsForEdge, &
    !$acc       advCellsForEdge,advCoefs,advCoefs3rd,tmpTracer,highOrderFlx,maxLevelCell, &
    !$acc       tracerCur)

#ifndef ORIG

call Tic(crono1)

do n =1,100
    !$acc parallel loop gang vector collapse(2) &
    !$acc          present(normalThicknessFlux,highOrderAdvectionMask,nAdvCellsForEdge,&
    !$acc                  advCellsForEdge,advCoefs,advCoefs3rd,tmpTracer,highOrderFlx) &
    !$acc          private(edgeFlx)
do iEdge = 1, nEdges
   do k = 1, nVertLevels
      ! Compute 3rd or 4th fluxes where requested.
      coef2 = normalThicknessFlux(k,iEdge)*highOrderAdvectionMask(k,iEdge)
      csgn = sign(1.0,normalThicknessFlux(k,iEdge))
      do i = 1, nAdvCellsForEdge(iEdge)
            coef1 = advCoefs       (i,iEdge)
            coef3 = advCoefs3rd    (i,iEdge)*coef3rdOrder
            edgeFlx(i) = tmpTracer(i,k,iEdge) * &
                   coef2 * (coef1 + coef3*csgn)
      end do ! i loop over nAdvCellsForEdge
      
      do i = 1,nAdvCellsForEdge(iEdge)
         highOrderFlx(k,iEdge) = highOrderFlx(k,iEdge) + edgeFlx(i)
      end do
   enddo
end do ! iEdge loop  
enddo

    !$acc wait
call Toc(crono1, 'nested')

#else

call Tic(crono1)

do n =1,100
    !$acc parallel loop gang vector &
    !$acc          present(normalThicknessFlux,highOrderAdvectionMask,nAdvCellsForEdge,&
    !$acc                  advCellsForEdge,advCoefs,advCoefs3rd,tracerCur,highOrderFlx) &
    !$acc          private(wgtTmp,sgnTmp,flxTmp)
do iEdge = 1, nEdges
  ! compute some common intermediate factors
  do k = 1, nVertLevels
     wgtTmp(k) = normalThicknessFlux   (k,iEdge)* &
                 highOrderAdvectionMask(k,iEdge)
     sgnTmp(k) = sign(1.0, &
                      normalThicknessFlux(k,iEdge))
     flxTmp(k) = 0.0
  end do

  ! Compute 3rd or 4th fluxes where requested.
  do i = 1, nAdvCellsForEdge(iEdge)
    iCell = advCellsForEdge(i,iEdge)
    coef1 = advCoefs       (i,iEdge)
    coef3 = advCoefs3rd    (i,iEdge)*coef3rdOrder
    do k = 1, maxLevelCell(iCell)
      flxTmp(k) = flxTmp(k) + tracerCur(k,iCell)* &
                  wgtTmp(k)*(coef1 + coef3*sgnTmp(k))
    end do ! k loop
  end do ! i loop over nAdvCellsForEdge
  
  do k=1,nVertLevels
     highOrderFlx(k,iEdge) = flxTmp(k)
  end do
enddo
enddo

    !$acc wait
call Toc(crono1, 'orig')

#endif

    !$acc exit data delete(normalThicknessFlux,highOrderAdvectionMask,nAdvCellsForEdge, &
    !$acc       advCellsForEdge,advCoefs,advCoefs3rd,tmpTracer,highOrderFlx,maxLevelCell, &
    !$acc       tracerCur)

contains

! taken from https://www.um.es/gmar/staff/gomez/computation/2016/06/28/timers.html

subroutine Tic(self)
    class (Timer), intent(inout) :: self
    integer :: start, rate

    call system_clock(count_rate=rate)
    call system_clock(start)
    self%start=start
    self%rate=rate
  end subroutine
  
  subroutine Toc(self,pref)
    class (Timer), intent(in) :: self
    character(*), intent(in) :: pref
    integer :: finish

    if(self%rate<0) then
      print*, 'Call to ''Tac'' subroutine must come after call to ''Tic'''
      stop
    endif
    call system_clock(finish)
    print*, pref, ': Elapsed time in seconds:', float(finish-self%start)/self%rate
  end subroutine
end program

