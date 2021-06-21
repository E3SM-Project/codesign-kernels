!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
! nested
!
! This program contains several forms of a particular nested loop
! kernel that is common within the MPAS ocean model. It is meant to
! test how this kernel can be optimized across CPU, GPU architectures.
!
!***********************************************************************

program nested

   implicit none
   use timerMod

   integer, parameter :: &
      RKIND = selected_real_kind(12) ! default double precision

   integer ::       &! model size variables
      nIters,       &! number of loop iterations to run
      nEdges,       &! number of edges
      nCells,       &! number of cells
      nVertLevels,  &! number of vertical levels
      nAdv           ! max num cells contributing to edge for advection

   integer :: &
      i,j,k,n,iEdge,iCell ! various loop iterators

   real (kind=RKIND), parameter :: &
      coef3rdOrder = 2.14 ! an advection parameter

   integer, dimension(:), allocatable :: &
      nAdvCellsForEdge, &! num cells contributing to each edge 
      minLevelCell,     &! min vert level at cell center
      maxLevelCell       ! max vert level at cell center

   integer, dimension(:,:), allocatable :: &
      advCellsForEdge  ! address of each cell contrib to edge

   real (kind=RKIND), dimension(:), allocatable :: &
      wgtTmp, flxTmp, sgnTmp, &! temp vert vectors for high-order flux
      edgeFlx                  ! temp for contrib to edge fluxes

   real (kind=RKIND), dimension(:,:), allocatable :: &
      tracerCur,             &! current tracer values
      normalThicknessFlux,   &! thickness flux normal to edge
      advMaskHighOrder,      &! mask for using high-order advection
      highOrderFlx,          &! temp for holding high-order adv flux
      advCoefs,              &! precomputed adv coefficients
      advCoefs3rd             ! precomputed adv coefficients 3rd-order

   real (kind=RKIND), dimension(:,:,:), allocatable :: &
      tmpTracer               ! temp for tracers contributing to edges

   real (kind=RKIND) :: &
      coef1,coef2,coef3,csgn, &! various local temps 
      randNum           ! random number
 
   type (Timer) :: &! various timers for timing loop forms
      timerOrig,   &! timer for original CPU-optimized form
      timerGPU,    &! timer for a GPU-optimized form
      timerPortbl1  ! timer for a performance portable loop form

   ! End preamble
   !-------------
   ! Begin code

   ! Set model size
   ! Read this from namelist eventually
   nIters = 100
   nEdges = 16*1600    ! number of edges
   nCells = 4*700      ! number of cells
   nVertLevels = 100   ! number of vertical levels
   nAdv = 10           ! cells contributing to an edge for advection
 
   ! Allocate various arrays
   allocate(nAdvCellsForEdge(nEdges), &
            minLevelCell    (nCells), &
            maxLevelCell    (nCells), &
            advCellsForEdge(nadv, nEdges), &
            tmpTracer     (nadv,nVertLevels,nEdges), 
            tracerCur          (nVertLevels,ncells), &
            normalThicknessFlux(nVertLevels,nEdges), &
            highOrderFlx       (nVertLevels,nEdges), &
            advMaskHighOrder(nVertLevels,nEdges), & 
            advCoefs   (nadv,nEdges), &
            advCoefs3rd(nadv,nEdges), &
            wgtTmp(nVertLevels), &
            flxTmp(nVertLevels), &
            sgnTmp(nVertLevels), &
            edgeFlx(nAdv))

   !--------------------------------------------------------------------
   ! Initialize various arrays. For this kernel test, the actual value
   ! does not matter, but the index arrays should be somewhat
   ! representative. Random indices for the horizontal are used
   ! and would be a worst case scenario. For the vertical, we
   ! want something that represents at some level a topography
   ! in which not all the cells are active.
   !--------------------------------------------------------------------

   do iCell=1,nCells
      minLevelCell(iCell) = 1 ! usual case - can be >1 for ice shelves

      ! For bottom level, we want the depth to be somewhere between
      ! 3 and the max nVertLevels with a good fraction (half?) at
      ! the max depth so...
      call random_number(randNumber)
      k = nint(randNumber*nVertLevels*2.0)
      maxLevelCell(iCell) = min(max(3,k),nVertLevels)
   end do

   ! Set tracer values to random number
   do iCell = 1,nCells
   do k = 1,nVertLevels
      if (k >= minVertLevels(iCell) .and. &
          k <= maxVertLevels(iCell)) then
         call random_number(randNumber)
         tracerCur(i,j) = 15.0_RKIND * randNumber
      else
         tracerCur(i,j) = 0.0_RKIND
      endif
   end do
   end do

   ! Create a random connectivity between cells and edges
   ! Also create random values for advective weights
   do iEdge=1,nEdges
      nAdvCellsForEdge(iEdge) = nadv
      do i = 1,nadv
         call random_number(randNumber)
         advCellsForEdge(i,iEdge) = int(nCells * randNumber) + 1
         call random_number(randNumber)
         advCoefs   (i,iEdge) = 20.0_RKIND * randNumber
         call random_number(randNumber)
         advCoefs3rd(i,iEdge) = 21.0_RKIND * randNumber
      end do
   end do
   
   ! Initialize highOrderFlx and set normalThickFlux to random val
   do iEdge = 1,nEdges
   do k = 1,nVertLevels
      highOrderFlx(k,iEdge) = 0.0_RKIND
      call random_number(randNumber)
      normalThicknessFlux(k,iEdge) = 15.0_RKIND*(0.5_RKIND - randNumber)
      advMaskHighOrder   (k,iEdge) = 1.0_RKIND
   end do
   end do

   do iEdge = 1,nEdges
   do k = 1,nVertLevels
   do n = 1,nadv
      iCell = advCellsForEdge(n,iEdge)
      tmpTracer(n,k,iEdge) = tracerCur(k,iCell)
   end do
   end do
   end do

   ! Transfer data to device
   !$acc enter data copyin(nAdvCellsForEdge, advCellsForEdge, &
   !$acc    minLevelCell, maxLevelCell, tmpTracer, tracerCur, &
   !$acc    normalThicknessFlux, highOrderFlx, advMaskHighOrder, &
   !$acc    advCoefs, advCoefs3rd, wgtTmp, flxTmp, sgnTmp, edgeFlx)

   !--------------------------------------------------------------------
   ! First loop form - original CPU optimized version
   !--------------------------------------------------------------------

   timerOrig = timerCreate('Original CPU')
   call timerStart(timerOrig)

   do n =1,nIters
      !$acc parallel loop gang vector &
      !$acc          present(normalThicknessFlux,advMaskHighOrder,nAdvCellsForEdge,&
      !$acc                  advCellsForEdge,advCoefs,advCoefs3rd,tracerCur,highOrderFlx) &
      !$acc          private(wgtTmp,sgnTmp,flxTmp)
      do iEdge = 1, nEdges
         ! compute some common intermediate factors
         do k = 1, nVertLevels
            wgtTmp(k) = normalThicknessFlux   (k,iEdge)* &
                        advMaskHighOrder(k,iEdge)
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
      end do ! edge loop
   end do ! iteration loop

   call timerStop(timerOrig)
   call timerPrint(timerOrig)

   !--------------------------------------------------------------------
   ! Second loop form - parallel and GPU optimized
   !--------------------------------------------------------------------

   timerGPU  = timerCreate('GPU Optimized')
   call timerStart(timerGPU)

   do n =1,nIters
      !$acc parallel loop gang vector collapse(2) &
      !$acc          present(normalThicknessFlux,advMaskHighOrder,nAdvCellsForEdge,&
      !$acc                  advCellsForEdge,advCoefs,advCoefs3rd,tmpTracer,highOrderFlx) &
      !$acc          private(edgeFlx)
      do iEdge = 1, nEdges
      do k = 1, nVertLevels
         ! Compute 3rd or 4th fluxes where requested.
         coef2 = normalThicknessFlux(k,iEdge)*advMaskHighOrder(k,iEdge)
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
      end do ! vertical loop
      end do ! iEdge loop  
   end do ! iteration loop

   call timerStop(timerGPU)
   call timerPrint(timerGPU)

   !--------------------------------------------------------------------
   !timerPortbl1 = timerCreate('Portable form 1')

   !--------------------------------------------------------------------
   ! Clean up and deallocate
   !--------------------------------------------------------------------

   !$acc exit data delete(nAdvCellsForEdge, advCellsForEdge, &
   !$acc    minLevelCell, maxLevelCell, tmpTracer, tracerCur, &
   !$acc    normalThicknessFlux, highOrderFlx, advMaskHighOrder, &
   !$acc    advCoefs, advCoefs3rd, wgtTmp, flxTmp, sgnTmp, edgeFlx)

   deallocate(nAdvCellsForEdge, advCellsForEdge, minLevelCell, &
              maxLevelCell, tmpTracer, tracerCur, &
              normalThicknessFlux, highOrderFlx, advMaskHighOrder, &
              advCoefs, advCoefs3rd, wgtTmp, flxTmp, sgnTmp, edgeFlx)

!***********************************************************************

end program nested

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
