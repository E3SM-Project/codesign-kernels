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

   use timerMod
#ifndef NO_MPI
   use mpi
#endif
#ifdef USE_CKE
   use cke_mod
#endif
   implicit none

   integer, parameter :: &
      RKIND = selected_real_kind(12) ! default double precision

   integer ::       &! model size variables
      ierr,         &! error flag
      nIters,       &! number of loop iterations to run
      nEdges,       &! number of edges
      nCells,       &! number of cells
      nVertLevels,  &! number of vertical levels
      nAdv           ! max num cells contributing to edge for advection

   namelist /nested_nml/ nIters, nEdges, nCells, nVertLevels, nAdv

   integer :: &
      i,k,n,iEdge,iCell, &! various loop iterators
      kmin, kmax            ! loop limits

   real (kind=RKIND), parameter :: &
      coef3rdOrder = 2.14, &! an advection parameter
      errTol       = 1.d-10 ! error tolerance

   integer, dimension(:), allocatable :: &
      nAdvCellsForEdge, &! num cells contributing to each edge 
      minLevelCell,     &! min vert level at cell center
      maxLevelCell       ! max vert level at cell center

   integer, dimension(:,:), allocatable :: &
      advCellsForEdge  ! address of each cell contrib to edge

   real (kind=RKIND) :: &
      relErr,           &! relative error in result
      refVal,           &! reference value for correctness
      edgeFlx            ! temp for contrib to edge fluxes

#ifdef USE_OMPOFFLOAD
   ! current OpenMP offload doesn't support array private so
   ! need to promote these temps
   real (kind=RKIND), dimension(:,:), allocatable :: &
      wgtTmp, sgnTmp   ! temp vert vectors for high-order flux
#else
   real (kind=RKIND), dimension(:), allocatable :: &
      wgtTmp, sgnTmp, flxTmp   ! temp vert vectors for high-order flux
#endif

   real (kind=RKIND), dimension(:,:), allocatable :: &
      tracerCur,             &! current tracer values
      normalThicknessFlux,   &! thickness flux normal to edge
      advMaskHighOrder,      &! mask for using high-order advection
      cellMask,              &! mask for active cells
      highOrderFlx,          &! temp for holding high-order adv flux
      refFlx,                &! reference value for high-order adv flux
      advCoefs,              &! precomputed adv coefficients
      advCoefs3rd             ! precomputed adv coefficients 3rd-order

   real (kind=RKIND) :: &
      coef1,coef2,coef3,csgn, &! various local temps 
      randNum           ! random number
 
   type (Timer) :: &! various timers for timing loop forms
      timerData,   &! timer for data transfers to device
      timerOrig,   &! timer for original CPU-optimized form
      timerGPU,    &! timer for a GPU-optimized form
      timerGPU2, timer_cke

#ifdef F90_PACK_SIZE
   integer, parameter :: packn = F90_PACK_SIZE
#else
   integer, parameter :: packn = 1
#endif
   integer :: nvlpk, nvldim, k0
   real(RKIND) :: csgnpk(packn), coef2pk(packn), edgeFlxPk(packn)
#define pkslc(k0) ((k0)+1):((k0)+packn)
   
   ! End preamble
   !-------------
   ! Begin code

#ifndef NO_MPI
   call MPI_Init(ierr)
#endif

   ! Read model size info from namelist input, overwriting defaults
   open (15, file='nested.nml')
   read (15, nml=nested_nml)
   close(15)

   nvlpk = (nVertLevels + packn - 1)/packn
   nvldim = nvlpk*packn

   ! Allocate various arrays
   allocate(nAdvCellsForEdge(nEdges), &
            minLevelCell    (nCells), &
            maxLevelCell    (nCells), &
            advCellsForEdge(nadv, nEdges), &
            tracerCur          (nvldim,ncells), &
            cellMask           (nvldim,ncells), &
            normalThicknessFlux(nvldim,nEdges), &
            highOrderFlx       (nvldim,nEdges), &
            refFlx             (nvldim,nEdges), &
            advMaskHighOrder(nvldim,nEdges), & 
            advCoefs   (nadv,nEdges), &
            advCoefs3rd(nadv,nEdges))
#ifdef USE_OMPOFFLOAD
   ! current OMP offload doesn't support private arrays so
   ! need to promote these to 2-d arrays
   allocate(wgtTmp(nVertLevels,nEdges), &
            sgnTmp(nVertLevels,nEdges))
#else
   allocate(wgtTmp(nVertLevels), &
            sgnTmp(nVertLevels), &
            flxTmp(nVertLevels))
#endif

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
      call random_number(randNum)
      k = nint(randNum*nVertLevels*2.0)
      maxLevelCell(iCell) = min(max(3,k),nVertLevels)
   end do

   ! Set tracer values to random number
   do iCell = 1,nCells
   do k = 1,nVertLevels
      if (k >= minLevelCell(iCell) .and. &
          k <= maxLevelCell(iCell)) then
         call random_number(randNum)
         tracerCur(k,iCell) = 15.0_RKIND * randNum
         cellMask (k,iCell) = 1.0_RKIND
      else
         tracerCur(k,iCell) = 0.0_RKIND
         cellMask (k,iCell) = 0.0_RKIND
      endif
   end do
   end do

   ! Create a random connectivity between cells and edges
   ! Also create random values for advective weights
   do iEdge=1,nEdges
      nAdvCellsForEdge(iEdge) = nadv
      do i = 1,nadv
         call random_number(randNum)
         advCellsForEdge(i,iEdge) = int(nCells * randNum) + 1
         call random_number(randNum)
         advCoefs   (i,iEdge) = 20.0_RKIND * randNum
         call random_number(randNum)
         advCoefs3rd(i,iEdge) = 21.0_RKIND * randNum
      end do
   end do
   
   ! Initialize highOrderFlx and set normalThickFlux to random val
   do iEdge = 1,nEdges
   do k = 1,nVertLevels
      highOrderFlx(k,iEdge) = 0.0_RKIND
      call random_number(randNum)
      normalThicknessFlux(k,iEdge) = 15.0_RKIND*(0.5_RKIND - randNum)
      advMaskHighOrder   (k,iEdge) = 1.0_RKIND
   end do
   end do

   !--------------------------------------------------------------------
   ! Use original form on CPU for the reference value to check
   ! correctness
   !--------------------------------------------------------------------

   do iEdge = 1, nEdges
      ! compute some common intermediate factors
#ifdef USE_OMPOFFLOAD
      do k = 1, nVertLevels
         wgtTmp(k,iEdge) = normalThicknessFlux   (k,iEdge)* &
                     advMaskHighOrder(k,iEdge)
         sgnTmp(k,iEdge) = sign(1.0_RKIND, &
                          normalThicknessFlux(k,iEdge))
         refFlx(k,iEdge) = 0.0_RKIND
      end do
#else
      do k = 1, nVertLevels
         wgtTmp(k) = normalThicknessFlux   (k,iEdge)* &
                     advMaskHighOrder(k,iEdge)
         sgnTmp(k) = sign(1.0_RKIND, &
                          normalThicknessFlux(k,iEdge))
         refFlx(k,iEdge) = 0.0_RKIND
      end do
#endif

      ! Compute 3rd or 4th fluxes where requested.
      do i = 1, nAdvCellsForEdge(iEdge)
         iCell = advCellsForEdge(i,iEdge)
         kmin  = minLevelCell(iCell)
         kmax  = maxLevelCell(iCell)
         coef1 = advCoefs       (i,iEdge)
         coef3 = advCoefs3rd    (i,iEdge)*coef3rdOrder
         do k = kmin, kmax
#ifdef USE_OMPOFFLOAD
            refFlx(k,iEdge) = refFlx(k,iEdge) + tracerCur(k,iCell)* &
                        wgtTmp(k,iEdge)*(coef1 + coef3*sgnTmp(k,iEdge))
#else
            refFlx(k,iEdge) = refFlx(k,iEdge) + tracerCur(k,iCell)* &
                        wgtTmp(k)*(coef1 + coef3*sgnTmp(k))
#endif
         end do ! k loop
      end do ! i loop over nAdvCellsForEdge
  
   end do ! edge loop

   !--------------------------------------------------------------------
   ! Transfer data to device
   !--------------------------------------------------------------------

   timerData = timerCreate('Data transfer')
#ifdef USE_CKE_dont_skip
   ! if above is just USE_CKE, then skip all but the C++ perf test
   print *,'skip F90 code'
#else
   call timerStart(timerData)

#ifdef USE_OPENACC
   !$acc enter data copyin(nAdvCellsForEdge, advCellsForEdge, &
   !$acc    minLevelCell, maxLevelCell, tracerCur, cellMask, &
   !$acc    normalThicknessFlux, highOrderFlx, advMaskHighOrder, &
   !$acc    advCoefs, advCoefs3rd, wgtTmp, sgnTmp, flxTmp)
#endif
#ifdef USE_OMPOFFLOAD
   !$omp target enter data map(to:nAdvCellsForEdge, advCellsForEdge, &
   !$omp    minLevelCell, maxLevelCell, tracerCur, cellMask, &
   !$omp    normalThicknessFlux, highOrderFlx, advMaskHighOrder, &
   !$omp    advCoefs, advCoefs3rd, wgtTmp, sgnTmp)
#endif
   call timerStop(timerData)

   !--------------------------------------------------------------------
   ! First loop form - original CPU optimized version
   !--------------------------------------------------------------------

   timerOrig = timerCreate('Original CPU')
   call timerStart(timerOrig)

   do n =1,nIters
#ifdef USE_OPENACC
      !$acc parallel loop gang vector &
      !$acc    present(normalThicknessFlux, advMaskHighOrder, &
      !$acc            nAdvCellsForEdge, advCellsForEdge, &
      !$acc            advCoefs, advCoefs3rd, tracerCur, highOrderFlx) &
      !$acc    private(wgtTmp, sgnTmp, flxTmp)
#endif
#ifdef USE_OMPOFFLOAD
      !$omp target teams &
      !$omp    map(to: normalThicknessFlux, advMaskHighOrder, &
      !$omp            nAdvCellsForEdge, advCellsForEdge, &
      !$omp            advCoefs, advCoefs3rd, tracerCur, &
      !$omp            wgtTmp, sgnTmp) &
      !$omp    map(from: highOrderFlx)
      !$omp distribute parallel do
#endif
      do iEdge = 1, nEdges
         ! compute some common intermediate factors
#ifdef USE_OMPOFFLOAD
         do k = 1, nVertLevels
            wgtTmp(k,iEdge) = normalThicknessFlux   (k,iEdge)* &
                        advMaskHighOrder(k,iEdge)
            sgnTmp(k,iEdge) = sign(1.0_RKIND, &
                             normalThicknessFlux(k,iEdge))
            highOrderFlx(k,iEdge) = 0.0_RKIND
         end do
#else
         do k = 1, nVertLevels
            wgtTmp(k) = normalThicknessFlux   (k,iEdge)* &
                        advMaskHighOrder(k,iEdge)
            sgnTmp(k) = sign(1.0_RKIND, &
                             normalThicknessFlux(k,iEdge))
            flxTmp(k) = 0.0_RKIND
         end do
#endif

         ! Compute 3rd or 4th fluxes where requested.
         do i = 1, nAdvCellsForEdge(iEdge)
            iCell = advCellsForEdge(i,iEdge)
            kmin  = minLevelCell(iCell)
            kmax  = maxLevelCell(iCell)
            coef1 = advCoefs       (i,iEdge)
            coef3 = advCoefs3rd    (i,iEdge)*coef3rdOrder
#ifdef USE_OMPOFFLOAD
            do k = kmin, kmax
               highOrderFlx(k,iEdge) = highOrderFlx(k,iEdge) + tracerCur(k,iCell)* &
                           wgtTmp(k,iEdge)*(coef1 + coef3*sgnTmp(k,iEdge))
            end do ! k loop
#else
            do k = kmin, kmax
               !highOrderFlx(k,iEdge) = highOrderFlx(k,iEdge) + tracerCur(k,iCell)* &
               !            wgtTmp(k)*(coef1 + coef3*sgnTmp(k))
               flxTmp(k) = flxTmp(k) + tracerCur(k,iCell)* &
                           wgtTmp(k)*(coef1 + coef3*sgnTmp(k))
            enddo ! k loop
#endif
         end do ! i loop over nAdvCellsForEdge
#ifndef USE_OMPOFFLOAD
         do k = 1,nVertLevels
            highOrderFlx(k,iEdge) = flxTmp(k)
         end do
#endif
  
      end do ! edge loop
#ifdef USE_OMPOFFLOAD
      !$omp end distribute parallel do
      !$omp end target teams
#endif
   end do ! iteration loop

   call timerStop(timerOrig)
   call timerPrint(timerOrig)

   !--------------------------------------------------------------------
   ! Check for correctness and reset result
   !--------------------------------------------------------------------

   call timerStart(timerData)
#ifdef USE_OPENACC
   !$acc update host(highOrderFlx)
#endif
#ifdef USE_OMPOFFLOAD
   !$omp target update from(highOrderFlx)
#endif
   call timerStop(timerData)
   do iEdge=1,nEdges
   do k=1,nVertLevels
      refVal = refFlx(k,iEdge)
      relErr = abs(highOrderFlx(k,iEdge) - refVal)
      if (refVal /= 0.0_RKIND) relErr = relErr/abs(refVal)
      if (relErr > errTol) then
         print *,'Error computing highOrderFlx, loop 1: ', &
                  k,iEdge,highOrderFlx(k,iEdge),refVal
      endif

      highOrderFlx(k,iEdge) = 0.0_RKIND
   end do
   end do
   call timerStart(timerData)
#ifdef USE_OPENACC
   !$acc update device(highOrderFlx)
#endif
#ifdef USE_OMPOFFLOAD
   !$omp target update to(highOrderFlx)
#endif
   call timerStop(timerData)

   !--------------------------------------------------------------------
   ! Second loop form - parallel and GPU optimized
   !--------------------------------------------------------------------

   timerGPU  = timerCreate('GPU Optimized')
   call timerStart(timerGPU)

   do n =1,nIters
#ifdef USE_OPENACC
      !$acc parallel loop gang vector collapse(2) &
      !$acc    present(normalThicknessFlux, advMaskHighOrder, &
      !$acc            nAdvCellsForEdge, advCellsForEdge, cellMask, &
      !$acc            advCoefs, advCoefs3rd, tracerCur, highOrderFlx)
#endif
#ifdef USE_OMPOFFLOAD
      !$omp target teams &
      !$omp    map(to: normalThicknessFlux, advMaskHighOrder, &
      !$omp            nAdvCellsForEdge, advCellsForEdge, cellMask, &
      !$omp            advCoefs, advCoefs3rd, tracerCur) &
      !$omp    map(from: highOrderFlx)
      !$omp distribute parallel do collapse(2) &
      !$omp    private(iCell, coef1, coef2, coef3, edgeFlx,csgn)
#endif
      do iEdge = 1, nEdges
      do k = 1, nVertLevels
         ! Compute 3rd or 4th fluxes where requested.
         coef2 = normalThicknessFlux(k,iEdge)*advMaskHighOrder(k,iEdge)
         csgn = sign(1.0_RKIND,normalThicknessFlux(k,iEdge))
         edgeFlx = 0.0_RKIND
         do i = 1, nAdvCellsForEdge(iEdge)
            iCell = advCellsForEdge(i,iEdge)
            coef1 = advCoefs       (i,iEdge)
            coef3 = advCoefs3rd    (i,iEdge)*coef3rdOrder
            edgeFlx = edgeFlx + tracerCur(k,iCell)*cellMask(k,iCell)* &
                      coef2 * (coef1 + coef3*csgn)
         end do ! i loop over nAdvCellsForEdge
      
         highOrderFlx(k,iEdge) = edgeFlx
      end do ! vertical loop
      end do ! iEdge loop  
#ifdef USE_OMPOFFLOAD
      !$omp end distribute parallel do
      !$omp end target teams
#endif
   end do ! iteration loop

   call timerStop(timerGPU)
   call timerPrint(timerGPU)

   !--------------------------------------------------------------------
   ! Check for correctness and reset result
   !--------------------------------------------------------------------

   call timerStart(timerData)
#ifdef USE_OPENACC
   !$acc update host(highOrderFlx)
#endif
#ifdef USE_OMPOFFLOAD
   !$omp target update from(highOrderFlx)
#endif
   call timerStop(timerData)
   do iEdge=1,nEdges
   do k=1,nVertLevels
      refVal = refFlx(k,iEdge)
      relErr = abs(highOrderFlx(k,iEdge) - refVal)
      if (refVal /= 0.0_RKIND) relErr = relErr/abs(refVal)
      if (relErr > errTol) then
         print *,'Error computing highOrderFlx, loop 2: ', &
                  k,iEdge,highOrderFlx(k,iEdge),refVal
      endif

      highOrderFlx(k,iEdge) = 0.0_RKIND
   end do
   end do
   call timerStart(timerData)
#ifdef USE_OPENACC
   !$acc update device(highOrderFlx)
#endif
#ifdef USE_OMPOFFLOAD
   !$omp target update to(highOrderFlx)
#endif
   call timerStop(timerData)

   !--------------------------------------------------------------------
   ! Third loop form - parallel and GPU optimized with k tiling
   !--------------------------------------------------------------------

   timerGPU2  = timerCreate('GPU Optimized with k tiling')
   call timerStart(timerGPU2)

   do n =1,nIters
#ifdef USE_OPENACC
      !$acc parallel loop gang vector collapse(2) &
      !$acc    present(normalThicknessFlux, advMaskHighOrder, &
      !$acc            nAdvCellsForEdge, advCellsForEdge, cellMask, &
      !$acc            advCoefs, advCoefs3rd, tracerCur, highOrderFlx) &
      !$acc    private(coef2pk, csgnpk, edgeFlxPk)
#endif
#ifdef USE_OMPOFFLOAD
      !$omp target teams &
      !$omp    map(to: normalThicknessFlux, advMaskHighOrder, &
      !$omp            nAdvCellsForEdge, advCellsForEdge, cellMask, &
      !$omp            advCoefs, advCoefs3rd, tracerCur) &
      !$omp    map(from: highOrderFlx)
      !$omp distribute parallel do collapse(2) &
      !$omp    private(iCell, coef1, coef2pk, coef3, edgeFlxPk,csgnpk)
#endif
      do iEdge = 1, nEdges
      do k = 0, nvlpk-1
         k0 = packn*k
         ! Compute 3rd or 4th fluxes where requested.
         coef2pk = normalThicknessFlux(pkslc(k0),iEdge)*advMaskHighOrder(pkslc(k0),iEdge)
         csgnpk = sign(1.0_RKIND,normalThicknessFlux(pkslc(k0),iEdge))
         edgeFlxPk = 0.0_RKIND
         do i = 1, nAdvCellsForEdge(iEdge)
            iCell = advCellsForEdge(i,iEdge)
            coef1 = advCoefs       (i,iEdge)
            coef3 = advCoefs3rd    (i,iEdge)*coef3rdOrder
            edgeFlxPk = edgeFlxPk + tracerCur(pkslc(k0),iCell)*cellMask(pkslc(k0),iCell)* &
                 coef2pk*(coef1 + coef3*csgnpk)
         end do ! i loop over nAdvCellsForEdge
         highOrderFlx(pkslc(k0),iEdge) = edgeFlxPk
      end do ! vertical loop
      end do ! iEdge loop  
#ifdef USE_OMPOFFLOAD
      !$omp end distribute parallel do
      !$omp end target teams
#endif
   end do ! iteration loop

   call timerStop(timerGPU2)
   call timerPrint(timerGPU2)

   !--------------------------------------------------------------------
   ! Check for correctness and reset result
   !--------------------------------------------------------------------

   call timerStart(timerData)
#ifdef USE_OPENACC
   !$acc update host(highOrderFlx)
#endif
#ifdef USE_OMPOFFLOAD
   !$omp target update from(highOrderFlx)
#endif
   call timerStop(timerData)
   do iEdge=1,nEdges
   do k=1,nVertLevels
      refVal = refFlx(k,iEdge)
      relErr = abs(highOrderFlx(k,iEdge) - refVal)
      if (refVal /= 0.0_RKIND) relErr = relErr/abs(refVal)
      if (relErr > errTol) then
         print *,'Error computing highOrderFlx, loop GPU 2: ', &
                  k,iEdge,highOrderFlx(k,iEdge),refVal
      endif

      highOrderFlx(k,iEdge) = 0.0_RKIND
   end do
   end do
   call timerStart(timerData)
#ifdef USE_OPENACC
   !$acc update device(highOrderFlx)
#endif
#ifdef USE_OMPOFFLOAD
   !$omp target update to(highOrderFlx)
#endif
   call timerStop(timerData)
#endif

   !--------------------------------------------------------------------
   ! C++/Kokkos form
   !--------------------------------------------------------------------
#ifdef USE_CKE
   call kokkos_init()

   call timerStart(timerData)
   call cke_init(nIters, nEdges, nCells, nVertLevels, nvldim, nAdv, &
        nAdvCellsForEdge, minLevelCell, maxLevelCell, advCellsForEdge, &
        tracerCur, normalThicknessFlux, advMaskHighOrder, cellMask, &
        advCoefs, advCoefs3rd, coef3rdOrder)
   call timerStop(timerData)

   do i = 1,2
      timer_cke = timerCreate('C++/Kokkos ' // char(48+i))
      call timerStart(timer_cke)
      select case(i)
      case (1); call cke_impl1_run()
      case (2); call cke_impl2_run()
      ! other impls go here
      end select
      call timerStop(timer_cke)
      call timerPrint(timer_cke)

      call timerStart(timerData)
      call cke_get_results(nEdges, nVertLevels, highOrderFlx)
      call timerStop(timerData)

      iCell = 0;
      do iEdge=1,nEdges
         do k=1,nVertLevels
            refVal = refFlx(k,iEdge)
            relErr = abs(highOrderFlx(k,iEdge) - refVal)
            if (refVal /= 0.0_RKIND) relErr = relErr/abs(refVal)
            if ((isnan(relErr) .or. relErr > errTol) .and. iCell < 10) then
               print *,'Error computing highOrderFlx, C++/Kokkos impl', &
                    i,k,iEdge,highOrderFlx(k,iEdge),refVal
               iCell = iCell + 1
            endif
            highOrderFlx(k,iEdge) = 0.0_RKIND
         end do
      end do
   end do

   call cke_cleanup()
   call kokkos_finalize()
#endif

   !--------------------------------------------------------------------
   ! Clean up and deallocate
   !--------------------------------------------------------------------

   call timerStart(timerData)

#ifdef USE_OPENACC
   !$acc exit data delete(nAdvCellsForEdge, advCellsForEdge, &
   !$acc    minLevelCell, maxLevelCell, tracerCur, cellMask, &
   !$acc    normalThicknessFlux, highOrderFlx, advMaskHighOrder, &
   !$acc    advCoefs, advCoefs3rd, wgtTmp, sgnTmp)
#endif
#ifdef USE_OMPOFFLOAD
   !$omp target exit data map(delete: nAdvCellsForEdge, advCellsForEdge, &
   !$omp    minLevelCell, maxLevelCell, tracerCur, cellMask, &
   !$omp    normalThicknessFlux, highOrderFlx, advMaskHighOrder, &
   !$omp    advCoefs, advCoefs3rd, wgtTmp, sgnTmp)
#endif

   call timerStop(timerData)
   call timerPrint(timerData)

   deallocate(nAdvCellsForEdge, advCellsForEdge, minLevelCell, &
              maxLevelCell, tracerCur, cellMask, refFlx, &
              normalThicknessFlux, highOrderFlx, advMaskHighOrder, &
              advCoefs, advCoefs3rd, wgtTmp, sgnTmp)

#ifndef NO_MPI
   call MPI_Finalize(ierr)
#endif

!***********************************************************************

end program nested

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
