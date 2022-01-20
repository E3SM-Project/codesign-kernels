!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
! nested_vars
!
!***********************************************************************

#define C_WRAP(var) c_##var = c_loc(var)

module nested_vars

   use timerMod
   use iso_c_binding
   implicit none

   public
   
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
 
   type(c_ptr), BIND(C, NAME='c_maxLevelCell')        :: c_maxLevelCell
   type(c_ptr), BIND(C, NAME='c_minLevelCell')        :: c_minLevelCell
   type(c_ptr), BIND(C, NAME='c_nAdvCellsForEdge')    :: c_nAdvCellsForEdge
   type(c_ptr), BIND(C, NAME='c_advCellsForEdge')     :: c_advCellsForEdge
   type(c_ptr), BIND(C, NAME='c_tracerCur')           :: c_tracerCur
   type(c_ptr), BIND(C, NAME='c_normalThicknessFlux') :: c_normalThicknessFlux
   type(c_ptr), BIND(C, NAME='c_advMaskHighOrder')    :: c_advMaskHighOrder
   type(c_ptr), BIND(C, NAME='c_cellMask')            :: c_cellMask
   type(c_ptr), BIND(C, NAME='c_highOrderFlx')        :: c_highOrderFlx
   type(c_ptr), BIND(C, NAME='c_refFlx')              :: c_refFlx
   type(c_ptr), BIND(C, NAME='c_advCoefs')            :: c_advCoefs
   type(c_ptr), BIND(C, NAME='c_advCoefs3rd')         :: c_advCoefs3rd

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

contains

   subroutine alloc_vars
   
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

#ifdef USE_YAKL
      C_WRAP(minLevelCell)
      C_WRAP(maxLevelCell)
      C_WRAP(nAdvCellsForEdge)
      C_WRAP(advCellsForEdge)
      C_WRAP(tracerCur)
      C_WRAP(normalThicknessFlux)
      C_WRAP(advMaskHighOrder)
      C_WRAP(cellMask)
      C_WRAP(highOrderFlx)
      C_WRAP(refFlx)
      C_WRAP(advCoefs)
      C_WRAP(advCoefs3rd)
#endif

   end subroutine alloc_vars
   
!***********************************************************************
end module nested_vars
