module cke_mod
  implicit none

  interface

     subroutine kokkos_init() bind(c)
     end subroutine kokkos_init

     subroutine kokkos_finalize() bind(c)
     end subroutine kokkos_finalize

     subroutine cke_init(nIters, nEdges, nCells, nVertLevels, nvldim, nAdv, &
          nAdvCellsForEdge, minLevelCell, maxLevelCell, advCellsForEdge, &
          tracerCur, normalThicknessFlux, advMaskHighOrder, cellMask, &
          advCoefs, advCoefs3rd, coef3rdOrder, highOrderFlx) bind(c)
       use iso_c_binding, only: c_int, c_double
       integer(c_int), value, intent(in) :: &
            nIters, nEdges, nCells, nVertLevels, nvldim, nAdv
       real(c_double), value, intent(in) :: &
            coef3rdOrder
       integer(c_int), intent(in) :: &
            nAdvCellsForEdge(nEdges), minLevelCell(nCells), maxLevelCell(nCells), &
            advCellsForEdge(nAdv,nEdges)
       real(c_double), dimension(nvldim,nCells), intent(in) :: &
            tracerCur, cellMask
       real(c_double), dimension(nvldim,nEdges), intent(in) :: &
            normalThicknessFlux, advMaskHighOrder
       real(c_double), dimension(nAdv,nEdges), intent(in) :: &
            advCoefs, advCoefs3rd
       real(c_double), dimension(nvldim,nEdges), intent(inout) :: &
            highOrderFlx
     end subroutine cke_init

     subroutine cke_impl1_run() bind(c)
     end subroutine cke_impl1_run
     subroutine cke_impl2_run() bind(c)
     end subroutine cke_impl2_run

     subroutine cke_get_results(nEdges, nvldim, highOrderFlx) bind(c)
       use iso_c_binding, only: c_int, c_double
       integer(c_int), value, intent(in) :: &
            nEdges, nvldim
       real(c_double), dimension(nvldim,nEdges), intent(out) :: &
            highOrderFlx
     end subroutine cke_get_results

     subroutine cke_cleanup() bind(c)
     end subroutine cke_cleanup

  end interface

end module cke_mod
