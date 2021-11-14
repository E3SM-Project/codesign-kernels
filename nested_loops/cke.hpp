#ifndef INCLUDE_CKE_HPP
#define INCLUDE_CKE_HPP

typedef int Int;
typedef double Real;

extern "C" {
  void kokkos_init();
  void kokkos_finalize();

  void cke_init(
    const Int nIters, const Int nEdges, const Int nCells, const Int nVertLevels,
    const Int nvldim, const Int nAdv, const Int* nAdvCellsForEdge, const Int* minLevelCell,
    const Int* maxLevelCell, const Int* advCellsForEdge, const Real* tracerCur,
    const Real* normalThicknessFlux, const Real* advMaskHighOrder, const Real* cellMask,
    const Real* advCoefs, const Real* advCoefs3rd, const Real coef3rdOrder,
    Real* highOrderFlx);
  void cke_get_results(const Int nEdges, const Int nVertLevels, Real* highOrderFlx);
  void cke_cleanup();

  void cke_impl1_run();
  void cke_impl2_run();
}

#endif
