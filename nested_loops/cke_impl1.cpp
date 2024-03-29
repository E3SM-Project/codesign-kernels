#include "cke.hpp"
#include "cke_impl.hpp"

using namespace cke;

// C++/Kokkos/EKAT demo implementation 1.

namespace {

void run (const Data& d) {
  const auto f1 = KOKKOS_LAMBDA(const int idx) {
    int iCell, k;
    d.get_iCell_kPack_idxs(idx, iCell, k);
    d.tracerCur(iCell,k) *= d.cellMask(iCell,k);
  };
  Kokkos::parallel_for(d.get_rpolicy_iCell_kPack(), f1);
  Kokkos::fence();
  const auto f2 = KOKKOS_LAMBDA(const int idx) {
    int iEdge, k;
    d.get_iEdge_kPack_idxs(idx, iEdge, k);
    const auto coef2 = d.normalThicknessFlux(iEdge,k)*d.advMaskHighOrder(iEdge,k);
    Data::Pr csgn; {
      const auto ntf = d.normalThicknessFlux(iEdge,k);
      vector_simd for (int s = 0; s < Data::packn; ++s)
        csgn[s] = ntf[s] < 0 ? -1 : 1;
    }
    Data::Pr edgeFlx(0);
    const auto iend = d.nAdvCellsForEdge(iEdge);
    for (int i = 0; i < iend; ++i) {
      const auto coef1 = d.advCoefs(iEdge,i);
      const auto coef3 = d.advCoefs3rd(iEdge,i)*d.coef3rdOrder;
      const auto iCell = d.advCellsForEdge(iEdge,i);
      edgeFlx += d.tracerCur(iCell,k)*coef2*(coef1 + coef3*csgn);
    }
    d.highOrderFlx(iEdge,k) = edgeFlx;
  };
  Kokkos::parallel_for(d.get_rpolicy_iEdge_kPack(), f2);
  Kokkos::fence();
}

} // namespace

void cke_impl1_run () {
  const auto d = cke::get_Data_singleton();
  assert(d);
  for (int iter = 0; iter < d->nIters; ++iter)
    run(*d);
}
