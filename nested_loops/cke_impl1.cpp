#include "cke.hpp"
#include "cke_impl.hpp"

//sec C++/Kokkos/EKAT demo implementation 1.

namespace cke {

struct NestedLoopKernel {
  struct Data d;
  NestedLoopKernel (const Data& d_) : d(d_) {}
  KOKKOS_FORCEINLINE_FUNCTION void operator() (const Int iEdge, const Int k) const {
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
      edgeFlx += d.tracerCur(iCell,k)*d.cellMask(iCell,k)*coef2*(coef1 + coef3*csgn);
    }
    d.highOrderFlx(iEdge,k) = edgeFlx;
  }
};

void run (const Data& d) {
  NestedLoopKernel nlk(d);
  parfor_iEdge_kPack(d, nlk);
  Kokkos::fence();
}

} // namespace cke

void cke_impl1_run () {
  const auto d = cke::get_Data_singleton();
  assert(d);
  for (int iter = 0; iter < d->nIters; ++iter)
    cke::run(*d);
}
