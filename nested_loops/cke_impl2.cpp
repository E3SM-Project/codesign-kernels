#include "cke.hpp"
#include "cke_impl.hpp"

// C++/Kokkos/EKAT demo implementation 2. This code demonstrates hierarchical
// parallelism, including suing shared scratch memory. It isn't very fast
// because for this particular kernel, hierarchical ||ism is not well motivated.
// A little bit of work might speed it up a little, though.

using namespace cke;

namespace {

struct NestedLoopKernel {
  struct Data d;
  NestedLoopKernel (const Data& d_) : d(d_) {}

  KOKKOS_INLINE_FUNCTION size_t get_col_size () const {
    return d.nvlpk*sizeof(Data::Pr);
  }

  size_t team_shmem_size (const int team_size) const {
    return 3*get_col_size();
  }

  KOKKOS_FORCEINLINE_FUNCTION
  void operator() (const Data::TeamMember& team) const {
    const int iEdge = team.league_rank();
    const auto col_sz = get_col_size();
    auto* const wgtTmp = static_cast<Data::Pr*>(team.team_shmem().get_shmem(3*col_sz));
    auto* const sgnTmp = wgtTmp + d.nvlpk;
    auto* const flxTmp = sgnTmp + d.nvlpk;
    const auto f1 = [&] (const int k) {
      wgtTmp[k] = d.normalThicknessFlux(iEdge,k)*d.advMaskHighOrder(iEdge,k);
      Data::Pr csgn; {
        const auto ntf = d.normalThicknessFlux(iEdge,k);
        vector_simd for (int s = 0; s < Data::packn; ++s)
          csgn[s] = ntf[s] < 0 ? -1 : 1;
      }
      sgnTmp[k] = csgn;
      flxTmp[k] = 0;
    };
    // Technically need a vector-range loop or single inside the team thread
    // loop, but not needed if the vector extent is guaranteed to be 1.
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, d.nvlpk), f1);
    team.team_barrier();
    const auto iend = d.nAdvCellsForEdge(iEdge);
    for (int i = 0; i < iend; ++i) {
      const auto coef1 = d.advCoefs(iEdge,i);
      const auto coef3 = d.advCoefs3rd(iEdge,i)*d.coef3rdOrder;
      const auto iCell = d.advCellsForEdge(iEdge,i);
      const int kbeg = d.minLevelCell(iCell)/Data::packn;
      const int kend = d.maxLevelCell(iCell)/Data::packn;
      const auto f2 = [&] (const int kos) {
        const int k = kbeg + kos;
        flxTmp[k] += (d.tracerCur(iCell,k)*d.cellMask(iCell,k)*
                      wgtTmp[k]*(coef1 + coef3*sgnTmp[k]));
      };
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, kend-kbeg+1), f2);
    }
    team.team_barrier();
    const auto f3 = [&] (const int k) {
      d.highOrderFlx(iEdge,k) = flxTmp[k];
    };
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, d.nvlpk), f3);
  }
};

void run (const Data& d) {
  NestedLoopKernel nlk(d);
  Kokkos::parallel_for(d.get_tpolicy_iEdge(), nlk);
  Kokkos::fence();
}

} // namespace

void cke_impl2_run () {
  const auto d = cke::get_Data_singleton();
  assert(d);
  for (int iter = 0; iter < d->nIters; ++iter)
    run(*d);
}
