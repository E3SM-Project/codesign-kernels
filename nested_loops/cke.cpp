#include "cke.hpp"
#include "cke_impl.hpp"

static bool in_charge_of_kokkos = false;

void kokkos_init () {
  in_charge_of_kokkos = ! Kokkos::is_initialized();
  ekat::initialize_ekat_session(false);
}

void kokkos_finalize () {
  if (in_charge_of_kokkos && Kokkos::is_initialized())
    Kokkos::finalize();
}

namespace cke {

#ifdef CKE_SAME_MEMORY
static constexpr bool layout_right = std::is_same<Data::Layout, Kokkos::LayoutRight>::value;
#endif

// Initialize a View<const Pack<Real,packn>**> from raw(1:d1,1:d2), where dim 2
// has the fast index.
template <typename Scalar, typename V> static void
initvpk (const Scalar* raw, const int d1, const int d2, const std::string& name, V& v,
         typename std::enable_if<V::value_type::packtag>::type* = 0) {
  // If not on the GPU and other conditions are met, we can use the same memory
  // in F90 and C++.
#ifdef CKE_SAME_MEMORY
  if ( ! ekat::OnGpu<typename V::execution_space>::value &&
       layout_right && d2 % V::value_type::n == 0) {
    v = V(reinterpret_cast<const typename V::value_type*>(raw), d1, d2/V::value_type::n);
    return;
  }
#endif
  // Get the number of packs that cover the scalar length.
  const int d2pk = ekat::PackInfo<V::value_type::n>::num_packs(d2);
  // Allocate the view as writeable.
  const auto vnc = typename V::non_const_type(name, d1, d2pk);
  // For convenience, take a scalar view of the original Pack<Scalar> view.
  const auto svnc = scalarize(vnc);
  // Copy the F90 data to a host mirror of the scalar view.
  const auto h = Kokkos::create_mirror_view(svnc);
  for (int i = 0; i < d1; ++i)
    for (int j = 0; j < d2; ++j)
      h(i,j) = raw[d2*i+j];
  // Copy the data to device.
  Kokkos::deep_copy(svnc, h);
  // Set the possibly read-only view with this data.
  v = vnc;
}

// Simpler forms of the above: scalar 1D and 2D views.

template <typename Scalar, typename V> static
void initv (const Scalar* raw, const int d1, const std::string& name, V& v,
            const Scalar delta = 0) {
#ifdef CKE_SAME_MEMORY
  if ( ! ekat::OnGpu<typename V::execution_space>::value && delta == 0) {
    v = V(raw, d1);
    return;
  }
#endif
  const auto vnc = typename V::non_const_type(name, d1);
  const auto h = Kokkos::create_mirror_view(vnc);
  for (int i = 0; i < d1; ++i) h(i) = raw[i] + delta;
  Kokkos::deep_copy(vnc, h);
  v = vnc;
}

template <typename Scalar, typename V> static
void initv (const Scalar* raw, const int d1, const int d2, const std::string& name,
            V& v, const Scalar delta = 0) {
#ifdef CKE_SAME_MEMORY
  if ( ! ekat::OnGpu<typename V::execution_space>::value &&
       layout_right && delta == 0) {
    v = V(raw, d1, d2);
    return;
  }
#endif
  const auto vnc = typename V::non_const_type(name, d1, d2);
  const auto h = Kokkos::create_mirror_view(vnc);
  for (int i = 0; i < d1; ++i)
    for (int j = 0; j < d2; ++j)
      h(i,j) = raw[d2*i+j] + delta;
  Kokkos::deep_copy(vnc, h);
  v = vnc;
}

void Data::init (
  const Int nIters_, const Int nEdges_, const Int nCells_, const Int nVertLevels_,
  const Int nvldim_, const Int nAdv_, const Int* nAdvCellsForEdge_, const Int* minLevelCell_,
  const Int* maxLevelCell_, const Int* advCellsForEdge_, Real* tracerCur_,
  const Real* normalThicknessFlux_, const Real* advMaskHighOrder_, const Real* cellMask_,
  const Real* advCoefs_, const Real* advCoefs3rd_, const Real coef3rdOrder_,
  Real* highOrderFlx_)
{
  nIters = nIters_; nEdges = nEdges_; nCells = nCells_; nVertLevels = nVertLevels_;
  nvldim = nvldim_; nAdv = nAdv_; coef3rdOrder = coef3rdOrder_;
  nvlpk = ekat::PackInfo<Data::packn>::num_packs(nVertLevels);

  initv(nAdvCellsForEdge_, nEdges, "nAdvCellsForEdge", nAdvCellsForEdge);
  initv(minLevelCell_, nCells, "minLevelCell", minLevelCell, -1);
  initv(maxLevelCell_, nCells, "maxLevelCell", maxLevelCell, -1);
  initv(advCellsForEdge_, nEdges, nAdv, "advCellsForEdge", advCellsForEdge, -1);
  initv(advCoefs_, nEdges, nAdv, "advCoefs", advCoefs);
  initv(advCoefs3rd_, nEdges, nAdv, "advCoefs3rd", advCoefs3rd);
  initvpk(tracerCur_, nCells, nvldim, "tracerCur", tracerCur);
  initvpk(cellMask_, nCells, nvldim, "cellMask", cellMask);
  initvpk(normalThicknessFlux_, nEdges, nvldim, "normalThicknessFlux", normalThicknessFlux);
  initvpk(advMaskHighOrder_, nEdges, nvldim, "advMaskHighOrder", advMaskHighOrder);

  const int npack = ekat::PackInfo<packn>::num_packs(nVertLevels);
#ifdef CKE_SAME_MEMORY
  if ( ! ekat::OnGpu<ExeSpace>::value && layout_right && nvldim % Pr::n == 0)
    highOrderFlx = Apr2(reinterpret_cast<Pr*>(highOrderFlx_), nEdges, npack);
  else
#endif
    highOrderFlx = Apr2("highOrderFlx", nEdges, npack);
}

static Data::Ptr g_data;

Data::Ptr get_Data_singleton () { return g_data; }

} // namespace cke

void cke_init (
  const Int nIters, const Int nEdges, const Int nCells, const Int nVertLevels,
  const Int nvldim, const Int nAdv, const Int* nAdvCellsForEdge, const Int* minLevelCell,
  const Int* maxLevelCell, const Int* advCellsForEdge, Real* tracerCur,
  const Real* normalThicknessFlux, const Real* advMaskHighOrder, const Real* cellMask,
  const Real* advCoefs, const Real* advCoefs3rd, const Real coef3rdOrder,
  Real* highOrderFlx)
{
  cke::g_data = std::make_shared<cke::Data>();
  cke::g_data->init(nIters, nEdges, nCells, nVertLevels, nvldim, nAdv,
                    nAdvCellsForEdge, minLevelCell, maxLevelCell, advCellsForEdge,
                    tracerCur, normalThicknessFlux, advMaskHighOrder, cellMask,
                    advCoefs, advCoefs3rd, coef3rdOrder, highOrderFlx);
}

void cke_get_results (const Int, const Int, Real* highOrderFlx) {
  const auto d = cke::g_data;
  assert(d);
  const auto nvldim = d->nvldim;
#ifdef CKE_SAME_MEMORY
  if ( ! ekat::OnGpu<cke::Data::ExeSpace>::value && cke::layout_right &&
       nvldim % cke::Data::Pr::n == 0)
    return;
#endif
  const auto shof = scalarize(d->highOrderFlx);
  const auto h = Kokkos::create_mirror_view(shof);
  Kokkos::deep_copy(h, shof);
  for (int i = 0; i < d->nEdges; ++i)
    for (int j = 0; j < nvldim; ++j)
      highOrderFlx[nvldim*i+j] = h(i,j);
  Kokkos::deep_copy(d->highOrderFlx, 0);
}

void cke_cleanup () { cke::g_data = nullptr; }
