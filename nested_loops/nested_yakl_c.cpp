#include <iostream>
#include <YAKL.h>
#include <fstream>
#include <random>

#define DRAND   (rand() / rdenom)
#define SIGN(a,b)   std::abs(a) * ( ( (b) < 0 ) ? -1 : 1)
#define YAKL_LOCAL_NS(ns, var) auto & var = *(ns::var)
#define YAKL_LOCAL_VAR(var) auto & var = ::var

typedef yakl::Array<int,1,yakl::memHost,yakl::styleFortran>     h_int_1d_t;
typedef yakl::Array<int,2,yakl::memHost,yakl::styleFortran>     h_int_2d_t;
typedef yakl::Array<double,2,yakl::memHost,yakl::styleFortran>  h_double_2d_t;

typedef yakl::Array<int,1,yakl::memDefault,yakl::styleFortran>     d_int_1d_t;
typedef yakl::Array<int,2,yakl::memDefault,yakl::styleFortran>     d_int_2d_t;
typedef yakl::Array<double,1,yakl::memDefault,yakl::styleFortran>  d_double_1d_t;
typedef yakl::Array<double,2,yakl::memDefault,yakl::styleFortran>  d_double_2d_t;


template<typename R, int Mem, int Dims>
std::vector<int> getBounds(yakl::Array<R,Dims,Mem,yakl::styleFortran> & arr)
{
    std::vector<int>       bnds;

    auto rnk = arr.get_rank();
    for ( int i = 0; i < rnk; ++i )
    {
        bnds.push_back(arr.extent(i));
    }
    return bnds;
}

template <typename R, typename...T>
yakl::Array<R,sizeof...(T),yakl::memDefault,yakl::styleFortran> *
yakl_create_array(T...dims)
{
    typedef yakl::Array<R,sizeof...(T),yakl::memDefault,yakl::styleFortran> ret_type;

    ret_type * w_var_p = new ret_type("var_p", dims...);

    return w_var_p;
}

template <typename R, typename...T>
yakl::Array<R,sizeof...(T),yakl::memDefault,yakl::styleFortran> *
yakl_wrap_array(R * var_p, T...dims)
{
    typedef yakl::Array<R,sizeof...(T),yakl::memHost,yakl::styleFortran> arr_type;
    typedef yakl::Array<R,sizeof...(T),yakl::memDefault,yakl::styleFortran> ret_type;

    arr_type h_var("h_var_p", var_p, dims...);
    ret_type * w_var_p = new ret_type("var_p", dims...);

    h_var.deep_copy_to(*w_var_p);

    return w_var_p;
}


template <typename R, int N>
void
yakl_update_host(yakl::Array<R,N,yakl::memDefault,yakl::styleFortran> * d_var,
                 R * h_var_p)
{
    typedef yakl::Array<R,N,yakl::memHost,yakl::styleFortran>  host_type;

    auto & rw_var = *d_var;
    host_type   h_var("h_var_p", h_var_p, getBounds(rw_var));
    rw_var.deep_copy_to(h_var);
}

template <typename R, int N>
void
yakl_update_device(yakl::Array<R,N,yakl::memDefault,yakl::styleFortran> * d_var,
                 R * h_var_p)
{
    typedef yakl::Array<R,N,yakl::memHost,yakl::styleFortran>  host_type;

    auto & rw_var = *d_var;
    host_type   h_var("h_var_p", h_var_p, getBounds(rw_var));
    h_var.deep_copy_to(rw_var);
}

namespace
{
int   nvldim,       
      nEdges,       // number of edges
      nCells,       // number of cells
      nVertLevels,  // number of vertical levels
      nAdv          // max num cells contributing to edge for advection
      ;
      
double
      coef3rdOrder = 2.14; // an advection parameter

}

namespace nested
{
d_int_1d_t      * maxLevelCell,
                * minLevelCell,
                * nAdvCellsForEdge
                ;

d_int_2d_t      * advCellsForEdge
                ;

d_double_1d_t   * areaCell,
                * inverseAreaCell,
                * dcEdge,
                * dvEdge
                ;
                
d_double_2d_t   
                * wgtTmp,
                * sgnTmp,
                * tracerCur,             // current tracer values
                * normalThicknessFlux,   // thickness flux normal to edge
                * advMaskHighOrder,      // mask for using high-order advection
                * cellMask,              // mask for active cells
                * highOrderFlx,          // temp for holding high-order adv flux
                * refFlx,                // reference value for high-order adv flux
                * advCoefs,              // precomputed adv coefficients
                * advCoefs3rd            // precomputed adv coefficients 3rd-order
                ;
}


extern "C" int * c_maxLevelCell = nullptr;
extern "C" int * c_minLevelCell = nullptr;
extern "C" int * c_nAdvCellsForEdge = nullptr;
extern "C" int * c_advCellsForEdge = nullptr;
extern "C" double * c_advCoefs = nullptr;
extern "C" double * c_advCoefs3rd = nullptr;
extern "C" double * c_cellMask = nullptr;
extern "C" double * c_advMaskHighOrder = nullptr;
extern "C" double * c_tracerCur = nullptr;
extern "C" double * c_normalThicknessFlux = nullptr;
extern "C" double * c_highOrderFlx = nullptr;

extern "C"
void yakl_init()
{
    yakl::init();
}

extern "C"
void yakl_init_arrays( int fnvldim,
                int fnEdges,
                int fnCells,
                int fnVertLevels,
                int fnAdv,
                double fcoef3rdOrder)
{
    using namespace nested;

    nvldim = fnvldim;
    nEdges = fnEdges;
    nCells = fnCells;
    nVertLevels = fnVertLevels;
    nAdv = fnAdv;
    coef3rdOrder = fcoef3rdOrder;
    
    maxLevelCell = yakl_wrap_array(c_maxLevelCell, nCells);
    minLevelCell = yakl_wrap_array(c_minLevelCell, nCells);
    
    nAdvCellsForEdge = yakl_wrap_array(c_nAdvCellsForEdge, nEdges);
    advCellsForEdge = yakl_wrap_array(c_advCellsForEdge, nAdv, nEdges);
    advCoefs = yakl_wrap_array(c_advCoefs, nAdv, nEdges);
    advCoefs3rd = yakl_wrap_array(c_advCoefs3rd, nAdv, nEdges);
    cellMask = yakl_wrap_array(c_cellMask, nvldim, nCells);
    advMaskHighOrder = yakl_wrap_array(c_advMaskHighOrder, nvldim, nEdges);
    
    tracerCur = yakl_create_array<double>(nvldim,nCells);
    normalThicknessFlux = yakl_create_array<double>(nvldim,nEdges);
    highOrderFlx = yakl_create_array<double>(nvldim,nEdges);
    wgtTmp = yakl_create_array<double>(nVertLevels,nEdges);
    sgnTmp = yakl_create_array<double>(nVertLevels,nEdges);
}


extern "C"
void yakl_original_loop(double * h_highOrderFlx, double * h_tracerCur,
                        double * h_normalThicknessFlux)
{
    yakl_update_device(nested::tracerCur, h_tracerCur);
    yakl_update_device(nested::normalThicknessFlux, h_normalThicknessFlux);

    YAKL_LOCAL_NS(nested, wgtTmp);
    YAKL_LOCAL_NS(nested, sgnTmp);
    YAKL_LOCAL_NS(nested, minLevelCell);
    YAKL_LOCAL_NS(nested, maxLevelCell);
    YAKL_LOCAL_NS(nested, tracerCur);
    YAKL_LOCAL_NS(nested, highOrderFlx);
    YAKL_LOCAL_NS(nested, normalThicknessFlux);
    YAKL_LOCAL_NS(nested, advMaskHighOrder);
    YAKL_LOCAL_NS(nested, nAdvCellsForEdge);
    YAKL_LOCAL_NS(nested, advCellsForEdge);
    YAKL_LOCAL_NS(nested, advCoefs);
    YAKL_LOCAL_NS(nested, advCoefs3rd);
    YAKL_LOCAL_VAR(nVertLevels);
    YAKL_LOCAL_VAR(coef3rdOrder);

    yakl::fortran::parallel_for( yakl::fortran::Bounds<1>({1,nEdges}) ,
    YAKL_LAMBDA(int iEdge)
    {
         // compute some common intermediate factors
         for (int k = 1; k <= nVertLevels; ++k )
         {
            wgtTmp(k,iEdge) = normalThicknessFlux(k,iEdge)*
                        advMaskHighOrder(k,iEdge);
            sgnTmp(k,iEdge) = SIGN(1.0,
                             normalThicknessFlux(k,iEdge));
            highOrderFlx(k,iEdge) = 0.0;
         }

         // Compute 3rd or 4th fluxes where requested.
         for ( int i = 1; i <= nAdvCellsForEdge(iEdge); ++i )
         {
            int iCell = advCellsForEdge(i,iEdge);
            int kmin  = minLevelCell(iCell);
            int kmax  = maxLevelCell(iCell);
            double coef1 = advCoefs       (i,iEdge);
            double coef3 = advCoefs3rd    (i,iEdge)*coef3rdOrder;
            for ( int k = kmin; k <= kmax; ++k )
            {
               highOrderFlx(k,iEdge) = highOrderFlx(k,iEdge) + tracerCur(k,iCell)*
                           wgtTmp(k,iEdge)*(coef1 + coef3*sgnTmp(k,iEdge));
            }
         }
    });
    
    yakl_update_host(nested::highOrderFlx, h_highOrderFlx);
}

extern "C" 
void yakl_check(double * h_cellMask, double * h_tracerCur)
{
    yakl_update_device(nested::tracerCur, h_tracerCur);

    YAKL_LOCAL_NS(nested, cellMask);
    YAKL_LOCAL_NS(nested, tracerCur);

    h_double_2d_t   hcmask("hcmask", h_cellMask, nvldim, nCells);
    h_double_2d_t   tracer("tracer", h_tracerCur, nvldim, nCells);
    for ( int n = 1; n <= nCells; ++n )
    {
        for ( int k = 1; k <= nVertLevels; ++k )
        {
            if ( hcmask(k,n) != cellMask(k,n) )
            {
                fprintf(stderr, "cellMask: Expected %lf, but got %lf, for k,n = %d, %d\n",
                    hcmask(k,n), cellMask(k,n), k, n);
            }

            if ( tracer(k,n) != tracerCur(k,n) )
            {
                fprintf(stderr, "tracerCur: Expected %lf, but got %lf, for k,n = %d, %d\n",
                    tracer(k,n), tracerCur(k,n), k, n);
            }
        }
    }
}

extern "C"
void yakl_gpu_optimized(double * h_highOrderFlx, double * h_tracerCur,
                        double * h_normalThicknessFlux)
{
    yakl_update_device(nested::tracerCur, h_tracerCur);
    yakl_update_device(nested::normalThicknessFlux, h_normalThicknessFlux);

    YAKL_LOCAL_NS(nested, tracerCur);
    YAKL_LOCAL_NS(nested, highOrderFlx);
    YAKL_LOCAL_NS(nested, normalThicknessFlux);
    YAKL_LOCAL_NS(nested, advMaskHighOrder);
    YAKL_LOCAL_NS(nested, nAdvCellsForEdge);
    YAKL_LOCAL_NS(nested, advCellsForEdge);
    YAKL_LOCAL_NS(nested, advCoefs);
    YAKL_LOCAL_NS(nested, advCoefs3rd);
    YAKL_LOCAL_NS(nested, cellMask);
    YAKL_LOCAL_VAR(nVertLevels);
    YAKL_LOCAL_VAR(coef3rdOrder);

    
    yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nEdges},{1,nVertLevels}) ,
    YAKL_LAMBDA(int iEdge, int k)
    {
         // Compute 3rd or 4th fluxes where requested.
         double coef2 = normalThicknessFlux(k,iEdge)*advMaskHighOrder(k,iEdge);
         double csgn = SIGN(1.0,
                             normalThicknessFlux(k,iEdge));
         double edgeFlx = 0.0;
         for ( int i = 1; i <= nAdvCellsForEdge(iEdge); ++i )
         {
            int iCell = advCellsForEdge(i,iEdge);
            double coef1 = advCoefs       (i,iEdge);
            double coef3 = advCoefs3rd    (i,iEdge)*coef3rdOrder;
            edgeFlx = edgeFlx + tracerCur(k,iCell)*cellMask(k,iCell)*
                      coef2 * (coef1 + coef3*csgn);
         }
         highOrderFlx(k,iEdge) = edgeFlx;
    });

    yakl_update_host(nested::highOrderFlx, h_highOrderFlx);
}
