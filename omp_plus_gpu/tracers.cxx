#include "YAKL.h"
#include <vector>
#include <omp.h>
#include <random>
#include <chrono>
#include <ctime>
#include <cmath>

#define YAKL_LOCAL(var) auto & var = *(mpas::var)
#define YAKL_LOCAL_TID(var, tid) auto & var = *(mpas::var[tid].field)
#define YAKL_UPD_LOCAL_TID(var, tid) var = *(mpas::var[tid].field)
#define SIGN(a,b) ( (b) >= 0 ) ? std::abs(a) : -std::abs(a)
#define DMIN(a,b) std::min<double>(a,b)
#define DMAX(a,b) std::max<double>(a,b)

typedef double real;

typedef yakl::Array<real,1,yakl::memHost,yakl::styleFortran>  h_double_1d_t;
typedef yakl::Array<real,1,yakl::memDefault,yakl::styleFortran>  d_double_1d_t;
typedef yakl::Array<int,1,yakl::memHost,yakl::styleFortran>  h_int_1d_t;
typedef yakl::Array<int,1,yakl::memDefault,yakl::styleFortran>  d_int_1d_t;

typedef yakl::Array<real,2,yakl::memHost,yakl::styleFortran>  h_double_2d_t;
typedef yakl::Array<real,2,yakl::memDefault,yakl::styleFortran>  d_double_2d_t;
typedef yakl::Array<int,2,yakl::memDefault,yakl::styleFortran>  d_int_2d_t;
typedef yakl::Array<int,2,yakl::memHost,yakl::styleFortran>  h_int_2d_t;

typedef yakl::Array<real,3,yakl::memHost,yakl::styleFortran>  h_double_3d_t;
typedef yakl::Array<real,3,yakl::memDefault,yakl::styleFortran>  d_double_3d_t;

typedef struct {
    d_double_2d_t  * field;
} field2d_t;


std::vector<yakl::Stream>   streams;

int ntracers = 42;
int nthreads = 1;
int nVertLevels = 100;
//int nEdges = 25600;
//int nCells = 2800;
int nEdges = 201600;
int nCells = 22400;

// taken from https://gist.github.com/mcleary/b0bf4fa88830ff7c882d
class Timer
{
public:
    void start()
    {
        m_StartTime = std::chrono::system_clock::now();
        m_bRunning = true;
    }
    
    void stop()
    {
        m_EndTime = std::chrono::system_clock::now();
        m_bRunning = false;
    }
    
    double elapsedMilliseconds()
    {
        std::chrono::time_point<std::chrono::system_clock> endTime;
        
        if(m_bRunning)
        {
            endTime = std::chrono::system_clock::now();
        }
        else
        {
            endTime = m_EndTime;
        }
        
        return std::chrono::duration_cast<std::chrono::milliseconds>(endTime - m_StartTime).count();
    }
    
    double elapsedSeconds()
    {
        return elapsedMilliseconds() / 1000.0;
    }

private:
    std::chrono::time_point<std::chrono::system_clock> m_StartTime;
    std::chrono::time_point<std::chrono::system_clock> m_EndTime;
    bool                                               m_bRunning = false;
};

namespace mpas {

int nAdv = 10;

double          coef3rdOrder = 0.25;
double          eps = 1.e-10;

d_int_1d_t
   * nAdvCellsForEdge,
   * minLevelCell,
   * maxLevelCell,
   * minLevelEdgeBot,
   * maxLevelEdgeTop
   ;

d_int_2d_t
   * advCellsForEdge,
   * cellsOnEdge
   ;

d_double_1d_t  
   * dvEdge
   ;

d_double_2d_t  
   * normalThicknessFlux,
   * advMaskHighOrder,
   * advCoefs,
   * advCoefs3rd
   ;

std::vector<field2d_t>   tracers;
std::vector<field2d_t>   wgtTmp;
std::vector<field2d_t>   sgnTmp;
std::vector<field2d_t>   tracerMin;
std::vector<field2d_t>   tracerMax;
std::vector<field2d_t>   highOrderFlx;
std::vector<field2d_t>   lowOrderFlx;
};

template <typename R, typename...T>
yakl::Array<R,sizeof...(T),yakl::memDefault,yakl::styleFortran> *
yakl_create_array(const std::string & name, T...dims)
{
    typedef yakl::Array<R,sizeof...(T),yakl::memDefault,yakl::styleFortran> ret_type;

    ret_type * w_var_p = new ret_type(name.c_str(), dims...);

    return w_var_p;
}

template <typename...T>
yakl::Array<real,sizeof...(T),yakl::memDefault,yakl::styleFortran> *
yakl_create_real(const std::string & name, T...Dims)
{
    return yakl_create_array<real>(name, Dims...);
}

template <typename...T>
yakl::Array<int,sizeof...(T),yakl::memDefault,yakl::styleFortran> *
yakl_create_int(const std::string & name, T...Dims)
{
    return yakl_create_array<int>(name, Dims...);
}

template <typename R, int N>
void
yakl_update_device(yakl::Array<R,N,yakl::memDefault,yakl::styleFortran> * d_var,
                 R * h_var_p,yakl::Stream stream = yakl::Stream() )
{
    typedef yakl::Array<R,N,yakl::memHost,yakl::styleFortran>  host_type;

    auto & rw_var = *d_var;
    host_type   h_var("h_var_p", h_var_p, getBounds(rw_var));
    h_var.deep_copy_to(rw_var, stream);
}

template <typename R, int N>
void
yakl_update_device(yakl::Array<R,N,yakl::memDefault,yakl::styleFortran> * d_var,
                 yakl::Array<R,N,yakl::memHost,yakl::styleFortran> * h_var_p,yakl::Stream stream = yakl::Stream() )
{
    h_var_p->deep_copy_to(d_var, stream);
}

void alloc_fields() {

   using namespace mpas;

   int maxNAdvCells = 1;

   nthreads = omp_get_max_threads();

   tracers.resize(ntracers);
   tracerMin.resize(nthreads);
   tracerMax.resize(nthreads);
   wgtTmp.resize(nthreads);
   sgnTmp.resize(nthreads);
   highOrderFlx.resize(nthreads);
   lowOrderFlx.resize(nthreads);

   cellsOnEdge = yakl_create_int("cellsOnEdge", 2, nEdges);
   nAdvCellsForEdge = yakl_create_int("nAdvCellsForEdge", nEdges);
   minLevelCell = yakl_create_int("minLevelCell", nCells);
   maxLevelCell = yakl_create_int("maxLevelCell", nCells);
   minLevelEdgeBot = yakl_create_int("minLevelEdgeBot", nEdges);
   maxLevelEdgeTop = yakl_create_int("maxLevelEdgeTop", nEdges);
   dvEdge = yakl_create_real("dvEdge", nEdges);
   normalThicknessFlux = yakl_create_real("normalThicknessFlux", nVertLevels, nEdges);
   advMaskHighOrder = yakl_create_real("advMaskHighOrder", nVertLevels, nEdges);


   for ( int n = 0; n < nthreads; n++ )
   {
       streams.push_back(yakl::create_stream());
      auto nstr = std::to_string(n);
      wgtTmp[n].field = yakl_create_real(std::string("wgtTmp") + nstr, nVertLevels, nEdges);
      sgnTmp[n].field = yakl_create_real(std::string("sgnTmp") + nstr, nVertLevels, nEdges);
      lowOrderFlx[n].field = yakl_create_real(std::string("lowOrderFlx") + nstr, nVertLevels, nEdges);
      highOrderFlx[n].field = yakl_create_real(std::string("highOrderFlx") + nstr, nVertLevels, nEdges);
   }

   for ( int n = 0; n < ntracers; n++ ) {
      tracers[n].field = yakl_create_real("tracers", nVertLevels, nCells);
   }

   advCoefs = yakl_create_real("advCoefs", nAdv, nEdges);
   advCoefs3rd = yakl_create_real("advCoefs3rd", nAdv, nEdges);
   advCellsForEdge = yakl_create_int("advCellsForEdge", nAdv, nEdges);

}


void init_fields() {
   YAKL_SCOPE(nAdv, mpas::nAdv);
   YAKL_SCOPE(nVertLevels, ::nVertLevels);
   YAKL_LOCAL(minLevelEdgeBot);
   YAKL_LOCAL(maxLevelEdgeTop);

   unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
   std::default_random_engine rgen(seed);

   std::uniform_real_distribution<double> distro(0.0,1.0);

   h_int_1d_t  h_minLevelCell("h_minLevelCell", nCells);
   h_int_1d_t  h_maxLevelCell("h_maxLevelCell", nCells);
   h_int_2d_t  h_advCellsForEdge("temp", nAdv, nEdges);
   h_int_2d_t  h_cellsOnEdge("h_cellsOnEdge", 2, nEdges);
   h_double_1d_t  h_dvEdge("h_dvEdge", nEdges);
   h_double_2d_t  h_advCoefs("h_advCoefs", nAdv, nEdges);
   h_double_2d_t  h_advCoefs3rd("h_advCoefs3rd", nAdv, nEdges);
   h_double_2d_t h_normalThicknessFlux("h_normalThicknessFlux", nVertLevels, nEdges);

   YAKL_LOCAL(minLevelCell);

	/*
   for ( int iCell = 1; iCell <= nCells; ++iCell) {
      h_minLevelCell(iCell) = 1;
   }
   h_minLevelCell.deep_copy_to(*mpas::minLevelCell);
	*/

   yakl::fortran::parallel_for( nCells ,
   YAKL_LAMBDA(int iCell) {
      minLevelCell(iCell) = 1;
   });
	
   yakl::fortran::parallel_for( yakl::fortran::Bounds<1>(nEdges) ,
   YAKL_LAMBDA(int iEdge) {
      minLevelEdgeBot(iEdge) = 1;
      maxLevelEdgeTop(iEdge) = nVertLevels;
   });

   for ( int iCell = 1; iCell <= nCells; ++iCell) {
      int k = int(distro(rgen) * nVertLevels*2.0);
      h_maxLevelCell(iCell) = std::min(std::max(3,k),nVertLevels);
   }
   h_maxLevelCell.deep_copy_to(*mpas::maxLevelCell);

   // Create a random connectivity between cells and edges
   // Also create random values for advective weights
   for ( int iEdge = 1; iEdge <= nEdges; ++iEdge ) {
   	  h_dvEdge(iEdge) = 10.0 * distro(rgen);
   	  h_cellsOnEdge(1,iEdge) = int(nCells * distro(rgen)) + 1;
   	  h_cellsOnEdge(2,iEdge) = int(nCells * distro(rgen)) + 1;
      for ( int i = 1; i <= nAdv; ++i ) {
         h_advCellsForEdge(i,iEdge) = int(nCells * distro(rgen)) + 1;
          assert( h_advCellsForEdge(i,iEdge) >= 1);
          assert( h_advCellsForEdge(i,iEdge) <= nCells);
         h_advCoefs   (i,iEdge) = 20.0 * distro(rgen);
         h_advCoefs3rd(i,iEdge) = 21.0 * distro(rgen);
      }

      for ( int k = 1; k <= nVertLevels; ++k ) {
         h_normalThicknessFlux(k,iEdge) = 15.0*(0.5 - distro(rgen));
      }
   }

   YAKL_LOCAL(advMaskHighOrder);
   
   yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nEdges},{1,nVertLevels}) ,
   YAKL_LAMBDA(int iEdge,int k) {
      advMaskHighOrder(k,iEdge) = 1.0;
   });

   h_dvEdge.deep_copy_to(*mpas::dvEdge);
   h_advCellsForEdge.deep_copy_to(*mpas::advCellsForEdge);
   h_cellsOnEdge.deep_copy_to(*mpas::cellsOnEdge);
   h_advCoefs.deep_copy_to(*(mpas::advCoefs));
   h_advCoefs3rd.deep_copy_to(*(mpas::advCoefs3rd));
   h_normalThicknessFlux.deep_copy_to(*(mpas::normalThicknessFlux));

   YAKL_LOCAL(nAdvCellsForEdge);
    yakl::fortran::parallel_for( yakl::fortran::Bounds<1>(nEdges) ,
    YAKL_LAMBDA(int iEdge) {
      nAdvCellsForEdge(iEdge) = nAdv;
    });

}


void compute_tracer_flx(d_double_2d_t & tracerCur, int tid) {

   YAKL_LOCAL_TID(wgtTmp, tid);
   YAKL_LOCAL_TID(sgnTmp, tid);
   YAKL_LOCAL_TID(lowOrderFlx, tid);
   YAKL_LOCAL_TID(highOrderFlx, tid);
   YAKL_LOCAL(nAdvCellsForEdge);
   YAKL_LOCAL(cellsOnEdge);
   YAKL_LOCAL(advMaskHighOrder);
   YAKL_LOCAL(advCellsForEdge);
   YAKL_LOCAL(minLevelCell);
   YAKL_LOCAL(maxLevelCell);
   YAKL_LOCAL(advCoefs);
   YAKL_LOCAL(advCoefs3rd);
   YAKL_LOCAL(normalThicknessFlux);
   YAKL_LOCAL(minLevelEdgeBot);
   YAKL_LOCAL(maxLevelEdgeTop);
   YAKL_LOCAL(dvEdge);
   YAKL_SCOPE(coef3rdOrder, mpas::coef3rdOrder);

    /*
    yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nCells},{1,nVertLevels}) ,
    YAKL_LAMBDA(int iCell,int k)
    {
        tracerCur(k,iCell) = 0.0;
    }, yakl::LaunchConfig<>());
    yakl::fence();
    return;
    
    yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nEdges},{1,nVertLevels}) ,
                                YAKL_LAMBDA(int iEdge,int k)
    {
        wgtTmp(k,iEdge) = normalThicknessFlux(k,iEdge)*advMaskHighOrder(k,iEdge);
        sgnTmp(k,iEdge) = SIGN(1.0, normalThicknessFlux(k,iEdge));
       highOrderFlx(k,iEdge) = 0.0;
   }, yakl::LaunchConfig<>());
    //yakl::fence();

    yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nEdges},{1,nVertLevels}) ,
    YAKL_LAMBDA(int iEdge,int k)
    {
        for ( int i = 1; i <= nAdvCellsForEdge(iEdge); ++i )
        {
            int iCell = advCellsForEdge(i,iEdge);
            highOrderFlx(k,iEdge) = highOrderFlx(k,iEdge) + tracerCur(k,iCell)*
                       wgtTmp(k,iEdge)*(1 + sgnTmp(k,iEdge));
        }
    }, yakl::LaunchConfig<>());
    return;
     */

    
    yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nEdges},{1,nVertLevels}) ,
    YAKL_LAMBDA(int iEdge,int k)
    {
        // compute some common intermediate factors
        wgtTmp(k,iEdge) = normalThicknessFlux(k,iEdge)*advMaskHighOrder(k,iEdge);
        sgnTmp(k,iEdge) = SIGN(1.0, normalThicknessFlux(k,iEdge));
        highOrderFlx(k,iEdge) = 0.0;
    }, yakl::DefaultLaunchConfig().set_stream(streams[tid]));

    yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nEdges},{1,nVertLevels}) ,
    YAKL_LAMBDA(int iEdge,int k)
    {
        for ( int i = 1; i <= nAdvCellsForEdge(iEdge); ++i )
        {
            int iCell = advCellsForEdge(i,iEdge);
            double coef1 = advCoefs       (i,iEdge);
            double coef3 = advCoefs3rd    (i,iEdge)*coef3rdOrder;
            if ( (k >= minLevelCell(iCell)) && (k <= maxLevelCell(iCell)) )
            {
               highOrderFlx(k,iEdge) = highOrderFlx(k,iEdge) + tracerCur(k,iCell)*
                          wgtTmp(k,iEdge)*(coef1 + coef3*sgnTmp(k,iEdge));
            }
        }
    }, yakl::DefaultLaunchConfig().set_stream(streams[tid]));

    yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nEdges},{1,nVertLevels}) ,
    YAKL_LAMBDA(int iEdge,int k)
    {
        if ( (k >= minLevelEdgeBot(iEdge)) &&  (k<= maxLevelEdgeTop(iEdge)) )
        {
            int cell1 = cellsOnEdge(1, iEdge);
            int cell2 = cellsOnEdge(2, iEdge);

            // Compute 2nd order fluxes where needed.
            // Also compute low order upwind horizontal flux (monotonic and diffused)
            // Remove low order flux from the high order flux
            // Store left over high order flux in highOrderFlx array
            double tracerWeight = (1.0 - advMaskHighOrder(k,iEdge))
                         * (dvEdge(iEdge) * 0.5)
                         * normalThicknessFlux(k, iEdge);

            lowOrderFlx(k,iEdge) = dvEdge(iEdge) * 
               (DMAX(0.0,normalThicknessFlux(k,iEdge))*tracerCur(k,cell1)
              + DMIN(0.0,normalThicknessFlux(k,iEdge))*tracerCur(k,cell2));

            highOrderFlx(k, iEdge) = highOrderFlx(k, iEdge)
                                   + tracerWeight * (tracerCur(k, cell1)
                                                   + tracerCur(k, cell2));

            highOrderFlx(k,iEdge) = highOrderFlx(k,iEdge)
                                  -  lowOrderFlx(k,iEdge);
        }
    }, yakl::DefaultLaunchConfig().set_stream(streams[tid]));

}

//
// extracted kernel of this loop
//

/**
         do iEdge = 1, nEdges
           cell1 = cellsOnEdge(1, iEdge)
           cell2 = cellsOnEdge(2, iEdge)

           ! compute some common intermediate factors
           do k = 1, nVertLevels
              wgtTmp(k) = normalThicknessFlux(k,iEdge)* &
                          advMaskHighOrder(k,iEdge)
              sgnTmp(k) = sign(1.0_RKIND, &
                               normalThicknessFlux(k,iEdge))
              flxTmp(k) = 0.0_RKIND
           end do

           ! Compute 3rd or 4th fluxes where requested.
           do i = 1, nAdvCellsForEdge(iEdge)
              iCell = advCellsForEdge(i,iEdge)
              coef1 = advCoefs       (i,iEdge)
              coef3 = advCoefs3rd    (i,iEdge)*coef3rdOrder
              do k = minLevelCell(iCell), maxLevelCell(iCell)
                 flxTmp(k) = flxTmp(k) + tracerCur(k,iCell)* &
                             wgtTmp(k)*(coef1 + coef3*sgnTmp(k))
              end do ! k loop
           end do ! i loop over nAdvCellsForEdge

           do k=1,nVertLevels
              highOrderFlx(k,iEdge) = flxTmp(k)
            !if ( abs(highOrderFlx(k,iEdge) - tmp_highOrderFlx(k,iEdge)) > 1e-12) &
             !print *,mpas_myrank,": highOrderFlx expected ", highOrderFlx(k,iEdge), ", but got ",tmp_highOrderFlx(k,iEdge), &
              !   " for k,iEdge = ",k,iEdge
           end do

           ! Compute 2nd order fluxes where needed.
           ! Also compute low order upwind horizontal flux (monotonic)
           ! Remove low order flux from the high order flux
           ! Store left over high order flux in highOrderFlx array
           do k = minLevelEdgeBot(iEdge), maxLevelEdgeTop(iEdge)
              tracerWeight = (1.0_RKIND - advMaskHighOrder(k,iEdge)) &
                           * (dvEdge(iEdge) * 0.5_RKIND)             &
                           * normalThicknessFlux(k, iEdge)

              lowOrderFlx(k,iEdge) = dvEdge(iEdge) * &
               (max(0.0_RKIND,normalThicknessFlux(k,iEdge))*tracerCur(k,cell1) &
              + min(0.0_RKIND,normalThicknessFlux(k,iEdge))*tracerCur(k,cell2))

              highOrderFlx(k,iEdge) = highOrderFlx(k,iedge) &
                                    + tracerWeight*(tracerCur(k,cell1) &
                                                  + tracerCur(k,cell2))

              highOrderFlx(k,iEdge) = highOrderFlx(k,iEdge) &
                                    -  lowOrderFlx(k,iEdge)
           end do ! k loop
        end do ! iEdge loop
*/


int main(int, const char **)
{
	yakl::init();
	
   alloc_fields();
   init_fields();

	std::cerr << " Running with " << omp_get_max_threads() << " OMP threads..." << std::endl;
    Timer timer;
    
    timer.start();
   #pragma omp parallel for
   for ( int tr = 0; tr < ntracers; ++tr ) {
      int tid = omp_get_thread_num();
	std::cerr << " tid, tr = " << tid << " " << tr << std::endl;
      compute_tracer_flx(*(mpas::tracers[tr].field), tid);
   }

   yakl::fence();
    timer.stop();
    
    std::cout << "Seconds: " << timer.elapsedSeconds() << std::endl;

		//yakl::finalize();

   return 0;
}
