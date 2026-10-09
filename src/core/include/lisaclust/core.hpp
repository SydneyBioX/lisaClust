// lisaClust's numeric core: plain C++17, free of R and Python objects. It is shared verbatim by the R
// package (src/bindings.cpp) and the Python package lisaclust (src/bindings.cpp there), which copies
// this folder with tools/sync_core.sh.
#ifndef LISACLUST_CORE_HPP
#define LISACLUST_CORE_HPP

#include <cstddef>
#include <vector>

namespace lisaclust {

struct Pt {
  double x, y;
};
using Ring = std::vector<Pt>;

// Area of the intersection of the window `rings` with the npoly-gon that spatstat.geom::disc() builds
// around each centre, i.e. area(intersect.owin(disc(r, centre, npoly = npoly), W)). Rings are in
// spatstat orientation (outer anticlockwise, holes clockwise). The same code as spicyR's.
void discAreas(const double* x, const double* y, std::size_t n, double r, int npoly,
               const std::vector<Ring>& rings, double* out);

// Local indicators of one image, as lisaClust's inhomLocalK() defines them.
//   x, y, type (0-based, < K) and wt (the inverse-density weight of each cell as a neighbour), n cells.
//   Rs: the radii with a leading 0, increasing; labelVal[k]: radius k + 1 as written in the column names.
//   lam[J]: cells of type J per unit area. edge: n x (Rs.size() - 1), column-major, the share of each
//   cell's disc inside the window.
//   Lfunction: false for K, (sum - E) / sqrt(E); true for L, sqrt(sum) - sqrt(E).
//   includeSelf: count each cell among its own neighbours, at distance 0 (as Patrick et al. 2023 define the LISA).
// For cell i, radius k and type J, sum is the weighted count of type-J neighbours within Rs[k + 1],
// accumulated over the radii at which the image has any pair, and E = Rs[k + 1]^2 pi e lam[J], where e is
// the cell's edge share at that radius when it has a type-J neighbour in that distance band, else 1.
struct LocalCurves {
  int n = 0, nb = 0, K = 0;
  std::vector<double> value;       // n x K x nb, index (k * K + J) * n + i
  std::vector<char> cellPresent;   // cells with a neighbour within the largest radius
  std::vector<char> binPresent;    // distance bands with a pair
  std::vector<char> typePresent;   // types that are a neighbour of some cell
};
LocalCurves localCurves(const double* x, const double* y, const int* type, int n, int K,
                        const std::vector<double>& Rs, const std::vector<double>& labelVal,
                        const double* wt, const std::vector<double>& lam, const double* edge,
                        bool Lfunction, bool includeSelf);

// The label of the nearest training point to each query point. Points at exactly the nearest distance
// vote; ties between labels go to the smallest label. Labels are 0-based.
std::vector<int> nearestLabels(const double* tx, const double* ty, const int* label, int nTrain,
                               const double* qx, const double* qy, int nQuery);

// Convex hull, anticlockwise, without repeated end point (Andrew's monotone chain).
Ring convexHull(const double* x, const double* y, std::size_t n);

// Distance from each point to the boundary of the window `rings`.
void distanceToBoundary(const double* x, const double* y, std::size_t n, const std::vector<Ring>& rings,
                        double* out);

}  // namespace lisaclust

#endif
