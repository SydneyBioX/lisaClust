#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "lisaclust/core.hpp"

// The local indicators of inhomLocalK(), which before lisaClust 1.21.1 were computed in R with
// spatstat.geom::closepairs(), cut(), dplyr and data.table. This reproduces that computation:
// - ordered pairs (i, j), i != j, with distance <= the largest radius, the distance put in band k when
//   Rs[k] < d <= Rs[k + 1] (cut(d, Rs, include.lowest = TRUE)); with includeSelf (from 1.21.2) also the pair
//   (i, i) at distance 0, so that a cell counts itself;
// - the pairs' weights wt[j] summed per (i, band, type of j), then accumulated over the bands of the
//   image (data.table's CJ() crosses only the cells, bands and types that occur among the pairs, and
//   both sums accumulate in long double, as data.table's sum() and base cumsum() do);
// - the expectation uses the cell's edge share in that band when the cell has a pair of that type in
//   the band, and 1 otherwise (the value of the unmatched rows of the join).
lisaclust::LocalCurves lisaclust::localCurves(const double* x, const double* y, const int* type, int n, int K,
                                              const std::vector<double>& Rs, const std::vector<double>& labelVal,
                                              const double* wt, const std::vector<double>& lam, const double* edge,
                                              bool Lfunction, bool includeSelf) {
  LocalCurves out;
  const int nb = static_cast<int>(Rs.size()) - 1;
  out.n = n;
  out.nb = nb;
  out.K = K;
  out.cellPresent.assign(n, 0);
  out.binPresent.assign(std::max(nb, 0), 0);
  out.typePresent.assign(K, 0);
  if (n == 0 || nb < 1) return out;
  const std::size_t size = static_cast<std::size_t>(n) * nb * K;
  std::vector<long double> S(size, 0.0L);
  std::vector<char> has(size, 0);
  auto at = [&](int i, int k, int J) { return (static_cast<std::size_t>(k) * K + J) * n + i; };

  const double rmax = Rs[nb], r2max = rmax * rmax;
  // Uniform grid with cells of side h, so every neighbour within rmax lies in the block of
  // (2 * reach + 1)^2 grid cells around a point. A zero largest radius still finds coincident points.
  const double x0 = *std::min_element(x, x + n), x1 = *std::max_element(x, x + n);
  const double y0 = *std::min_element(y, y + n), y1 = *std::max_element(y, y + n);
  // The grid has at most about 2048 cells a side, however small the radius is relative to the image.
  const double span = std::max(x1 - x0, y1 - y0);
  const double h = std::max({rmax / 2, span / 2048, 1e-12});
  const int reach = std::max(1, static_cast<int>(std::ceil(rmax / h)));
  const int gx = static_cast<int>((x1 - x0) / h) + 1;
  const int gy = static_cast<int>((y1 - y0) / h) + 1;
  std::vector<int> cellOf(n), start(static_cast<std::size_t>(gx) * gy + 1, 0), order(n);
  for (int i = 0; i < n; ++i) {
    const int cx = std::min(static_cast<int>((x[i] - x0) / h), gx - 1);
    const int cy = std::min(static_cast<int>((y[i] - y0) / h), gy - 1);
    cellOf[i] = cy * gx + cx;
    ++start[cellOf[i] + 1];
  }
  for (std::size_t c = 1; c < start.size(); ++c) start[c] += start[c - 1];
  {
    std::vector<int> fill(start.begin(), start.end() - 1);
    for (int i = 0; i < n; ++i) order[fill[cellOf[i]]++] = i;
  }

  auto add = [&](int i, int j, int k) {
    const std::size_t a = at(i, k, type[j]);
    S[a] += wt[j];
    has[a] = 1;
    out.cellPresent[i] = 1;
    out.binPresent[k] = 1;
    out.typePresent[type[j]] = 1;
  };
  for (int cy = 0; cy < gy; ++cy) {
    for (int cx = 0; cx < gx; ++cx) {
      const int c = cy * gx + cx;
      for (int si = start[c]; si < start[c + 1]; ++si) {
        const int i = order[si];
        // Each unordered pair once: later points of this grid cell, and the grid cells after it.
        for (int ny = cy; ny <= std::min(cy + reach, gy - 1); ++ny) {
          const int lo = ny == cy ? cx : std::max(cx - reach, 0);
          const int hi = std::min(cx + reach, gx - 1);
          for (int nx = lo; nx <= hi; ++nx) {
            const int nc = ny * gx + nx;
            for (int sj = (nc == c ? si + 1 : start[nc]); sj < start[nc + 1]; ++sj) {
              const int j = order[sj];
              const double dx = x[j] - x[i], dy = y[j] - y[i];
              const double d2 = dx * dx + dy * dy;
              if (d2 > r2max) continue;
              const double d = std::sqrt(d2);
              int k = 0;
              while (k < nb && d > Rs[k + 1]) ++k;
              if (k == nb) continue;
              add(i, j, k);
              add(j, i, k);
            }
          }
        }
      }
    }
  }

  // Each cell is one of its own neighbours, at distance 0 (Patrick et al. 2023 sum over every cell of the type,
  // the cell itself included).
  if (includeSelf)
    for (int i = 0; i < n; ++i) add(i, i, 0);

  out.value.assign(size, std::numeric_limits<double>::quiet_NaN());
  for (int J = 0; J < K; ++J) {
    if (!out.typePresent[J]) continue;
    for (int i = 0; i < n; ++i) {
      if (!out.cellPresent[i]) continue;
      long double cum = 0.0L;
      for (int k = 0; k < nb; ++k) {
        if (!out.binPresent[k]) continue;
        const std::size_t a = at(i, k, J);
        cum += static_cast<double>(S[a]);
        const double sum = static_cast<double>(cum);
        const double e = has[a] ? edge[static_cast<std::size_t>(k) * n + i] : 1.0;
        const double E = labelVal[k] * labelVal[k] * M_PI * e * lam[J];
        out.value[a] = Lfunction ? std::sqrt(sum) - std::sqrt(E) : (sum - E) / std::sqrt(E);
      }
    }
  }
  return out;
}

// Nearest training point by a uniform grid search, growing the search ring until it can hold no closer
// point. All points at exactly the nearest distance vote; ties go to the smallest label.
std::vector<int> lisaclust::nearestLabels(const double* tx, const double* ty, const int* label, int nTrain,
                                          const double* qx, const double* qy, int nQuery) {
  std::vector<int> out(nQuery, -1);
  if (nTrain == 0) return out;
  const double x0 = *std::min_element(tx, tx + nTrain), x1 = *std::max_element(tx, tx + nTrain);
  const double y0 = *std::min_element(ty, ty + nTrain), y1 = *std::max_element(ty, ty + nTrain);
  const double w = std::max(x1 - x0, 1e-12), hgt = std::max(y1 - y0, 1e-12);
  // about one training point per grid cell, with at most about 1024 cells a side
  const double h = std::max({std::sqrt(w * hgt / nTrain), std::max(w, hgt) / 1024, 1e-12});
  const int gx = static_cast<int>(w / h) + 1, gy = static_cast<int>(hgt / h) + 1;
  std::vector<int> start(static_cast<std::size_t>(gx) * gy + 1, 0), order(nTrain), cellOf(nTrain);
  auto gridX = [&](double v) { return std::min(std::max(static_cast<int>((v - x0) / h), 0), gx - 1); };
  auto gridY = [&](double v) { return std::min(std::max(static_cast<int>((v - y0) / h), 0), gy - 1); };
  for (int i = 0; i < nTrain; ++i) {
    cellOf[i] = gridY(ty[i]) * gx + gridX(tx[i]);
    ++start[cellOf[i] + 1];
  }
  for (std::size_t c = 1; c < start.size(); ++c) start[c] += start[c - 1];
  {
    std::vector<int> fill(start.begin(), start.end() - 1);
    for (int i = 0; i < nTrain; ++i) order[fill[cellOf[i]]++] = i;
  }
  int maxLabel = 0;
  for (int i = 0; i < nTrain; ++i) maxLabel = std::max(maxLabel, label[i]);
  std::vector<int> votes(maxLabel + 1, 0), voted;
  for (int q = 0; q < nQuery; ++q) {
    const double px = qx[q], py = qy[q];
    const int cx = gridX(px), cy = gridY(py);
    // distance from the query to the outside of its own grid cell block of radius ring
    double best = std::numeric_limits<double>::infinity();
    voted.clear();
    for (int ring = 0;; ++ring) {
      for (int ny = cy - ring; ny <= cy + ring; ++ny) {
        if (ny < 0 || ny >= gy) continue;
        const bool edgeRow = ny == cy - ring || ny == cy + ring;
        for (int nx = cx - ring; nx <= cx + ring; nx += (edgeRow || ring == 0) ? 1 : 2 * ring) {
          if (nx < 0 || nx >= gx) continue;
          const int c = ny * gx + nx;
          for (int s = start[c]; s < start[c + 1]; ++s) {
            const int i = order[s];
            const double dx = tx[i] - px, dy = ty[i] - py, d2 = dx * dx + dy * dy;
            if (d2 < best) {
              best = d2;
              for (int v : voted) votes[v] = 0;
              voted.clear();
            }
            if (d2 == best) {
              if (votes[label[i]]++ == 0) voted.push_back(label[i]);
            }
          }
        }
      }
      // Points outside the searched block are at least this far from the query.
      const double gap = std::min({px - (x0 + (cx - ring) * h), (x0 + (cx + ring + 1) * h) - px,
                                   py - (y0 + (cy - ring) * h), (y0 + (cy + ring + 1) * h) - py});
      const bool covered = cx - ring <= 0 && cy - ring <= 0 && cx + ring >= gx - 1 && cy + ring >= gy - 1;
      if (covered || (best < std::numeric_limits<double>::infinity() && gap > 0 && gap * gap > best)) break;
    }
    int lab = -1, most = 0;
    for (int v : voted)
      if (votes[v] > most || (votes[v] == most && v < lab)) { most = votes[v]; lab = v; }
    for (int v : voted) votes[v] = 0;
    out[q] = lab;
  }
  return out;
}
