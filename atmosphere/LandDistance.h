#ifndef LAND_DISTANCE_H
#define LAND_DISTANCE_H

#include <cmath>
#include <queue>
#include <utility>
#include <vector>

namespace AtomLand {
    // Distance [m] of every column (index j * km + k) to the nearest ocean column (i_topography == 0), 0 over ocean.
    // Multi-source Dijkstra on the lat-lon grid (8 neighbours, periodic in longitude), ocean columns = sources.
    // A copy of the block in initWaterWapour (ATM_RH_LAND, 2026-10-02), for ATM_MC_T_ADD_LAND_L (2026-10-03).
    inline std::vector<double> distanceToOcean(const std::vector<std::vector<int> >& i_topography, int jm, int km) {
        std::vector<double> d_ocean;
        const double a = 6.371e6, dlat = M_PI / (jm - 1), dlon = 2.0 * M_PI / (km - 1);
        d_ocean.assign((size_t)jm * km, 1e30);
        using QE = std::pair<double, int>;
        std::priority_queue<QE, std::vector<QE>, std::greater<QE>> pq;
        for (int j = 0; j < jm; j++) for (int k = 0; k < km; k++)
            if (i_topography[j][k] == 0) { d_ocean[(size_t)j * km + k] = 0.0; pq.push({0.0, j * km + k}); }
        while (!pq.empty()) {
            auto [dd, id] = pq.top(); pq.pop();
            if (dd > d_ocean[id]) continue;
            const int j = id / km, k = id % km;
            for (int dj = -1; dj <= 1; dj++) for (int dk = -1; dk <= 1; dk++) {
                if (!dj && !dk) continue;
                const int jj = j + dj; if (jj < 0 || jj >= jm) continue;
                const int kk = ((k + dk) % (km - 1) + (km - 1)) % (km - 1);
                const double cl = cos((90.0 - 0.5 * (j + jj) * 180.0 / (jm - 1)) * M_PI / 180.0);
                const double step = a * std::sqrt((dj * dlat) * (dj * dlat) + (dk * dlon * cl) * (dk * dlon * cl));
                const size_t nid = (size_t)jj * km + kk;
                if (dd + step < d_ocean[nid]) { d_ocean[nid] = dd + step; pq.push({d_ocean[nid], (int)nid}); }
            }
        }
        for (int j = 0; j < jm; j++) d_ocean[(size_t)j * km + km - 1] = d_ocean[(size_t)j * km];   // seam copy
        return d_ocean;
    }
}
#endif
