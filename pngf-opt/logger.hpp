// Logger helper class.
//
// Copyright (c) 2025, 2026, Constantine Sideris (sideris@stanford.edu) and Jui-Hung Sun
// (juihungs@usc.edu)
// 
// This program is free software: you can redistribute it and/or modify it under the terms 
// of the GNU Affero General Public License as published by the Free Software Foundation, 
// either version 3 of the License, or (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful, but WITHOUT ANY 
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
// PARTICULAR PURPOSE. See the GNU Affero General Public License for more details.
// 
// You should have received a copy of the GNU Affero General Public License along with 
// this program. If not, see <https://www.gnu.org/licenses/>. 
//
// 
// 
// Helper class to log information about the design being optimized
// during each iteration. The data is stored in memory and only
// exported to a file at the end using the "save" function.
//
// Constantine Sideris (sideris@stanford.edu), Jui-Hung Sun (juihungs@usc.edu)
//

class Logger
{
    int rng_seed;
    std::vector<int> flip_history;
    std::vector<precision> perf;
    std::vector<int> flip_decision;
    std::vector<int> tilemap_start, tilemap_end;
    int num_perf_metrics;
    int tot_flips;
    int num_flips;
    int tilemap_size;

    public:
        Logger(int seed, int tilemap_sz, int max_flips, int num_metrics)
        {
            rng_seed = seed;
            tilemap_size = tilemap_sz;
            flip_history.resize(max_flips);
            flip_decision.resize(max_flips);
            tot_flips = max_flips;
            num_flips = 0;
            num_perf_metrics = num_metrics;
            perf.resize(max_flips*num_metrics);
            tilemap_start.resize(tilemap_sz);
            tilemap_end.resize(tilemap_sz);
        }

        void record_tilemap_init (std::vector<int> tile_map)
        {
            for (int i = 0; i < tilemap_size; i++)
                tilemap_start[i] = tile_map[i];            
        }
        void record_tilemap_final (std::vector<int> tile_map)
        {
            for (int i = 0; i < tilemap_size; i++)
                tilemap_end[i] = tile_map[i];            
        }

        void record_flip (int flip_index, int flip_dec, std::vector<precision> flip_perf)
        {
            if (num_flips>tot_flips)
                return; // ran out of buffer space to record flips!
            flip_history[num_flips] = flip_index;
            flip_decision[num_flips] = flip_dec;
            for (int i = 0; i < num_perf_metrics; i++)
                perf[num_flips*num_perf_metrics+i] = flip_perf[i];
            num_flips++;
        }

        void save (std::string fname)
        {
            std::ofstream ofil(fname);

            ofil << rng_seed << " " << num_flips << std::endl;
            for (int j = 0; j < tilemap_size; j++)
            {
                if (tilemap_start[j])
                    ofil << "x";
                else
                    ofil << "."; 
            }
            ofil << std::endl;
            for (int i = 0; i < num_flips; i++)
            {
                ofil << i << " ";
                ofil << flip_history[i] << " ";
                ofil << flip_decision[i] << " ";
                for (int j = 0; j < num_perf_metrics; j++)
                    ofil << perf[i*num_perf_metrics+j] << " ";
                ofil << std::endl;                
            }
            for (int j = 0; j < tilemap_size; j++)
            {
                if (tilemap_end[j])
                    ofil << "x";
                else
                    ofil << "."; 
            }
            ofil << std::endl;

            ofil.close();
        }
};