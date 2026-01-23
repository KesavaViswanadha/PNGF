// Magic Hat random number drawing class.
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
// Helper class to draw a set of random numbers without replacement.
// The constructor initializes a MT19937 RNG with either a provided
// seed or generates a random RNG seed when the provided seed is 0.
//
// Known issues: The std::uniform_int_distribution class seems to behave
// differently on Mac vs. Linux-- the same RNG seed will lead to different
// uniform random values being drawn on a Mac and a Linux machine. To
// reproduce the exact Substrate Antenna (SA) device example reported in
// the manuscript, a Linux machine should be used with the default RNG seed.
//
// Constantine Sideris (sideris@stanford.edu), Jui-Hung Sun (juihungs@usc.edu)
//

class MagicHat
{

public:
    std::random_device rd;  // Seed source.
    std::mt19937 gen; // Mersenne Twister engine.
    std::uniform_int_distribution<> dis;
    int rng_seed;

    int total_entries;
    int num_left;
    std::vector<int> entry_idx;

    // Initialize the "magic" hat.
    // Inputs:
    // seed: Specify the RNG seed, if 0, compute a random one.
    // run_index: If running many optimizations in parallel, this
    // specifies the ID of the current run. It is used to help
    // randomize the initial RNG seed further. Can specify 0 if
    // not needed.
    // num_entries: Number of entries to put in the "hat".
    MagicHat (int seed, int run_index, int num_entries)
    {
        if (seed==0)
        {
            rng_seed = rd();
            gen.seed(rng_seed);
            int trnd = 0;
            for (int i = 0; i < run_index+1; i++)
                trnd += rd();
            rng_seed = rng_seed ^ (trnd * 0xc0decafe) ^ (run_index * 0x12345678);
            gen.seed(rng_seed);
        } else
            rng_seed = seed;
        std::cout << "RNG seed: " << rng_seed << std::endl;
        gen.seed(rng_seed);
        dis = std::uniform_int_distribution<>(0, num_entries-1); // Uniform distribution in [0, N].
        total_entries = num_entries;
        num_left = num_entries;

        entry_idx.resize(num_entries);
        for (int i = 0; i < num_entries; i++)
            entry_idx[i] = i;
    }

    // return a random binary bit (0 or 1):
    int rand_binary(void)
    {
        return gen()&1;
    }

    // Draw an entry from the "hat" without replacement:
    int draw(void)
    {
        int val, idx;

        if (num_left<=0)
            return -1;

        dis.param(std::uniform_int_distribution<int>::param_type(0, num_left-1));
        idx = dis(gen);
        val = entry_idx[idx];
        entry_idx[idx] = entry_idx[num_left-1];
        num_left--;

        return val;
    }

    // Reset the hat and fill it back up with all of the entries:
    int reset(void)
    {
        for (int i = 0; i < total_entries; i++)
            entry_idx[i] = i;
        num_left = total_entries;
        return num_left;
    }

    ~MagicHat()
    {
    }
};