// Helper class to accurately time operations.
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

#pragma once

#include <thread>
#include <chrono>


class acmetime
{
    public:
        double period;
        std::chrono::time_point<std::chrono::high_resolution_clock> s_time, e_time;

        acmetime(void)
        {
            // number of seconds per tick of the clock
            period = (double)std::chrono::high_resolution_clock::period::num /
                    (double)std::chrono::high_resolution_clock::period::den;
        }

        // start the timer
        void tik(void)
        {
            s_time = std::chrono::high_resolution_clock::now();
        }

        // stop the timer and return elapsed time in seconds
        double tok (void)
        {
            e_time = std::chrono::high_resolution_clock::now();
        
            long int diff = e_time.time_since_epoch().count() - s_time.time_since_epoch().count();

            return (double)diff * period;
        }

        ~acmetime()
        {
            return;
        }
};
