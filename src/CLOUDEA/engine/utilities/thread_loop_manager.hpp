/*
 * CLOUDEA - Software for solving PDEs using explicit methods.
 * Copyright (C) 2017  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */



#ifndef CLOUDEA_UTILITIES_THREAD_LOOP_MANAGER_HPP_
#define CLOUDEA_UTILITIES_THREAD_LOOP_MANAGER_HPP_

/**
   \file thread_loop_manager.hpp
   \brief ThreadLoopManager class header file.
   \author Konstantinos A. Mountris
   \date 06/08/2019
*/


#include <iostream>
#include <cmath>
#include <vector>


namespace CLOUDEA {

/** \addtogroup Utilities \{c*/


/**
 * \class ThreadLoopManager
 * \brief Class implemmenting a manager for multithreaded loops.
 */
class ThreadLoopManager {

private:
    std::vector<std::size_t> start_id_;

    std::vector<std::size_t> end_id_;

public:

    ThreadLoopManager() : start_id_(), end_id_() {}


    ~ThreadLoopManager() {}


    inline void AppendLoopInfo(std::size_t start_id, std::size_t end_id)
    {
        this->start_id_.emplace_back(start_id);
        this->end_id_.emplace_back(end_id);
    }


    inline void Reset()
    {
        this->start_id_.clear();
        this->end_id_.clear();
    }


    inline void SetLoopRanges(std::size_t entries_num, std::size_t threads_num)
    {
        this->Reset();

        // Iteration range per thread.
        if (threads_num <= 1 || threads_num >= entries_num) {
            // One thread through the whole range.
            this->AppendLoopInfo(0, entries_num);
        } else {

            // Multiple ranges across threads.
            auto entries_per_thread = entries_num / threads_num;
            auto remainder = entries_num % threads_num;

            auto loop_start = std::size_t{0};
            auto loop_end = entries_per_thread;
            for (std::size_t t = 0; t != threads_num-1; ++t) {
                this->AppendLoopInfo(loop_start, loop_end);
                loop_start += entries_per_thread;
                loop_end += entries_per_thread;
            }
            this->AppendLoopInfo(loop_start, loop_end+remainder);
        }
    }


    inline auto LoopStartId(std::size_t thread_id) const { return this->start_id_[thread_id]; }


    inline auto LoopEndId(std::size_t thread_id) const { return this->end_id_[thread_id]; }

};


/** \} End of Doxygen Groups */
} //end of namespace CLOUDEA

#endif //CLOUDEA_UTILITIES_THREAD_LOOP_MANAGER_HPP_