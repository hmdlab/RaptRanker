/*
MIT License

Copyright (c) 2019 Hamada Laboratory, and Ryoga Ishida

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#ifndef SIMPLE_TIMER_H
#define SIMPLE_TIMER_H

#include <stdio.h>
#include <iostream>
#include <chrono>
#include <string>
#include <vector>

class Stopwatch {
    public:
        explicit Stopwatch(bool display_result = true)
            :display_result_(display_result),
            start_time_(std::chrono::high_resolution_clock::now() ),
            end_time_(start_time_) {};

        ~Stopwatch(){
            if(display_result_) {
                if(end_time_ == start_time_){
                    end_time_ = std::chrono::high_resolution_clock::now();
                    std::cout << "Stopwatch : The measurement was finished automatically." << std::endl;
                }

                std::cout << "======Stopwatch : measurement result ======" << std::endl;
                if(!check_points_.empty()){
                    //split_times (start->check1, start->check2, ....)
                    std::cout << "Split times are :" << std::endl;
                    for (unsigned long i = 0; i < check_points_.size(); ++i) {
                        std::cout << "\t Split " << i << "\t" << check_points_[i].first << "\t";
                        output_elapsed_time(start_time_, check_points_[i].second);
                    }
                    std::cout << "\t Finish ";
                    output_elapsed_time(start_time_, end_time_);


                    //lap_times (start->check1, check1->check2, ....)
                    std::cout << "Lap times are :" << std::endl;
                    std::cout << "\t Lap 1(start ~ "<< check_points_[0].first << ")  ";
                    output_elapsed_time(start_time_, check_points_[0].second);
                    for (unsigned long i = 0; i < check_points_.size() - 1; ++i) {
                        std::cout << "\t Lap " << i+2 << "(" << check_points_[i].first << " ~ " << check_points_[i+1].first << ")  ";
                        output_elapsed_time(check_points_[i].second, check_points_[i+1].second);
                    }
                    std::cout << "\t FinalLap " << "(" << check_points_[check_points_.size()-1].first << " ~ Finish)  ";
                    output_elapsed_time(check_points_[check_points_.size()-1].second, end_time_);
                }


                std::cout << "Elapsed time is :" << std::endl << "\t";
                output_elapsed_time(start_time_, end_time_);
                std::cout << "===========================================" << std::endl;
            }
        };

        void add_check_point(std::string check_point_name){
            auto check_point = std::chrono::high_resolution_clock::now();
            check_points_.emplace_back(check_point_name, check_point);
        };

        void finish_measurement(){
            end_time_= std::chrono::high_resolution_clock::now();
        };

        std::chrono::high_resolution_clock::duration get_elapsed_duration(){
            if (end_time_ == start_time_){
                end_time_ = std::chrono::high_resolution_clock::now();
                std::cout << "Stopwatch : The measurement was finished automatically." << std::endl;
            }
            return(end_time_ - start_time_);
        }

        std::vector<std::chrono::high_resolution_clock::duration> get_split_durations(){
            //split_times (start->check1, start->check2, ....)
            if (end_time_ == start_time_){
                end_time_ = std::chrono::high_resolution_clock::now();
                std::cout << "Stopwatch : The measurement was finished automatically." << std::endl;
            }

            std::vector<std::chrono::high_resolution_clock::duration> split_durations{};

            if(check_points_.empty()){
                std::cerr << "Stopwatch : There are no check points." << std::endl;
            }else{
                for (const auto& each_check_points : check_points_) {
                    split_durations.emplace_back(each_check_points.second - start_time_);
                }
            }
            split_durations.emplace_back(end_time_ - start_time_);

            return(split_durations);
        };

        std::vector<std::chrono::high_resolution_clock::duration> get_lap_durations(){
        //lap_times (start->check1, check1->check2, ....)
        if (end_time_ == start_time_){
            end_time_ = std::chrono::high_resolution_clock::now();
            std::cout << "Stopwatch : The measurement was finished automatically." << std::endl;
        }

        std::vector<std::chrono::high_resolution_clock::duration> lap_durations{};

        if(check_points_.empty()){
            std::cerr << "Stopwatch : There are no check points." << std::endl;
        }else{
            lap_durations.emplace_back(check_points_[0].second - start_time_);
            for (unsigned long i = 0; i < check_points_.size() - 1; ++i) {
                lap_durations.emplace_back(check_points_[i+1].second - check_points_[i].second);
            }
        }
        lap_durations.emplace_back(end_time_ - check_points_[check_points_.size()-1].second);

        return(lap_durations);
    };


    private:
        bool display_result_;
        std::chrono::high_resolution_clock::time_point start_time_;
        std::chrono::high_resolution_clock::time_point end_time_;
        std::vector<std::pair<std::string, std::chrono::high_resolution_clock::time_point> > check_points_;

    void output_elapsed_time(std::chrono::high_resolution_clock::time_point start_time, std::chrono::high_resolution_clock::time_point end_time) {
            using namespace std::chrono;
            auto elapsed_time = end_time - start_time;


            auto nano_time = duration_cast<nanoseconds>(elapsed_time).count();
            if(nano_time > 1 && nano_time < 1000){
                std::cout << nano_time << "[ns] (" << (double)nano_time/1000 << "[µs])" << std::endl;
            }

            auto micro_time = duration_cast<microseconds>(elapsed_time).count();
            if(micro_time > 1 && micro_time < 1000){
                std::cout <<  micro_time << "[µs] (" << (double)micro_time/1000 << "[ms])" << std::endl;
            }

            auto milli_time = duration_cast<milliseconds>(elapsed_time).count();
            if(milli_time > 1 && milli_time < 60000){
                std::cout << milli_time << "[ms] (" << (double)milli_time/1000 << "[s])" << std::endl;
            }

            if(milli_time > 60000){
                auto sec_time = duration_cast<seconds>(elapsed_time).count() % 60;
                auto min_time = duration_cast<minutes>(elapsed_time).count() % 60;
                auto hour_time = duration_cast<hours>(elapsed_time).count();
                std::cout << hour_time << ":" << min_time << ":" << sec_time << std::endl;
            }
        };


};


#endif //SIMPLE_TIMER_H
