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

#ifndef RaptRanker_HPP
#define RaptRanker_HPP


#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <string>
#include <memory>

#include "sqlite3.h"
#include "CapR.h"
#include "SketchSort.hpp"
#include "SQLiteCpp.h"
#include "json11.hpp"

constexpr char key[5] = {'A', 'T', 'G', 'C', '\n'};
constexpr int key_length = 5;
constexpr int ROUND_BUFFER = 100;

//using namespace std;

class RaptRanker{
 public:
    explicit RaptRanker(std::string &parameter_file){
        parameter_file_name_ = parameter_file;
    }
    
    void Run();

    struct Input_file_info{
    public:
        int round_id_;
        std::string round_name_;
        std::string file_path_;
    };
    
    struct Parameter_info{
    public:
        int file_type_;
        std::string forward_primer_;
        std::string reverse_primer_;
        std::string add_forward_primer_;
        std::string add_reverse_primer_;
        int sequence_maximum_length_;
        int sequence_minimum_length_;
        int wide_length_;
        double nucleotide_weight_;
        double cosine_distance_;
        double missing_ratio_;
        int input_file_nums_;
        std::vector<Input_file_info> input_file_infos_;
        std::string experiment_dbfile_;
        std::string analysis_dbfile_;
        std::string analysis_output_path_;


        //for dev
        bool exFiltering_ = true;
        bool exCapR_ = true;
        bool exSketchSort_ = true;
        bool exKmer_ = true;
        bool export_scoreCSV_ = true;
        bool add_binding_ = false;
        std::string bindingfile_path_;
    };


    class WeightedUnionFind {
    private:
        std::vector<int> data_;
        std::vector<int> cluster_ids_;
        std::vector<double> dist_;
    public:
        explicit WeightedUnionFind(unsigned long size): data_(size,-1), cluster_ids_(size,-1),dist_(size,0.0) {};
        ~WeightedUnionFind() = default;

        void unite(int x, int y, double dist) {
            //x + dist = y
            dist += distance(x); dist -= distance(y);
            x = root(x); y = root(y);

            if (x != y) {
                if (data_[y] < data_[x]) {
                    std::swap(x, y); dist = -dist;
                }
                data_[x] += data_[y]; data_[y] = x;
                dist_[y] = dist;
            }
        }
        bool isSameTree(int x, int y) {
            return root(x) == root(y);
        }
        int root(int x) {
            if(data_[x] < 0){
                return x;
            }else{
                int parent = root(data_[x]);
                dist_[x] += dist_[data_[x]];
                return data_[x] = parent;
            }
//            return data_[x] < 0 ? x : data_[x]; //do not change the root
        }
        int size(int x) {
            return -data_[root(x)];
        }
        double distance(int x){
            root(x);
            return dist_[x];
        }

        //void make_cluster(SINGLE_MC *ht_MotifClusterinfo, SINGLE_PS &ht_PartSeqinfo);
        void make_clusterDB(const Parameter_info &param);
    };

    class Edge {
    public:
        Edge(int temp_member1, int temp_member2, double temp_cosdist)
                :member1_(temp_member1),member2_(temp_member2),cosdist_(temp_cosdist){};
        ~Edge() = default;
        int member1_, member2_;
        double cosdist_;
        bool operator<(const Edge& e) const {
            return cosdist_ < e.cosdist_;
        }
    };

    class Graph {
    public:
        ~Graph() = default;
        std::vector<Edge> edges_;

        void sort_edge() {
            sort(edges_.begin(), edges_.end());
        }
    };

private:
    
    std::string parameter_file_name_;
    

    void CalcMSD(std::vector<std::string> &mseq, const int low, const int high, std::vector<std::vector<int> > &simi, std::vector<std::string> &temp, const int n);
    void SetParameterJSON(const std::string &parameter_file, Parameter_info &param);
    void CalcMain(const Parameter_info &param);

    //DB
    void InputDatasDB(const Parameter_info &param);
    void CreateCycleSeqinfoDB(const Parameter_info &param);
    void FastaToSeq_DB(int file_index, const Parameter_info &param);
    void FastqToSeq_DB(int file_index, const Parameter_info &param);
    void AddSecondaryStructureDB(const Parameter_info &param);

    void CreateHTPartSeqinfoDB(const Parameter_info &param);
    void RunSketchsortDB(const Parameter_info &param);
    void CreateAllconnectHTMotifClusterDB(const Parameter_info &param);

    void CalcSeqFreqEnrich_DB(const Parameter_info &param);
    void CalcPartseqScoreDB(const Parameter_info &param);
    void ScoreCalculateAndInserter(const Parameter_info &param, const std::string& score_name, const std::string& query);
    void CalcClusterScoreDB(const Parameter_info &param);
    void CalcKmerScoreDB(const Parameter_info &param);
    void CalcSeqScoreDB(const Parameter_info &param);

    void ExportScoreCSV(const Parameter_info &param);
    void AddBindingFlagDB(const Parameter_info &param);
};

#endif /* RaptRanker_H */
