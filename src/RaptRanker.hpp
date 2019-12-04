//
//  RaptRanker.h
//  
//
//  Created by adinnovaiton on 2017/02/10.
//
//

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
        double alphabet_weight_;
        double cosine_distance_;
        double missing_ratio_;
        int input_file_nums_;
        //std::vector<std::string> input_file_names_;
        std::vector<Input_file_info> input_file_infos_;
        std::string experiment_dbfile_;
        std::string analysis_dbfile_;
        std::string analysis_output_path_;
        //std::string path_to_CapR;
        //std::string path_to_sketchsort;
        //std::string path_to_Weblogo;
        //int exWeblogo;


        //for dev
        bool exFiltering_ = true;
        bool exCapR_ = true;
        bool exSketchSort_ = true;
        bool exKmer_ = true;
        bool export_scoreCSV_ = true;
        bool add_binding_ = false;
        std::string bindingfile_path_;
    };

    /*struct Round_info{
    public:
        int rawdata_reads;
        int trimmed_5;
        int trimmed_3;
        int total_reads;
        int unique_reads;
    };*/

    /*class Seq_info{
    public:
        Seq_info(std::string &temp_seq, int temp_count) :seq_(std::move(temp_seq)), count_(temp_count){};
        ~Seq_info() = default;
        std::string seq_;
        int count_;
    };

    class HTSeq_info;
    class PartingSeq_info;
    class MotifCluster_info;


    class HTSeq_info{
    public:
        HTSeq_info(int temp_id, unsigned long rounds, std::string temp_seq)
                :seq_id_(temp_id),binding_(std::numeric_limits<double>::quiet_NaN()),
                 seq_(std::move(temp_seq))
        {
            std::vector<double> init(rounds,0.0);
            AveMC_profile_.reserve(rounds);
            BestMC_profile_.reserve(rounds);
            AveDegRead_profile_.reserve(rounds);
            BestDegRead_profile_.reserve(rounds);
            AveMC_profile_=init;
            BestMC_profile_=init;
            AveDegRead_profile_=init;
            BestDegRead_profile_=init;
            secondary_structure_profile_.reserve(temp_seq.size());
            alphabet_profile_.reserve(temp_seq.size());
        };
        ~HTSeq_info()=default;
        int seq_id_;
        double binding_;
        std::string seq_;
        std::vector<int> count_profile_;
        std::vector<double> AveMC_profile_;
        std::vector<double> BestMC_profile_;
        std::vector<double> AveDegRead_profile_;
        std::vector<double> BestDegRead_profile_;
        std::vector<std::vector<double> > secondary_structure_profile_;
        std::vector<std::vector<double> > alphabet_profile_;

        std::vector<PartingSeq_info*> part_ptrs_;
    };

    class PartingSeq_info{
    public:
        PartingSeq_info(int temp_ssid, int temp_pos, HTSeq_info* temp_origin)
                :ss_id_(temp_ssid), position_(temp_pos), origin_seq_(temp_origin), cluster_ptr_(nullptr){};
        ~PartingSeq_info() = default;
        int ss_id_;
        int position_;
        HTSeq_info* origin_seq_;
        MotifCluster_info* cluster_ptr_;
    };

    class MotifCluster_info{
    public:
        explicit MotifCluster_info(int temp_ss_id)
                :motif_cluster_id_(-1), root_ss_id_(temp_ss_id),
                 binding_true_(false), binding_false_(false), radius_(-1 * std::numeric_limits<double>::infinity()){};
        ~MotifCluster_info() = default;

        int motif_cluster_id_;
        int root_ss_id_;
        bool binding_true_; //true->there is
        bool binding_false_;
        double radius_;
        std::vector<double> motif_count_;

        std::vector<PartingSeq_info*> cluster_member_;

        void set_radius(double x){
            if(x > radius_) radius_ = x;
        }
    };

    typedef std::vector<Seq_info> SINGLE_SEQ;
    typedef std::vector<std::vector<Seq_info> > DOUBLE_SEQ;
    typedef std::vector<PartingSeq_info> SINGLE_PS;
    typedef std::vector<MotifCluster_info> SINGLE_MC;
    typedef std::vector<HTSeq_info> SINGLE_HT;
    typedef std::vector<Round_info> DOUBLE_ROUND;*/




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
        std::vector<Edge> edges_;  // 辺集合

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
