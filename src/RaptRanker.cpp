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

#include "RaptRanker.hpp"
#include "Stopwatch.hpp"

void RaptRanker::Run() {
    Parameter_info param;
    SetParameterJSON(parameter_file_name_, param);
    CalcMain(param);
}

void RaptRanker::SetParameterJSON(const std::string &parameter_file, Parameter_info &param) {
    std::ifstream input(parameter_file, std::ios::in);
    if (input.fail()) {
        std::cerr << "error : input parameter_file failed" << std::endl;
        exit(1);
    }
    std::string str((std::istreambuf_iterator<char>(input)),
                    std::istreambuf_iterator<char>());
    string err;
    auto json = json11::Json::parse(str, err);

    param.file_type_ = json["file_type"].int_value();
    param.forward_primer_ = json["forward_primer"].string_value();
    param.reverse_primer_ = json["reverse_primer"].string_value();
    param.add_forward_primer_ = json["add_forward_primer"].string_value();
    param.add_reverse_primer_ = json["add_reverse_primer"].string_value();
    param.sequence_maximum_length_ = json["sequence_maximum_length"].int_value();
    param.sequence_minimum_length_ = json["sequence_minimum_length"].int_value();
    param.wide_length_ = json["wide_length"].int_value();
    param.nucleotide_weight_ = json["nucleotide_weight"].number_value();
    param.cosine_distance_ = json["cosine_distance"].number_value();
    param.missing_ratio_ = json["missing_ratio"].number_value();
    param.input_file_nums_ = json["input_file_nums"].int_value();
    for (const auto &item : json["input_file_list"].array_items()) {
        Input_file_info tmp;
        tmp.round_id_ = item["round_id"].int_value();
        param.round_ids_.emplace_back(item["round_id"].int_value());
        if(item["round_id"].int_value() + 1 > param.score_stl_size_){
            param.score_stl_size_ = item["round_id"].int_value() + 1;
        }
        //convert "-" to "_" in round_name (SQL can't use "-")
        tmp.round_name_ = item["round_name"].string_value();
        std::replace(tmp.round_name_.begin(), tmp.round_name_.end(), '-', '_');
        tmp.file_path_ = item["file_path"].string_value();
        param.input_file_infos_.emplace_back(tmp);
    }
    param.experiment_dbfile_ = json["experiment_dbfile"].string_value();
    param.analysis_dbfile_ = json["analysis_dbfile"].string_value();
    param.analysis_output_path_ = json["analysis_output_path"].string_value();
    param.exKmer_ = json["exKmer"].bool_value();
    param.add_binding_ = json["add_binding"].bool_value();
    param.bindingfile_path_ = json["binding_file_path"].string_value();

    std::cout << "===== Parameter setting =====\n";
    std::cout << "file_type : " << param.file_type_ << "\n";
    std::cout << "forward_primer : " << param.forward_primer_ << "\n";
    std::cout << "reverse_primer : " << param.reverse_primer_ << "\n";
    std::cout << "add_forward_primer : " << param.add_forward_primer_ << "\n";
    std::cout << "add_reverse_primer : " << param.add_reverse_primer_ << "\n";
    std::cout << "sequence_maximum_length : " << param.sequence_maximum_length_ << "\n";
    std::cout << "sequence_minimum_length : " << param.sequence_minimum_length_ << "\n";
    std::cout << "wide_length : " << param.wide_length_ << "\n";
    std::cout << "nucleotide_weight : " << param.nucleotide_weight_ << "\n";
    std::cout << "cosine_distance : " << param.cosine_distance_ << "\n";
    std::cout << "missing_ratio : " << param.missing_ratio_ << "\n";
    std::cout << "input_file_nums : " << param.input_file_nums_ << "\n";
    std::cout << "input_files : \n";
    for (const auto &item : param.input_file_infos_) {
        std::cout << "    round_id : " << item.round_id_ << "\n";
        std::cout << "    round_name : " << item.round_name_ << "\n";
        std::cout << "    file_path : " << item.file_path_ << "\n\n";
    }
    std::cout << "experiment_dbfile : " << param.experiment_dbfile_ << "\n";
    std::cout << "analysis_output_path : " << param.analysis_output_path_ << "\n";
    std::cout << "analysis_dbfile : " << param.analysis_dbfile_ << "\n";
    std::cout << "=============================\n";
    std::cout << "===== extra options =====\n";
    std::cout << "exKmer : " << param.exKmer_ << "\n";
    std::cout << "addbinding : " << std::to_string(param.add_binding_) << "\n";
    std::cout << "bindingfile : " << param.bindingfile_path_ << "\n";
    std::cout << "=============================\n";
}

void RaptRanker::CalcMain(const Parameter_info &param) {
    Stopwatch timer(true);
    InputDatas(param);
    timer.add_check_point("DataInput");
    SecondaryStructurePrediction(param);
    timer.add_check_point("CapR");
    EnumerateSubSeqs(param);
    timer.add_check_point("MakeSubseq");
    RunSketchsort(param);
    timer.add_check_point("SketchSort");
    CalcMSF(param);
    timer.add_check_point("MakeCluster");
    CalcSeqFreqEnrich(param);
    timer.add_check_point("CalcFreqEnrich");
    CalcPartseqScore(param);
    timer.add_check_point("CalcSubseqScore");
    CalcClusterScore(param);
    timer.add_check_point("CalcClusterScore");
    if (param.exKmer_) {
        CalcKmerScore(param);
        timer.add_check_point("CalcKmerScore");
    }
    CalcSeqScoreDB(param);
    timer.add_check_point("CalcSeqScore");
    //add flags(optional)
    if (param.add_binding_) {
        AddBindingFlagDB(param);
        timer.add_check_point("AddBindFlag");
    }
    ExportScoreCSV(param);
    timer.add_check_point("ExportCSV");
    timer.finish_measurement();
}

void RaptRanker::InputDatas(const Parameter_info &param) {
    //DB preparing
    std::string dbfile = param.experiment_dbfile_;
    std::cout << "Getting ready \"" << dbfile << "\"..." << std::endl;
    try {
        SQLite::Database seqinfoDB(dbfile, SQLite::OPEN_READWRITE | SQLite::OPEN_CREATE);
        seqinfoDB.exec("PRAGMA journal_mode = WAL");
        seqinfoDB.exec(
                "CREATE TABLE IF NOT EXISTS all_seq("
                "seq_id                INTEGER  PRIMARY KEY  AUTOINCREMENT,"
                "sequence              TEXT     NOT NULL   UNIQUE)"
        );
        std::cout << "\t all_seq table is ready." << std::endl;

        seqinfoDB.exec(
                "CREATE TABLE IF NOT EXISTS binding_seq("
                "seq_id   INTEGER  NOT NULL  UNIQUE,"
                "flag     INTEGER  NOT NULL,"
                "FOREIGN KEY(seq_id) REFERENCES all_seq(seq_id)"
                ")"
        );
        std::cout << "\t binding_seq table is ready." << std::endl;

        seqinfoDB.exec(
                "CREATE TABLE IF NOT EXISTS seq_secondary_structure("
                "seq_id                INTEGER  UNIQUE,"
                "secondary_structure   TEXT     NOT NULL,"
                "FOREIGN KEY(seq_id) REFERENCES all_seq(seq_id)"
                ")"
        );
        std::cout << "\t seq_secondary_structure table is ready." << std::endl;

        seqinfoDB.exec(
                "CREATE TABLE IF NOT EXISTS round_info("
                "round_id       INTEGER  PRIMARY KEY,"
                "round_name     TEXT     NOT NULL   UNIQUE,"
                "rawdata_reads  INTEGER  ,"
                "trimmed_5      INTEGER  ,"
                "trimmed_3      INTEGER  ,"
                "total_reads    INTEGER  ,"
                "unique_reads   INTEGER  "
                ")"
        );
        std::cout << "\t round_info table is ready." << std::endl;

        std::string seq_score_template = "("
                                         "seq_id     INTEGER  NOT NULL,"
                                         "round_id   INTEGER  NOT NULL,"
                                         "value      INTEGER  NOT NULL,"
                                         "FOREIGN KEY(seq_id) REFERENCES all_seq(seq_id),"
                                         "UNIQUE(seq_id, round_id)"
                                         ")";
        seqinfoDB.exec("CREATE TABLE IF NOT EXISTS seq_count" + seq_score_template);
        seqinfoDB.exec("CREATE TABLE IF NOT EXISTS seq_freq" + seq_score_template);
        seqinfoDB.exec("CREATE TABLE IF NOT EXISTS seq_enrich" + seq_score_template);
        seqinfoDB.exec("CREATE INDEX IF NOT EXISTS seq_freq_value_index ON seq_freq(value)");
        seqinfoDB.exec("CREATE INDEX IF NOT EXISTS seq_count_value_index ON seq_count(value)");
        seqinfoDB.exec("CREATE INDEX IF NOT EXISTS seq_enrich_value_index ON seq_enrich(value)");
        std::cout << "\t seq score tables and indexes are ready." << std::endl;
    } catch (std::exception &e) {
        std::cerr << "Error in preparing experiment DB file." << std::endl;
        std::cerr << "exception:" << e.what() << std::endl;
        exit(1);
    }
    std::cout << "\"" << dbfile << "\" is ready." << std::endl;

    //Importing
    for (int file_index = 0; file_index < param.input_file_nums_; ++file_index) {
        //open DB
        std::string dbfile = param.experiment_dbfile_;

        int seq_count = 0;
        try {
            SQLite::Database seqinfoDB(dbfile, SQLite::OPEN_READWRITE);
            SQLite::Statement select_query(seqinfoDB, "SELECT COUNT(*) FROM round_info WHERE round_id = ?1");
            select_query.bind(1, param.input_file_infos_[file_index].round_id_);
            SQLite::Transaction transaction(seqinfoDB);
            //get seq_id and sequence not secondary-structure predicted yet.
            while (select_query.executeStep()) {
                seq_count = select_query.getColumn(0).getInt();
            }
            transaction.commit();
        } catch (std::exception &e) {
            std::cerr << "exception:" << e.what() << std::endl;
            exit(1);
        }

        if (seq_count == 0) {
            std::cout << "Start importing round_id=" << param.input_file_infos_[file_index].round_id_
                      << " (\"" << param.input_file_infos_[file_index].file_path_ << "\")..." << std::endl;
            switch (param.file_type_) {
                case 1:
                    ImportFasta(file_index, param);
                    break;
                case 2:
                    ImportFastq(file_index, param);
                    break;
                default:
                    std::cerr << "file_type is invalid. file_type must be 1 (FASTA) or 2 (FASTQ)." << std::endl;
                    exit(1);
            }
        } else {
            int rawdata = 0;
            int filtered = 0;
            int unique = 0;
            try {
                SQLite::Database seqinfoDB(dbfile, SQLite::OPEN_READWRITE);
                SQLite::Statement select_query(seqinfoDB, "SELECT rawdata_reads,total_reads,unique_reads FROM round_info WHERE round_id = ?1");
                select_query.bind(1, param.input_file_infos_[file_index].round_id_);
                SQLite::Transaction transaction(seqinfoDB);
                //get seq_id and sequence not secondary-structure predicted yet.
                while (select_query.executeStep()) {
                    rawdata = select_query.getColumn(0).getInt();
                    filtered = select_query.getColumn(1).getInt();
                    unique = select_query.getColumn(2).getInt();
                }
                transaction.commit();
            } catch (std::exception &e) {
                std::cerr << "exception:" << e.what() << std::endl;
                exit(1);
            }
            std::cout << "round_id=" << param.input_file_infos_[file_index].round_id_
                      << " (\"" << param.input_file_infos_[file_index].file_path_ << "\") is already imported."
                      << "\n\t inputted : " << rawdata << "\n\t filtered : " << filtered << "\n\t unique : " << unique
                      << std::endl;
        }
    }
}

void
RaptRanker::CalcMSD(std::vector<std::string> &mseq, const int low, const int high, std::vector<std::vector<int> > &simi,
                    std::vector<std::string> &temp, const int n = 0) {
    int cal_length = high - low + 1;
    int count[key_length];
    for (int i = 0; i < key_length; ++i) {
        count[i] = 0;
    }
    int *trace;
    trace = new int[cal_length]();
    for (int x = low; x <= high; x++) {
        for (int y = 0; y < key_length; y++) {
            if (mseq[x][n] == key[y]) {
                count[y] += 1;
                trace[x - low] = y;
                break;
            } else {
                if (y == key_length - 1) {
                    count[key_length - 1] += 1;
                    trace[x - low] = key_length - 1;
                }
            }
        }
    }
    for (int x = 1; x < key_length; x++) {
        count[x] += count[x - 1];
    }
    int END = count[key_length - 1] - count[key_length - 2];
    if (END > 1) {
        int non_low = low + count[key_length - 1] - END;
        int non_high = low + count[key_length - 1] - 1;
        std::vector<int> add;
        add.emplace_back(non_low);
        add.emplace_back(non_high - non_low + 1);
        simi.emplace_back(add);
        add.clear();
    }
    for (int x = cal_length - 1; x >= 0; x--) {
        int c = trace[x];
        int k = count[c];
        temp[low + k - 1] = mseq[low + x];
        count[c] -= 1;
    }
    delete[] trace;
    for (int x = 0; x < cal_length; x++) {
        mseq[low + x] = temp[low + x];
    }

    for (int x = 0; x < key_length - 1; x++) {
        int L = low + count[x];
        int H = low + count[x + 1] - 1;
        if (L < H) {
            CalcMSD(mseq, L, H, simi, temp, n + 1);
        }
    }
}

void RaptRanker::ImportFasta(int file_index, const Parameter_info &param) {
    std::ifstream input(param.input_file_infos_[file_index].file_path_);
    std::string line;
    int rawdata_reads = 0;
    int trimmed_5 = 0;
    int trimmed_3 = 0;
    int total_reads = 0;
    int round_id = param.input_file_infos_[file_index].round_id_;
    if (input.fail()) {
        std::cerr << "input " + param.input_file_infos_[file_index].file_path_ << " failed." << std::endl;
        exit(1);
    }
    std::string dbfile = param.experiment_dbfile_;

    try {
        SQLite::Database seqinfoDB(dbfile, SQLite::OPEN_READWRITE);
        //Reset the trigger of SQL if exists.
        std::string temp = "DROP TRIGGER IF EXISTS seq_count_insert";
        seqinfoDB.exec(temp);

        temp = "CREATE TRIGGER IF NOT EXISTS seq_count_insert AFTER INSERT ON all_seq "
               "BEGIN "
               "INSERT INTO seq_count(seq_id, round_id, value) VALUES(NEW.ROWID,"
               + std::to_string(round_id) + ", 1); "
                                            "END";
        seqinfoDB.exec(temp);

        SQLite::Statement select_query(seqinfoDB, "SELECT seq_id FROM all_seq WHERE sequence = ?1");
        SQLite::Statement allseq_insert_query(seqinfoDB, "INSERT INTO all_seq(sequence) VALUES(?1)");
        temp = "INSERT OR IGNORE INTO seq_count(seq_id, round_id, value) VALUES(?1,"
               + std::to_string(round_id) + ",0) ";
        SQLite::Statement seq_count_insert_query(seqinfoDB, temp);
        temp = "UPDATE seq_count SET value = value + 1 "
               "WHERE(seq_id == ?1 AND round_id == " + std::to_string(round_id) + ")";
        SQLite::Statement update_query(seqinfoDB, temp);

        SQLite::Transaction transaction(seqinfoDB);
        std::string sequence;
        while (getline(input, line)) {

            if (line.empty()) {
                continue;
            }

            if (line[0] == '>') {
                if (!sequence.empty()) {
                    ++rawdata_reads;

                    if (sequence.rfind(param.reverse_primer_) != std::string::npos) {
                        int rr = (int) sequence.rfind(param.reverse_primer_);
                        ++trimmed_3;
                        if (sequence.find(param.forward_primer_) != std::string::npos) {
                            int ff = (int) sequence.find(param.forward_primer_) + (int) param.forward_primer_.size();
                            ++trimmed_5;
                            int random_length = rr - ff;
                            if (param.sequence_minimum_length_ <= random_length &&
                                random_length <= param.sequence_maximum_length_) {
                                ++total_reads;
                                std::string add = sequence.substr((unsigned long) ff, (unsigned long) random_length);

                                select_query.bind(1, add);
                                int seq_id = 0;
                                while (select_query.executeStep()) {
                                    seq_id = select_query.getColumn(0).getInt();
                                }
                                select_query.reset();

                                if (seq_id == 0) {
                                    //the sequence is found first time.
                                    allseq_insert_query.bind(1, add);
                                    while (allseq_insert_query.executeStep()) {}
                                    allseq_insert_query.reset();
                                } else {
                                    //the sequence is found more than two times.
                                    seq_count_insert_query.bind(1, seq_id);
                                    while (seq_count_insert_query.executeStep()) {}
                                    seq_count_insert_query.reset();

                                    update_query.bind(1, seq_id);
                                    while (update_query.executeStep()) {}
                                    update_query.reset();
                                }
                            }
                        }
                    }
                }//if(!sequence.empty)
                sequence.clear();
            } else {
                sequence += line;
            }
        }//while(getline)
        transaction.commit();

        temp = "SELECT COUNT(*) FROM seq_count WHERE round_id == " + std::to_string(round_id);
        SQLite::Statement unique_count_query(seqinfoDB, temp);
        int unique_reads = -1;
        while (unique_count_query.executeStep()) {
            unique_reads = unique_count_query.getColumn(0).getInt();
        }

        temp = "INSERT INTO round_info(round_id, round_name, rawdata_reads, trimmed_5, trimmed_3, total_reads, unique_reads) "
               "VALUES(" + std::to_string(round_id)
               + ",\'" + param.input_file_infos_[file_index].round_name_ + "\',"
               + std::to_string(rawdata_reads) + ","
               + std::to_string(trimmed_5) + ","
               + std::to_string(trimmed_3) + ","
               + std::to_string(total_reads) + ","
               + std::to_string(unique_reads)
               + ")";
        seqinfoDB.exec(temp);

        if(total_reads==0 || unique_reads==0){
            std::cerr << "Error: RaptRanker could't find any sequence to be analysed.\n"
                         "Please check 'forward_primer' and 'reverse_primer' parameters are correct."
                      << std::endl;
            exit(1);
        }
        std::cout << "\t imported successfully."
                  << "\n\t inputted : " << rawdata_reads
                  << "\n\t filtered : " << total_reads
                  << "\n\t unique : " << unique_reads
                  << std::endl;
    } catch (std::exception &e) {
        std::cerr << "Error in importing fasta" << std::endl;
        std::cerr << "exception:" << e.what() << std::endl;
        exit(1);
    }
}

// FASTQ → DB
void RaptRanker::ImportFastq(int file_index, const Parameter_info &param) {
    std::ifstream input(param.input_file_infos_[file_index].file_path_);
    std::string line;
    int rawdata_reads = 0;
    int trimmed_5 = 0;
    int trimmed_3 = 0;
    int total_reads = 0;
    int round_id = param.input_file_infos_[file_index].round_id_;
    if (input.fail()) {
        std::cerr << "input " << param.input_file_infos_[file_index].file_path_ << " failed." << std::endl;
        exit(1);
    }
    std::string dbfile = param.experiment_dbfile_;

    try {
        SQLite::Database seqinfoDB(dbfile, SQLite::OPEN_READWRITE);

        std::string temp = "DROP TRIGGER IF EXISTS seq_count_insert";
        seqinfoDB.exec(temp);

        temp = "CREATE TRIGGER IF NOT EXISTS seq_count_insert AFTER INSERT ON all_seq "
               "BEGIN "
               "INSERT INTO seq_count(seq_id, round_id, value) VALUES(NEW.ROWID,"
               + std::to_string(round_id) + ", 1); "
                                            "END";
        seqinfoDB.exec(temp);

        SQLite::Statement select_query(seqinfoDB, "SELECT seq_id FROM all_seq WHERE sequence = ?1");
        SQLite::Statement allseq_insert_query(seqinfoDB, "INSERT INTO all_seq(sequence) VALUES(?1)");
        temp = "INSERT OR IGNORE INTO seq_count(seq_id, round_id, value) VALUES(?1,"
               + std::to_string(round_id) + ",0) ";
        SQLite::Statement seq_count_insert_query(seqinfoDB, temp);
        temp = "UPDATE seq_count SET value = value + 1 "
               "WHERE(seq_id == ?1 AND round_id == " + std::to_string(round_id) + ")";
        SQLite::Statement update_query(seqinfoDB, temp);

        SQLite::Transaction transaction(seqinfoDB);
        while (getline(input, line)) {
            getline(input, line);
            ++rawdata_reads;
            if (line.rfind(param.reverse_primer_) != std::string::npos) {
                int rr = (int) line.rfind(param.reverse_primer_);
                ++trimmed_3;
                if (line.find(param.forward_primer_) != std::string::npos) {
                    int ff = (int) line.find(param.forward_primer_) + (int) param.forward_primer_.size();
                    ++trimmed_5;
                    int random_length = rr - ff;
                    if (param.sequence_minimum_length_ <= random_length &&
                        random_length <= param.sequence_maximum_length_) {
                        ++total_reads;
                        std::string add = line.substr((unsigned long) ff, (unsigned long) random_length);

                        select_query.bind(1, add);
                        int seq_id = 0;
                        while (select_query.executeStep()) {
                            seq_id = select_query.getColumn(0).getInt();
                        }
                        select_query.reset();

                        if (seq_id == 0) {
                            //the sequence is found first time.
                            allseq_insert_query.bind(1, add);
                            while (allseq_insert_query.executeStep()) {}
                            allseq_insert_query.reset();
                        } else {
                            //the sequence is found more than two times.
                            seq_count_insert_query.bind(1, seq_id);
                            while (seq_count_insert_query.executeStep()) {}
                            seq_count_insert_query.reset();

                            update_query.bind(1, seq_id);
                            while (update_query.executeStep()) {}
                            update_query.reset();
                        }
                    }
                }
            }
            getline(input, line);
            getline(input, line);
        }
        transaction.commit();

        temp = "SELECT COUNT(*) FROM seq_count WHERE round_id == " + std::to_string(round_id);
        SQLite::Statement unique_count_query(seqinfoDB, temp);
        int unique_reads = -1;
        while (unique_count_query.executeStep()) {
            unique_reads = unique_count_query.getColumn(0).getInt();
        }
        temp = "INSERT INTO round_info(round_id, round_name, rawdata_reads, trimmed_5, trimmed_3, total_reads, unique_reads) "
               "VALUES(" + std::to_string(round_id)
               + ",\'" + param.input_file_infos_[file_index].round_name_ + "\',"
               + std::to_string(rawdata_reads) + ","
               + std::to_string(trimmed_5) + ","
               + std::to_string(trimmed_3) + ","
               + std::to_string(total_reads) + ","
               + std::to_string(unique_reads)
               + ")";
        seqinfoDB.exec(temp);

        if(total_reads==0 || unique_reads==0){
            std::cerr << "Error: RaptRanker could't find any sequence to be analysed.\n"
                         "Please check 'forward_primer' and 'reverse_primer' parameters are correct."
                      << std::endl;
            exit(1);
        }
        std::cout << "\t imported successfully."
                  << "\n\t inputted : " << rawdata_reads
                  << "\n\t filtered : " << total_reads
                  << "\n\t unique : " << unique_reads
                  << std::endl;
    } catch (std::exception &e) {
        std::cerr << "Error in importing fastq" << std::endl;
        std::cerr << "exception:" << e.what() << std::endl;
        exit(1);
    }
}

std::vector<std::string> split(const std::string &str, const char delim) {
    std::vector<std::string> res;
    size_t current = 0, found;
    while ((found = str.find_first_of(delim, current)) != std::string::npos) {
        res.emplace_back(std::string(str, current, found - current));
        current = found + 1;
    }
    res.emplace_back(std::string(str, current, str.size() - current));
    return res;
}
bool file_exists(const std::string &file_path) {
    std::ifstream ifs(file_path);
    return ifs.is_open();
}

void RaptRanker::SecondaryStructurePrediction(const Parameter_info &param) {
    std::cout << "Predicting secondary structures... \n"
                 "(This process may take a while. Even if it looks like it's freezed, please waiting for a while.)"
                 << std::endl;
    std::string dbfile = param.experiment_dbfile_;
    try {
        SQLite::Database seqinfoDB(dbfile, SQLite::OPEN_READWRITE);
        SQLite::Statement select_query(seqinfoDB, "SELECT * FROM all_seq LEFT OUTER JOIN seq_secondary_structure "
                                                  "ON (all_seq.seq_id = seq_secondary_structure.seq_id)"
                                                  "WHERE secondary_structure IS NULL");
        SQLite::Statement insert_query(seqinfoDB,
                                       "INSERT INTO seq_secondary_structure(seq_id,secondary_structure) VALUES(?1,?2)");
        SQLite::Transaction transaction(seqinfoDB);

        while (select_query.executeStep()) {
            //get seq_id and sequence not secondary-structure predicted yet.
            int seq_id = select_query.getColumn(0).getInt();
            std::string seq =
                    param.add_forward_primer_ + select_query.getColumn(1).getText() + param.add_reverse_primer_;
            //secondary-structure prediction by CapR
            CapR capr(200);
            std::string result = capr.Run_for_RaptRanker((int) param.add_forward_primer_.size(),
                                                         (int) param.add_reverse_primer_.size(), seq);
            //update results to DB
            insert_query.bind(1, seq_id);
            insert_query.bind(2, result);
            while (insert_query.executeStep()) {}
            insert_query.reset();
        }
        transaction.commit();
    } catch (std::exception &e){
        std::cerr << "Error in secondary structure prediction." << std::endl;
        std::cerr << "exception:" << e.what() << std::endl;
        exit(1);
    }
    std::cout << "Finished secondary structure prediction." << std::endl;
}

void RaptRanker::EnumerateSubSeqs(const Parameter_info &param) {
    std::cout << "Enumerating subsequences..." << std::endl;
    std::string dbfile = param.analysis_dbfile_;
    if (file_exists(dbfile)) {
        std::cerr << "\"" << dbfile << "\" is already exists.\n";
        exit(1);
    }
    try {
        SQLite::Database resultDB(dbfile, SQLite::OPEN_READWRITE | SQLite::OPEN_CREATE);
        resultDB.exec("PRAGMA journal_mode = WAL");
        resultDB.exec("CREATE TABLE IF NOT EXISTS all_subseq("
                      "subseq_id       INTEGER  PRIMARY KEY,"
                      "sub_sequence    TEXT     NOT NULL,"
                      "sub_structure   TEXT     NOT NULL,"
                      "seq_id           INTEGER  NOT NULL,"
                      "position         INTEGER  NOT NULL,"
                      "UNIQUE(seq_id,position)"
                      ")");
        std::cout << "\t all_subseq Table is ready." << std::endl;

        resultDB.exec("CREATE TABLE IF NOT EXISTS subseq_cluster("
                      "subseq_id  INTEGER  NOT NULL  UNIQUE,"
                      "cluster_id  INTEGER  NOT NULL,"
                      "FOREIGN KEY(subseq_id) REFERENCES all_subseq(subseq_id)"
                      ")");
        std::cout << "\t subseq_cluster Table is ready." << std::endl;

        std::string sub_score_template = "("
                                         "subseq_id  INTEGER  NOT NULL,"
                                         "round_id    INTEGER  NOT NULL,"
                                         "value       NUMERIC  NOT NULL,"
                                         "FOREIGN KEY(subseq_id) REFERENCES all_subseq(subseq_id),"
                                         "UNIQUE(subseq_id, round_id)"
                                         ")";
        resultDB.exec("CREATE TABLE IF NOT EXISTS subseq_count" + sub_score_template);
        resultDB.exec("CREATE TABLE IF NOT EXISTS subseq_freq" + sub_score_template);
        std::cout << "\t sub_sequence score Tables are ready." << std::endl;

        std::string cluster_score_template = "("
                                             "cluster_id  INTEGER  NOT NULL,"
                                             "round_id    INTEGER  NOT NULL,"
                                             "value       NUMERIC  NOT NULL,"
                                             "FOREIGN KEY(cluster_id) REFERENCES subseq_cluster(cluster_id),"
                                             "UNIQUE(cluster_id, round_id)"
                                             ")";
        resultDB.exec("CREATE TABLE IF NOT EXISTS cluster_count" + cluster_score_template);
        resultDB.exec("CREATE TABLE IF NOT EXISTS cluster_deg" + cluster_score_template);
        resultDB.exec("CREATE TABLE IF NOT EXISTS cluster_freq" + cluster_score_template);
        resultDB.exec("CREATE TABLE IF NOT EXISTS cluster_enrich" + cluster_score_template);
        std::cout << "\t cluster score Tables are ready." << std::endl;

        resultDB.exec(
                "CREATE TABLE IF NOT EXISTS seq_AMF("
                "seq_id      INTEGER  NOT NULL,"
                "round_id    INTEGER  NOT NULL,"
                "value       NUMERIC  NOT NULL,"
                "UNIQUE(seq_id, round_id)"
                ")"
        );
        resultDB.exec(
                "CREATE TABLE IF NOT EXISTS seq_AME("
                "seq_id      INTEGER  NOT NULL,"
                "round_id    INTEGER  NOT NULL,"
                "value       NUMERIC  NOT NULL,"
                "UNIQUE(seq_id, round_id)"
                ")"
        );
        std::cout << "\t sequence score Tables are ready." << std::endl;

        resultDB.exec("CREATE INDEX IF NOT EXISTS subsequence_index ON all_subseq(sub_sequence)");
        resultDB.exec("CREATE INDEX IF NOT EXISTS clusterid_index ON subseq_cluster(cluster_id)");
        resultDB.exec("CREATE INDEX IF NOT EXISTS subseq_freq_value_index ON subseq_freq(value)");
        resultDB.exec("CREATE INDEX IF NOT EXISTS subseq_count_value_index ON subseq_count(value)");
        resultDB.exec("CREATE INDEX IF NOT EXISTS cluster_freq_value_index ON cluster_freq(value)");
        resultDB.exec("CREATE INDEX IF NOT EXISTS cluster_deg_value_index ON cluster_deg(value)");
        resultDB.exec("CREATE INDEX IF NOT EXISTS cluster_count_value_index ON cluster_count(value)");
        resultDB.exec("CREATE INDEX IF NOT EXISTS cluster_enrich_value_index ON cluster_enrich(value)");
        std::cout << "\t indexes are ready." << std::endl;
    } catch (std::exception &e) {
        std::cerr << "Error in preparing analysis DB file." << std::endl;
        std::cerr << "exception:" << e.what() << std::endl;
        exit(1);
    }

    try {
        dbfile = param.analysis_dbfile_;
        SQLite::Database resultDB(dbfile, SQLite::OPEN_READWRITE);
        dbfile = param.experiment_dbfile_;
        SQLite::Database seqinfoDB(dbfile, SQLite::OPEN_READWRITE);
        SQLite::Statement select_query(seqinfoDB,
                                       "SELECT all_seq.seq_id,sequence,seq_secondary_structure.secondary_structure"
                                       " FROM all_seq NATURAL INNER JOIN seq_secondary_structure");
        SQLite::Statement insert_query(resultDB,
                                       "INSERT INTO all_subseq(seq_id, sub_sequence, sub_structure, position) VALUES(?1, ?2, ?3, ?4)");
        SQLite::Transaction seq_transaction(seqinfoDB);
        SQLite::Transaction result_transaction(resultDB);
        while (select_query.executeStep()) {
            int seq_id = select_query.getColumn(0).getInt();
            std::string sequence = select_query.getColumn(1).getText();
            std::vector<std::string> structure;
            std::string temp_str = select_query.getColumn(2).getText();
            std::stringstream input_sstream(temp_str);
            std::string line;
            while (getline(input_sstream, line)) {
                structure.emplace_back(line);
            }
            unsigned long position = 0;
            unsigned long last = sequence.length() - param.wide_length_;
            
            while (position <= last) {
                std::string subseq = sequence.substr(position, (unsigned long) param.wide_length_);
                std::string sub_structure;
                for (unsigned long i = position, n = position + param.wide_length_; i < n; ++i) {
                    sub_structure += structure[i] + "\n";
                }

                insert_query.bind(1, seq_id);
                insert_query.bind(2, subseq);
                insert_query.bind(3, sub_structure);
                insert_query.bind(4, (int) position + 1);

                while (insert_query.executeStep()) {}
                insert_query.reset();

                ++position;
            }
        }//while(select_query)

        seq_transaction.commit();
        result_transaction.commit();
    } catch (std::exception &e) {
        std::cerr << "Error in enumerating subsequences." << std::endl;
        std::cerr << "exception:" << e.what() << std::endl;
        exit(1);
    }
    std::cout << "Finished subsequences enumeration." << std::endl;
}

void RaptRanker::RunSketchsort(const Parameter_info &param) {
    std::cout << "Running SketchSort..." << std::endl;
    std::string dbfile = param.experiment_dbfile_;
    std::string inputSS = param.analysis_output_path_ + "inputSS_AllRound.txt"; //input file for SketchSort
    std::ofstream output(inputSS);
    try {
        SQLite::Database seqinfoDB(dbfile, SQLite::OPEN_READWRITE);
        SQLite::Statement select_query(seqinfoDB, "SELECT sequence,seq_secondary_structure.secondary_structure"
                                                  " FROM all_seq NATURAL INNER JOIN seq_secondary_structure");
        SQLite::Transaction transaction(seqinfoDB);
        while (select_query.executeStep()) {
            std::string sequence = select_query.getColumn(0).getText();
            std::vector<std::string> structure;
            std::string temp_str = select_query.getColumn(1).getText();
            std::stringstream input_sstream(temp_str);
            std::string line;
            while (getline(input_sstream, line)) {
                structure.emplace_back(line);
            }
            //each subseq
            unsigned long position = 0;
            unsigned long last = sequence.length() - param.wide_length_;
            while (position <= last) {
                std::string subseq = sequence.substr(position, (unsigned long) param.wide_length_);
                std::ostringstream output_sstream;
                for (int j = 0; j < param.wide_length_; ++j) {
                    //sequence profile
                    if (subseq[j] == 'A') {
                        output_sstream << param.nucleotide_weight_ << " 0 0 0 ";
                    } else if (subseq[j] == 'G') {
                        output_sstream << "0 " << param.nucleotide_weight_ << " 0 0 ";
                    } else if (subseq[j] == 'C') {
                        output_sstream << "0 0 " << param.nucleotide_weight_ << " 0 ";
                    } else if (subseq[j] == 'T' || subseq[j] == 'U') {
                        output_sstream << "0 0 0 " << param.nucleotide_weight_ << " ";
                    } else {
                        output_sstream << "0 0 0 0 ";
                        std::cerr << "RaptRanker only considers ACG(T/U). Other strings like 'N' may cause errors.";
                    }
                    //secondary structure profile
                    output_sstream << structure[position + j] << " ";
                }
                ++position;
                output << output_sstream.str() << "\n";
            }
        }//select_query

        output.close();
        transaction.commit();
    } catch (std::exception &e) {
        std::cerr << "Error in making input file for SketchSort" << std::endl;
        std::cerr << "exception:" << e.what() << std::endl;
        exit(1);
    }

    //excute SketchSort
    std::string outputSS = param.analysis_output_path_ + "outputSS_AllRound.txt"; //output file of SketchSort
    SketchSort sketchsort;
    sketchsort.run(inputSS.c_str(), outputSS.c_str(), 4, 1, (float) param.cosine_distance_, 3, true,
                   (float) param.missing_ratio_, true);
}

void RaptRanker::CalcMSF(const Parameter_info &param) {
    std::cout << "Clustering with Minimum Spanning Forest..." << std::endl;
    auto graph = new Graph;

    //input sketchsort_similar_set
    std::string outputSS = param.analysis_output_path_ + "outputSS_AllRound.txt"; //output file of SketchSort
    std::ifstream input(outputSS);
    std::string line;
    if (input.fail()) {
        std::cerr << "Error: failed at opening output of SketchSort." << std::endl;
    }
    while (getline(input, line)) {
        std::vector<std::string> split_line = split(line, ' ');
        graph->edges_.emplace_back(std::stoi(split_line[0]), std::stoi(split_line[1]), 1 - std::stod(split_line[2]));
        // member1, member2, cosdist = 2->Reverse　0->same
    }
    input.close();

    std::string dbfile = param.analysis_dbfile_;
    int records_num = 0; //number of subseq records
    try {
        SQLite::Database seqinfoDB(dbfile, SQLite::OPEN_READWRITE);
        SQLite::Transaction transaction(seqinfoDB);
        records_num = seqinfoDB.execAndGet("SELECT MAX(subseq_id) FROM all_subseq");
        transaction.commit();
    } catch (std::exception &e) {
        std::cerr << "Error in get MAX(subseq_id) for making UnionFind" << std::endl;
        std::cerr << "exception:" << e.what() << std::endl;
        exit(1);
    }

    WeightedUnionFind uf((unsigned long) records_num);
    //make Minimum Spanning Tree(s)  (minimum spanning forest)
    graph->sort_edge();
    for (const auto &each_edge:graph->edges_) {
        uf.unite(each_edge.member1_, each_edge.member2_, each_edge.cosdist_);
    }
    delete (graph);

    //update all distances
    for (int i = 0, n = records_num; i < n; ++i) {
        uf.distance(i);
    }
    uf.make_clusterDB(param);

    std::cout << "Finished clustering." << std::endl;
}

void RaptRanker::WeightedUnionFind::make_clusterDB(const Parameter_info &param) {
    //openDB
    std::string dbfile = param.analysis_dbfile_;
    try {
        SQLite::Database resultDB(dbfile, SQLite::OPEN_READWRITE);
        SQLite::Statement insert_query(resultDB,
                                       "INSERT OR IGNORE INTO subseq_cluster(subseq_id,cluster_id) VALUES(?1,?2)");
        SQLite::Transaction transaction(resultDB);

        //make clusters (define cluster_id for each root subseqs)
        int cluster_id = 1;
        for (int i = 0, n = (int) data_.size(); i < n; ++i) {
            if (data_[i] < 0) {
                cluster_ids_[i] = cluster_id;
                ++cluster_id;
            }
        }

        //add members
        for (int i = 0, n = (int) data_.size(); i < n; ++i) {
            cluster_id = cluster_ids_[root(i)];

            insert_query.bind(1, i + 1);
            insert_query.bind(2, cluster_id);
            while (insert_query.executeStep()) {}
            insert_query.reset();
        }

        transaction.commit();
    } catch (std::exception &e) {
        std::cerr << "Error in inserting into subseq_cluster table" << std::endl;
        std::cerr << "exception:" << e.what() << std::endl;
        exit(1);
    }
}

void RaptRanker::CalcSeqFreqEnrich(const Parameter_info &param) {
    std::cout << "Calculating Frequency and Enrichment..." << std::endl;
    try {
        std::string dbfile = param.experiment_dbfile_;
        SQLite::Database seqinfoDB(dbfile, SQLite::OPEN_READWRITE);

        SQLite::Statement id_select_query(seqinfoDB, "SELECT DISTINCT(all_seq.seq_id) FROM all_seq "
                                                     "LEFT OUTER JOIN seq_freq ON (all_seq.seq_id = seq_freq.seq_id)"
                                                     "WHERE value IS NULL");
//        SQLite::Statement max_round_id_select_query(seqinfoDB, "SELECT MAX(round_id) FROM round_info ");
        SQLite::Statement totalread_select_query(seqinfoDB, "SELECT round_id,total_reads FROM round_info");
        SQLite::Statement counts_select_query(seqinfoDB, "SELECT round_id, value FROM seq_count WHERE seq_id == ?1");
        SQLite::Statement freq_insert_query(seqinfoDB,
                                            "INSERT INTO seq_freq(seq_id,round_id,value) VALUES(?1, ?2, ?3)");
        SQLite::Statement enrich_insert_query(seqinfoDB,
                                              "INSERT INTO seq_enrich(seq_id,round_id,value) VALUES(?1, ?2, ?3)");

        SQLite::Transaction seqinfo_transaction(seqinfoDB);
        //get total_reads of each round
//        int size = -1;
//        while (max_round_id_select_query.executeStep()) {
//            size = max_round_id_select_query.getColumn(0).getInt() + 1;
//        }
//        max_round_id_select_query.reset();

//        std::vector<int> round_ids;
        std::vector<double> total_reads((unsigned long) param.score_stl_size_, 0.0);
        while (totalread_select_query.executeStep()) {
            int round = totalread_select_query.getColumn(0).getInt();
//            param.round_ids.emplace_back(round);
            total_reads[round] = totalread_select_query.getColumn(1).getDouble();
        }
        totalread_select_query.reset();

        //for each seq
        while (id_select_query.executeStep()) {
            int seq_id = id_select_query.getColumn(0).getInt();
            std::vector<double> Frequency((unsigned long) param.score_stl_size_, 0.0);

            //calc Frequency
            counts_select_query.bind(1, seq_id);
            while (counts_select_query.executeStep()) {
                int round = counts_select_query.getColumn(0).getInt();
                double count = counts_select_query.getColumn(1).getDouble();
                Frequency[round] = count / total_reads[round];
            }
            counts_select_query.reset();

            //insert Frequency
            for (const auto & round : param.round_ids_) {
                freq_insert_query.bind(1, seq_id);
                freq_insert_query.bind(2, round);
                freq_insert_query.bind(3, Frequency[round]);
                while (freq_insert_query.executeStep()) {}
                freq_insert_query.reset();
            }//for round

            //enrich scores
            for (int index = 1, m = param.round_ids_.size(); index < m; ++index) {
                int round = param.round_ids_[index];
                int prev_round = param.round_ids_[index - 1];

                if (Frequency[prev_round] != 0) {
                    enrich_insert_query.bind(1, seq_id);
                    enrich_insert_query.bind(2, round);
                    enrich_insert_query.bind(3, Frequency[round] / Frequency[prev_round]);
                    while (enrich_insert_query.executeStep()) {}
                    enrich_insert_query.reset();
                }
            }//for round
        }//for seq
        id_select_query.reset();
        seqinfo_transaction.commit();
    } catch (std::exception &e) {
        std::cerr << "exception:" << e.what() << std::endl;
        exit(1);
    }
    std::cout << "Finished Frequency and Enrichment calculation." << std::endl;
}

void RaptRanker::CalcPartseqScore(const Parameter_info &param) {
    std::cout << "Calculating subsequence scores..." << std::endl;

    //subseq count
    try {
        std::string dbfile = param.analysis_dbfile_;
        SQLite::Database resultDB(dbfile, SQLite::OPEN_READWRITE);
        dbfile = param.experiment_dbfile_;
        SQLite::Database seqinfoDB(dbfile, SQLite::OPEN_READWRITE);
        SQLite::Statement seq_count_select_query(seqinfoDB, "SELECT round_id, value FROM seq_count WHERE seq_id == ?1");
        SQLite::Statement subseq_id_select_query(resultDB, "SELECT subseq_id FROM all_subseq WHERE seq_id == ?1");
        SQLite::Statement subseq_count_insert_query(resultDB,
                                                    "INSERT INTO subseq_count(subseq_id,round_id,value) VALUES(?1,?2,?3)");
        SQLite::Transaction seqinfo_transaction(seqinfoDB);
        SQLite::Transaction result_transaction(resultDB);

        int last_seq_id = seqinfoDB.execAndGet("SELECT MAX(seq_id) FROM all_seq");
        for (int seq_id = 1; seq_id <= last_seq_id; ++seq_id) {
            std::vector<int> counts((unsigned long) param.score_stl_size_, 0);
            seq_count_select_query.bind(1, seq_id);

            while (seq_count_select_query.executeStep()) {
                int round = seq_count_select_query.getColumn(0).getInt();
                counts[round] = seq_count_select_query.getColumn(1).getInt();
            }
            seq_count_select_query.reset();

            subseq_id_select_query.bind(1, seq_id);
            while (subseq_id_select_query.executeStep()) {
                int sub_id = subseq_id_select_query.getColumn(0).getInt();
                for (const auto & round : param.round_ids_) {
                    subseq_count_insert_query.bind(1, sub_id);
                    subseq_count_insert_query.bind(2, round);
                    subseq_count_insert_query.bind(3, counts[round]);
                    while (subseq_count_insert_query.executeStep()) {}
                    subseq_count_insert_query.reset();
                }
            }
            subseq_id_select_query.reset();
        }//for seq
        seqinfo_transaction.commit();
        result_transaction.commit();
    } catch (std::exception &e) {
        std::cerr << "Error in calculation and inserting of subseq count" << std::endl;
        std::cerr << "exception:" << e.what() << std::endl;
        exit(1);
    }

    //subseq freq
    try {
        std::string dbfile = param.analysis_dbfile_;
        SQLite::Database resultDB(dbfile, SQLite::OPEN_READWRITE);
        SQLite::Statement total_count_select_query(resultDB,
                                                   "SELECT round_id,SUM(value) FROM subseq_count GROUP BY round_id");
        SQLite::Statement subseq_count_select_query(resultDB,
                                                    "SELECT round_id,value FROM subseq_count WHERE subseq_id == ?1");
        SQLite::Statement subseq_freq_insert_query(resultDB,
                                                   "INSERT INTO subseq_freq(subseq_id,round_id,value) VALUES(?1,?2,?3)");
        SQLite::Transaction result_transaction(resultDB);

        std::vector<double> total_counts((unsigned long) param.score_stl_size_, 0.0);
        while (total_count_select_query.executeStep()) {
            int round = total_count_select_query.getColumn(0).getInt();
            total_counts[round] = total_count_select_query.getColumn(1).getDouble();
        }
        total_count_select_query.reset();

        int last_subseq_id = resultDB.execAndGet("SELECT MAX(subseq_id) FROM all_subseq");
        for (int subseq_id = 1; subseq_id <= last_subseq_id; ++subseq_id) {
            std::vector<double> frequencies((unsigned long) param.score_stl_size_, 0.0);
            subseq_count_select_query.bind(1, subseq_id);

            while (subseq_count_select_query.executeStep()) {
                int round = subseq_count_select_query.getColumn(0).getInt();
                frequencies[round] = subseq_count_select_query.getColumn(1).getDouble() / total_counts[round];
            }//while subseq_count_select_query
            subseq_count_select_query.reset();

            //frequency
            for (const auto & round : param.round_ids_) {
                //insert Frequency
                subseq_freq_insert_query.bind(1, subseq_id);
                subseq_freq_insert_query.bind(2, round);
                subseq_freq_insert_query.bind(3, frequencies[round]);
                while (subseq_freq_insert_query.executeStep()) {}
                subseq_freq_insert_query.reset();
            }//for round
        }//for subseq

        result_transaction.commit();
    } catch (std::exception &e) {
        std::cerr << "Error in calculation and inserting of subseq frequency" << std::endl;
        std::cerr << "exception:" << e.what() << std::endl;
        exit(1);
    }

    std::cout << "Finished subsequence scores calculation." << std::endl;
}


void RaptRanker::ScoreCalculateAndInserter(const Parameter_info &param, const std::string &score_name,
                                           const std::string &query) {
    try {
        std::string dbfile = param.analysis_dbfile_;
        SQLite::Database resultDB(dbfile, SQLite::OPEN_READWRITE);
        SQLite::Statement insert_query(resultDB, query);
        SQLite::Transaction result_transaction(resultDB);

        while (insert_query.executeStep()) {}

        result_transaction.commit();
    } catch (std::exception &e) {
        std::cerr << "Error in calculation and inserting of " << score_name << std::endl;
        std::cerr << "exception:" << e.what() << std::endl;
        exit(1);
    }
}

void RaptRanker::CalcClusterScore(const Parameter_info &param) {
    std::cout << "Calculating cluster scores..." << std::endl;

    //cluster count
    std::string query = "INSERT INTO cluster_count(cluster_id,round_id,value)"
                        "SELECT subseq_cluster.cluster_id,round_id,SUM(value) "
                        "FROM subseq_cluster NATURAL INNER JOIN subseq_count "
                        "GROUP BY subseq_cluster.cluster_id, round_id";
    ScoreCalculateAndInserter(param, "cluster count", query);

    //cluster deg
    query = "INSERT INTO cluster_deg(cluster_id,round_id,value) "
            "SELECT subseq_cluster.cluster_id,round_id,COUNT(subseq_id) "
            "FROM subseq_cluster NATURAL INNER JOIN subseq_count "
            "WHERE value IS NOT NULL "
            "GROUP BY subseq_cluster.cluster_id, round_id";
    ScoreCalculateAndInserter(param, "cluster degree", query);

    //cluster freq
    query = "INSERT INTO cluster_freq(cluster_id,round_id,value) "
            "SELECT subseq_cluster.cluster_id,round_id,SUM(value) "
            "FROM subseq_cluster NATURAL INNER JOIN subseq_freq "
            "GROUP BY subseq_cluster.cluster_id, round_id";
    ScoreCalculateAndInserter(param, "cluster frequency", query);

    //cluster enrich
    try {
        std::string dbfile = param.analysis_dbfile_;
        SQLite::Database resultDB(dbfile, SQLite::OPEN_READWRITE);

        SQLite::Statement select_query(resultDB, "SELECT round_id,value FROM cluster_freq WHERE cluster_id = ?1");

        SQLite::Statement insert_query(resultDB,
                                       "INSERT INTO cluster_enrich(cluster_id,round_id,value) VALUES(?1,?2,?3)");


        SQLite::Transaction result_transaction(resultDB);

        int last_cluster_id = resultDB.execAndGet("SELECT MAX(cluster_id) FROM subseq_cluster");
        for (int cluster_id = 1; cluster_id <= last_cluster_id; ++cluster_id) {
            std::vector<double> cluster_freqs((unsigned long) param.score_stl_size_, 0.0);

            select_query.bind(1, cluster_id);

            while (select_query.executeStep()) {
                int round = select_query.getColumn(0).getInt();
                cluster_freqs[round] = select_query.getColumn(1).getDouble();
            }
            select_query.reset();

            for (int index = 1, m = param.round_ids_.size(); index < m; ++index) {
                int round = param.round_ids_[index];
                int prev_round = param.round_ids_[index - 1];

                if (cluster_freqs[prev_round] != 0) {
                    insert_query.bind(1, cluster_id);
                    insert_query.bind(2, round);
                    insert_query.bind(3, cluster_freqs[round] / cluster_freqs[prev_round]);
                    while (insert_query.executeStep()) {}
                    insert_query.reset();
                }
            }
        }//for cluster

        result_transaction.commit();
    } catch (std::exception &e) {
        std::cerr << "Error in calculation and inserting of cluster enrichment" << std::endl;
        std::cerr << "exception:" << e.what() << std::endl;
        exit(1);
    }

    std::cout << "Finished cluster scores calculation." << std::endl;
}

void RaptRanker::CalcKmerScore(const RaptRanker::Parameter_info &param) {
    std::cout << "Calculating k-mer scores..." << std::endl;

    try {
        std::string dbfile = param.analysis_dbfile_;
        SQLite::Database resultDB(dbfile, SQLite::OPEN_READWRITE);


        std::string kmer_score_template = "("
                                          "kmer_sequence  TEXT     NOT NULL,"
                                          "round_id       INTEGER  NOT NULL,"
                                          "value          NUMERIC  NOT NULL,"
                                          "FOREIGN KEY(kmer_sequence) REFERENCES all_subseq(sub_sequence),"
                                          "UNIQUE(kmer_sequence, round_id)"
                                          ")";
        resultDB.exec("CREATE TABLE IF NOT EXISTS kmer_count" + kmer_score_template);
        resultDB.exec("CREATE TABLE IF NOT EXISTS kmer_deg" + kmer_score_template);
        resultDB.exec("CREATE TABLE IF NOT EXISTS kmer_freq" + kmer_score_template);
        resultDB.exec("CREATE TABLE IF NOT EXISTS kmer_enrich" + kmer_score_template);
        resultDB.exec("CREATE INDEX IF NOT EXISTS kmer_freq_value_index ON kmer_freq(value)");
        resultDB.exec("CREATE INDEX IF NOT EXISTS kmer_deg_value_index ON cluster_deg(value)");
        resultDB.exec("CREATE INDEX IF NOT EXISTS kmer_count_value_index ON kmer_count(value)");
        resultDB.exec("CREATE INDEX IF NOT EXISTS kmer_enrich_value_index ON kmer_enrich(value)");
        std::cout << "\t kmer score tables and indexes are ready." << std::endl;
    } catch (std::exception &e) {
        std::cerr << "Error in preparation of DB for k-mer scores" << std::endl;
        std::cerr << "exception:" << e.what() << std::endl;
        exit(1);
    }

    //kmer_count
    std::string query = "INSERT INTO kmer_count(kmer_sequence,round_id,value)"
                        "SELECT sub_sequence,round_id,SUM(value) FROM all_subseq "
                        "NATURAL INNER JOIN subseq_count GROUP BY sub_sequence,round_id";
    ScoreCalculateAndInserter(param, "k-mer count", query);

    //kmer_deg
    query = "INSERT INTO kmer_deg(kmer_sequence,round_id,value)"
            "SELECT sub_sequence,round_id,COUNT(subseq_id) FROM all_subseq "
            "NATURAL INNER JOIN subseq_count WHERE value IS NOT NULL "
            "GROUP BY sub_sequence,round_id";
    ScoreCalculateAndInserter(param, "k-mer degree", query);

    //kmer_freq
    query = "INSERT INTO kmer_freq(kmer_sequence,round_id,value) "
            "SELECT sub_sequence,round_id,SUM(value) FROM all_subseq "
            "NATURAL INNER JOIN subseq_freq GROUP BY sub_sequence,round_id";
    ScoreCalculateAndInserter(param, "k-mer frequency", query);

    //kmer enrich
    try {
        std::string dbfile = param.analysis_dbfile_;
        SQLite::Database resultDB(dbfile, SQLite::OPEN_READWRITE);
        SQLite::Statement kmer_select_query(resultDB, "SELECT DISTINCT(sub_sequence) FROM all_subseq");
        SQLite::Statement score_select_query(resultDB, "SELECT round_id,value FROM kmer_freq WHERE kmer_sequence = ?1");
        SQLite::Statement insert_query(resultDB,"INSERT INTO kmer_enrich(kmer_sequence,round_id,value) VALUES(?1,?2,?3)");
        SQLite::Transaction result_transaction(resultDB);

        while (kmer_select_query.executeStep()) {
            std::vector<double> kmer_freqs((unsigned long) param.score_stl_size_, 0.0);
            std::string kmer_seq = kmer_select_query.getColumn(0).getText();

            score_select_query.bind(1, kmer_seq);
            while (score_select_query.executeStep()) {
                int round = score_select_query.getColumn(0).getInt();
                kmer_freqs[round] = score_select_query.getColumn(1).getDouble();
            }
            score_select_query.reset();

            for (int index = 1, m = param.round_ids_.size(); index < m; ++index) {
                int round = param.round_ids_[index];
                int prev_round = param.round_ids_[index - 1];

                if (kmer_freqs[prev_round] != 0) {
                    insert_query.bind(1, kmer_seq);
                    insert_query.bind(2, round);
                    insert_query.bind(3, kmer_freqs[round] / kmer_freqs[prev_round]);
                    while (insert_query.executeStep()) {}
                    insert_query.reset();
                }
            }
        }//for cluster

        result_transaction.commit();
    } catch (std::exception &e) {
        std::cerr << "Error in calculation and inserting of k-mer enrichment" << std::endl;
        std::cerr << "exception:" << e.what() << std::endl;
        exit(1);
    }
    std::cout << "Finished k-mer scores calculation." << std::endl;
}


void RaptRanker::CalcSeqScoreDB(const Parameter_info &param) {
    std::cout << "Calculating sequence scores..." << std::endl;

    //AMF
    std::string query = "INSERT INTO seq_AMF(seq_id,round_id,value)"
                        "SELECT seq_id,round_id,AVG(IFNULL(value,0.0)) FROM all_subseq "
                        "NATURAL INNER JOIN subseq_cluster "
                        "NATURAL INNER JOIN cluster_freq "
                        "GROUP BY seq_id,round_id";

    ScoreCalculateAndInserter(param, "AMF", query);


    //AME
    query = "INSERT INTO seq_AME(seq_id,round_id,value)"
            "SELECT seq_id,round_id,AVG(IFNULL(value,0.0)) FROM all_subseq "
            "NATURAL INNER JOIN subseq_cluster "
            "NATURAL INNER JOIN cluster_enrich "
            "GROUP BY seq_id,round_id";

    ScoreCalculateAndInserter(param, "AME", query);


    std::cout << "Finished sequence scores calculation." << std::endl;
}

void RaptRanker::ExportScoreCSV(const Parameter_info &param) {
    std::cout << "Exporting score csv..." << std::endl;

    std::string ofile = param.analysis_output_path_ + "score.csv";
    if (file_exists(ofile)) {
        std::cerr << "\"" << ofile << "\" is already exists.\n";
        exit(1);
    }

    try {
        std::ofstream output(ofile);

        std::string dbfile = param.analysis_dbfile_;
        SQLite::Database resultDB(dbfile, SQLite::OPEN_READWRITE);
        dbfile = param.experiment_dbfile_;
        SQLite::Database seqinfoDB(dbfile, SQLite::OPEN_READWRITE);

        SQLite::Statement round_name_select_query(seqinfoDB, "SELECT round_id,round_name FROM round_info");
        SQLite::Statement seq_select_query(seqinfoDB, "SELECT seq_id,sequence,flag FROM all_seq "
                                                      "NATURAL LEFT OUTER JOIN binding_seq");
        SQLite::Statement count_select_query(seqinfoDB, "SELECT round_id, value FROM seq_count WHERE seq_id = ?1");
        SQLite::Statement freq_select_query(seqinfoDB, "SELECT round_id, value FROM seq_freq WHERE seq_id = ?1");
        SQLite::Statement enrich_select_query(seqinfoDB, "SELECT round_id, value FROM seq_enrich WHERE seq_id = ?1");
        SQLite::Statement AME_select_query(resultDB, "SELECT round_id, value FROM seq_AME WHERE seq_id = ?1");

        SQLite::Transaction seqinfo_transaction(seqinfoDB);
        SQLite::Transaction result_transaction(resultDB);


        //make header
        std::ostringstream header;

        header << "seq_id,sequence,binding_flag";

        std::vector<std::string> round_names((unsigned long) param.score_stl_size_, "");
        while (round_name_select_query.executeStep()) {
            int round = round_name_select_query.getColumn(0).getInt();
            round_names[round] = round_name_select_query.getColumn(1).getText();
        }
        for (const auto & round : param.round_ids_) {
            header << "," << round_names[round] << "_" << "Count";
        }
        for (const auto & round : param.round_ids_) {
            header << "," << round_names[round] << "_" << "Frequency";
        }
        for (int index = 1, m = param.round_ids_.size(); index < m; ++index) {
            header << "," << round_names[param.round_ids_[index]] << "_" << "Enrichment";
        }
        for (int index = 1, m = param.round_ids_.size(); index < m; ++index) {
            header << "," << round_names[param.round_ids_[index]] << "_" << "AME";
        }
        output << header.str() << "\n";

        //export results
        while (seq_select_query.executeStep()) {
            std::ostringstream osstream;

            std::vector<int> counts((unsigned long) param.score_stl_size_, 0);
            std::vector<double> freqs((unsigned long) param.score_stl_size_, 0.0);
            std::vector<double> enrichs((unsigned long) param.score_stl_size_, 0.0);
            std::vector<double> AMEs((unsigned long) param.score_stl_size_, 0.0);

            int seq_id = seq_select_query.getColumn(0).getInt();
            std::string sequence = seq_select_query.getColumn(1).getText();
            std::string flag = seq_select_query.getColumn(2).getText();

            osstream << seq_id << "," << sequence << "," << flag;

            //get count
            count_select_query.bind(1, seq_id);
            while (count_select_query.executeStep()) {
                int round = count_select_query.getColumn(0).getInt();
                counts[round] = count_select_query.getColumn(1).getInt();
            }
            count_select_query.reset();

            //get freq
            freq_select_query.bind(1, seq_id);
            while (freq_select_query.executeStep()) {
                int round = freq_select_query.getColumn(0).getInt();
                freqs[round] = freq_select_query.getColumn(1).getDouble();
            }
            freq_select_query.reset();

            //get enrich
            enrich_select_query.bind(1, seq_id);
            while (enrich_select_query.executeStep()) {
                int round = enrich_select_query.getColumn(0).getInt();
                enrichs[round] = enrich_select_query.getColumn(1).getDouble();
            }
            enrich_select_query.reset();

            //get AME
            AME_select_query.bind(1, seq_id);
            while (AME_select_query.executeStep()) {
                int round = AME_select_query.getColumn(0).getInt();
                AMEs[round] = AME_select_query.getColumn(1).getDouble();
            }
            AME_select_query.reset();


            for (const auto & round : param.round_ids_) {
                osstream << "," << counts[round];
            }
            for (const auto & round : param.round_ids_) {
                osstream << "," << freqs[round];
            }
            for (int index = 1, m = param.round_ids_.size(); index < m; ++index) {
                osstream << "," << enrichs[param.round_ids_[index]];
            }
            for (int index = 1, m = param.round_ids_.size(); index < m; ++index) {
                osstream << "," << AMEs[param.round_ids_[index]];
            }

            output << osstream.str() << "\n";
        }//

        seqinfo_transaction.commit();
        result_transaction.commit();
    } catch (std::exception &e) {
        std::cerr << "Error in exporting score.csv" << std::endl;
        std::cerr << "exception:" << e.what() << std::endl;
        exit(1);
    }

    std::cout << "Finished score csv exporting." << std::endl;
}

void RaptRanker::AddBindingFlagDB(const Parameter_info &param) {
    std::cout << "Adding Binding_flag ..." << std::endl;

    //input samples
    std::ifstream input(param.bindingfile_path_);
    std::string line;
    std::vector<std::pair<std::string, int> > sample_data; //sample datas
    if (input.fail()) {
        std::cerr << "binding information data input failed\n";
    }
    while (getline(input, line)) {
        std::vector<std::string> split_line = split(line, ',');
        std::string sequence = split_line[0];
        int judge = stoi(split_line[1]);
        sample_data.emplace_back(make_pair(sequence, judge));
    }
    input.close();


    try {
        std::string dbfile = param.experiment_dbfile_;
        SQLite::Database seqinfoDB(dbfile, SQLite::OPEN_READWRITE);
        SQLite::Statement select_query(seqinfoDB, "SELECT seq_id FROM all_seq WHERE sequence = ?1");
        SQLite::Statement insert_query(seqinfoDB, "INSERT OR IGNORE INTO binding_seq(seq_id,flag) VALUES(?1,?2)");
        SQLite::Transaction transaction(seqinfoDB);

        for (const auto &each_sample : sample_data) {
            //add binding to ht_Seqinfo
            select_query.bind(1, each_sample.first);
            while (select_query.executeStep()) {
                int seq_id = select_query.getColumn(0).getInt();

                insert_query.bind(1, seq_id);
                insert_query.bind(2, each_sample.second);
                while (insert_query.executeStep()) {}
                insert_query.reset();
            }
            select_query.reset();
        }// for each_sample

        transaction.commit();
    } catch (std::exception &e) {
        std::cerr << "Error in adding Binding_flag" << std::endl;
        std::cerr << "exception:" << e.what() << std::endl;
        exit(1);
    }

    std::cout << "Finished Binding_flag addition." << std::endl;
}