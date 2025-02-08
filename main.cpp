/*
    
    *** WARNINGS ***
    WARNING : 읽을 txt 파일은 마지막 줄이 빈 줄이어서는 안됨 (segfault). 추후 막줄이 비어있으면 중지하도록 수정할 것.
    

    *** TODOLIST ***

    모든 반복성 객체 생성자들 초기화 리스트로 바꾸기
    [Done] 구조체 수정함에 따라 아웃풋 바꾸기
    [Done] make_lfc 함수 대대적으로 수정, 카운팅해야하는게 많고, NA값이 나올 일이 없어졌음
    [Done] threshold 퍼센티지는 받지 않도록 수정하기

    인풋 체커 만들기 & 예외처리
    하드코딩된 부분 전역변수나 #define으로 빼내기

    
    [Done] gene x name table에서 로그화 시 어느 요소를 각각 행/열로 잡는것이 유리한지?
    [Done] groupsum 계산할 때, if val != 0 일  때만 sum += 누산하기
    [Done] Gene 이름 저장하는 방식 ? 일단은 첫 열만 저장
    [Done] GENE 이름 저장하기
    [Done] Input sample 만들기


    *** LINUX DEPENDENCIES ***
    [Done] 하나의 단위작업을 멀티스레딩 또는 멀티프로세싱 (permutation별로 멀티스레딩)
    [Done] signal 등록 (SIGINT, SIGKILL)
    [Done] 파라미터를 argc와 argv로 받도록

*/

// ============================================================================================================================
// ============================================================================================================================

// Header

// ============================================================================================================================
// ============================================================================================================================

// 입출력 / 파일 관련
#include <iostream>
#include <fstream>
#include <sstream>

// 타입 관련
#include <cstring>
#include <string>
#include <vector>
#include <chrono>

// 연산 관련
// #include <algorithm>    // binary search
#include <random>       // shuffle
#include <cmath>        // log 
#include <numeric>      // accumulation 

// 멀티스레딩 관련
#include <signal.h>
#include <thread>
#include <mutex>

// nmap 관련
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>

#include "myfuncs.h"

using namespace std;


// ============================================================================================================================
// ============================================================================================================================

// Const Variables Configuration

// ============================================================================================================================
// ============================================================================================================================

#define DELIMITER       ','             // 나눌 DELIMITER

#define PRINT_LOG       true            // 로그 출력 여부

#define NA              log(1.0/0.0)    // NA값 (inf로 저장)

#define N_ESSEN_PARAM   6               // bin을 포함한 필수 파라미터 개수

const vector<string> DEFAULT_PARAM_LIST = {  // 파라미터 리스트
    "bin",                      // file name
    "Name_CC_table",            // Name : Case/Control Table Directory
    "Barcode_Name_table",       // Barcode : Name Table Directory
    "Gene_Barcode_table",       // Gene : Barcode Table Directory
    "Result_save_path",         // result save path
    "100",                      // Permutation Nombre
    "1",                        // cpm count threshold number
    };

// ============================================================================================================================
// ============================================================================================================================

// Global Variables

// ============================================================================================================================
// ============================================================================================================================

vector<string> PARAM_LIST;

int N_SAMPLE;                               // sample 개수 (사람수)
int N_GENE;                                 // Gene 개수 
int N_BARCODE;                              // Barcode 개수 (cell 수)

int CPM_COUNT_THRESHOLD;                    // cpm 후 RNA read의 임계값

int N_PERMUTATION;                          // permutation 총 횟수
int I_PERMUTATION;                          // permutation 현재 횟수
chrono::steady_clock::time_point time_permutation_start;    // permutation 시작 시간

vector<NameType> namecc_table;              // [ IDX : NAME : CASE/CONTROL]
vector<BarcodeType> barname_table;          // [ BARCODE : IDX(NAME) ]
GeneBarcodeTable genebarcode_table;         // [ Gene x Barcode ]   
vector<vector<double>> pseudobulk_table;    // [ Gene x Name ]
vector<LfcVal> lfc_table;                   // [ IFC TABLE per gene ]
mutex lfc_table_mtx;
vector<double> lfc_table_tmp;               // [ IFC TABLE per gene for imsi ]

vector<thread> vec_bulk_thread;             // pseudobulk 스레드 리스트
vector<thread> vec_perm_thread;             // permutation 스레드 리스트

// ============================================================================================================================
// ============================================================================================================================

// Function declaration

// ============================================================================================================================
// ============================================================================================================================

// 시그널 핸들러
void SignalHandler(int signum);

// ETA 출력기
void process_eta();

// 입력 처리 관련 함수
int input_checker(int argc, char* argv[]);

// 테이블 처리 관련 함수
void load_data_from_txt(string& path1, string& path2, string& path3);
void parse_name_cc_table(string& path);
void parse_barcode_name_table(string& path);
void parse_gene_barcode_table(string& path);
void processLine(const string& line, int n_line);

// PseudoBulk 관련 함수
void pseudobulk();
void pseudobulk_process_col(const vector<int>& idx_list, int n_sample, const int& n_gene);
int pseudobulk_validation(string& path, bool printlog);
void print_pseudobulk();

// PseudoBulk Table 결과값 가공 처리 함수
void logcpm();
void make_lfc();
void print_lfc();
void save_lfc(string path);

// Permutation 관련 함수
void pseudobulk_shuffle();
void make_lfc_tmp();
void process_permutation();

// 결과 출력 함수
void print_result_pval(string path);

// ============================================================================================================================
// ============================================================================================================================

// Function Definition

// ============================================================================================================================
// ============================================================================================================================

// 시그널 핸들러
void SignalHandler(int signum) {
    cout << "\nSIGINT Detected. Threads Terminating... ";
    for (thread& t : vec_bulk_thread) t.join();
    for (thread& t : vec_perm_thread) t.join();
    cout << "Done!" << endl;
    exit(signum);
}

// ETA 출력기
void process_eta() {

    // 10개 포인트 심음
    int eta_timepoint[10]; int eta_idx = 0;
    for (int i=0; i<10; i++) eta_timepoint[i] = (i+1)*0.1*N_PERMUTATION; 

    // busy waiting..
    while (I_PERMUTATION < N_PERMUTATION-1) {

        // 10% 도달 시
        if (I_PERMUTATION == eta_timepoint[eta_idx]) {
            eta_idx++;    
            int eta = ( chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now() - time_permutation_start).count() * (10-eta_idx) / eta_idx ) / 1000;
            printf("%d %% Done\tPermutation : %d / %d\tETA : %d seconds\n", eta_idx*10, I_PERMUTATION, N_PERMUTATION, eta);
        }
    }

}

// 입력 처리기
int input_checker(int argc, char* argv[]){

    cout << endl;

    // 필수 개수보다 적으면 오류
    if (argc < N_ESSEN_PARAM) {
        cerr << "Usage : " << argv[0] 
        << " <NameCC table file_name>" 
        << " <BarcodeName table file_name>" 
        << " <GeneBarcode table file_name>" 
        << " <Result Save Path>"
        << " <Permutation Number>"
        << " <optional : cpm count threshold = 1>"
        << endl;
        return -1;
    }

    // default parameter parsing
    int i = 0;
    while (i < DEFAULT_PARAM_LIST.size()) {
        if (i < argc) PARAM_LIST.emplace_back(argv[i]);
        else PARAM_LIST.emplace_back(DEFAULT_PARAM_LIST[i]);
        i++;
    }


    // Parameter Setting
    N_PERMUTATION = stoi(PARAM_LIST[5]);
    CPM_COUNT_THRESHOLD = stoi(PARAM_LIST[6]);


    // Parameter Validation : todo 하드코딩된 부분 빼낼 것
    // permutation 횟수는 1회 이상.
    if (!N_PERMUTATION) {
        cerr << "Parameter Error : Permutation number should be over 0." << endl;
        return -1;
    }

    // CPM 역치 범위는 (0, 100,000] 사이.
    if ( (CPM_COUNT_THRESHOLD <= 0) || (CPM_COUNT_THRESHOLD > 100'000) ) {
        cerr << "Parameter Error : Cpm count threshold number should be in range (0, 100,000]." << endl;
        return -1;
    }

    // Log
    printf("Input Taken.\n");
    printf("Permutation Number :\t%d\n", N_PERMUTATION);
    printf("Cpm count threshold :\t%d\n", CPM_COUNT_THRESHOLD);
    printf("Start Processing...\n\n");
    return 1;
    
}

// Wrapper
void load_data_from_txt(string& path1, string& path2, string& path3) {
    
    //cout << "[ IDX : NAME : CASE/CONTROL] TABLE Path : "; cin >> path; 
    // path = "./data/d_NameCC.txt";
    parse_name_cc_table(path1);
    cout << endl;

    //cout << "[ BARCODE : IDX(NAME) ] TABLE Path : "; cin >> path;
    // path = "./data/d_BarcodeName.txt";
    parse_barcode_name_table(path2);
    cout << endl;

    //cout << "[ Gene x Barcode ] Table Path : "; cin >> path;
    // path = "./data/d_GeneBarcode.txt";
    parse_gene_barcode_table(path3);
    cout << endl;

    return;
}

// Name : CaseControl Table
void parse_name_cc_table(string &path) {

    // ==========================================================================================
    // [ IDX : NAME : CASE/CONTROL] TABLE =======================================================
    // ==========================================================================================

    // 파일 열기 ===================================================================
    ifstream file(path);

    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << path << endl;
        exit(0);
    }

    if (PRINT_LOG) cout << path << " Opened. Reading file..\n"; // log


    // Name : CC row 파싱 =====================================================
    string line;
    while (getline(file, line)) {

        vector<string> tmpline = split(line, DELIMITER); // 쉼표 단위로 분리 (csv 호환)
        CaseControl tmpCC = (tmpline[1] == "case") ? CASE : CONTROL; // case control 지정
        NameType tmpname(tmpline[0], tmpCC); // 타입 선언하고
        namecc_table.push_back(tmpname); // 벡터에 넣음
    }
    file.close(); // 닫기


    // 결과 출력 =====================================================
    if (PRINT_LOG) {
        for (int i = 0; i < namecc_table.size(); i++) {
            cout << i << " : " << namecc_table[i].m_name << " : " << namecc_table[i].m_casecontrol << endl;
        }
    }

}

// Barcode : Name Table
void parse_barcode_name_table(string& path) {

    // ==========================================================================================
    // [ BARCODE : IDX(NAME) ] TABLE ============================================================
    // ==========================================================================================


    // 파일 열기 ===================================================================
    ifstream file(path);

    if (!file.is_open()) {
        cerr << "Failed to open file: " << path << endl;
        exit(0);
    }

    if (PRINT_LOG) cout << path << " Opened. Reading file..\n"; // log


    // BARCODE : IDX(NAME) row 파싱 =====================================================    
    string line;
    while (getline(file, line)) {
        vector<string> tmpline = split(line, DELIMITER); // 쉼표 단위로 분리 (csv 호환)
        BarcodeType tmpbarcode(tmpline[0], get_index(namecc_table, tmpline[1])); // 그 정보를 받아 객체 생성
        barname_table.push_back(tmpbarcode); // 벡터에 넣음
    }
    file.close();


    // BARCODE : IDX(NAME) row 파싱 =====================================================    
    sort(barname_table.begin(), barname_table.end()); // sort 불필요시 주석처리할 것.


    // 결과 출력 ========================================================================
    if (PRINT_LOG) {
        for (int i = 0; i < 5; i++) {
            cout << barname_table[i].m_barcode << " : " << barname_table[i].m_idx << "\t(" << namecc_table[barname_table[i].m_idx].m_name << ")" << endl;
        }
        cout << barname_table.size() << " rows" << endl;
    }

}

// Gene : Barcode Table
void parse_gene_barcode_table(string& path) {

    // ==========================================================================================
    // [ Gene x Barcode ] Table
    // ==========================================================================================

    chrono::steady_clock::time_point start_time = chrono::steady_clock::now();


    // Open the file for reading ======================================
    int fd = open(path.c_str(), O_RDONLY);
    if (fd == -1) {
        cerr << "Failed to open file: " << path << endl;
        exit(0);
    }

    // Get the size of the file
    off_t size = lseek(fd, 0, SEEK_END);
    if (size == (off_t)-1) {
        perror("lseek");
        exit(0);
    }
    
    // Map the file into memory
    void* data = mmap(NULL, size, PROT_READ, MAP_PRIVATE, fd, 0);
    if (data == MAP_FAILED) {
        perror("mmap");
        exit(0);
    }

    if (PRINT_LOG) cout << path << " Opened. Reading file..\n"; // log


    // 라인 개수 세기  ================================
    char* ptr = (char*)data;
    char* end = (char*)data + size;
    int line_counter = 0;
    while (ptr < end) {
        // Find the end of the line
        char* newline = (char*)memchr(ptr, '\n', end - ptr);
        if (newline == NULL) newline = end;
        line_counter++;
        ptr = newline + 1;
    }
    genebarcode_table.mv_gene_names.resize(line_counter-1); // gene의 개수는 총 라인-1 (바코드 줄 제외)

    // 본격적으로 줄별 처리  ================================
    ptr = (char*)data;
    end = (char*)data + size;

    // 첫 줄 (바코드 이름) 처리
    char* newline = (char*)memchr(ptr, '\n', end - ptr);
    string line(ptr, newline - ptr);
    genebarcode_table.mv_barcodes = split(line, DELIMITER);
    genebarcode_table.mvv_table.resize(genebarcode_table.mv_barcodes.size(), vector<GeneVal>()); // mvv table은 바코드 개수만큼.

    // 남은 줄 처리
    ptr = newline + 1;
    line_counter = 0;

    while (ptr < end) {
        // Find the end of the line
        char* newline = (char*)memchr(ptr, '\n', end - ptr);
        if (newline == NULL) newline = end;
        
        // !!! Process the lines !!!
        string line(ptr, newline - ptr);

        // 이 라인을 어떻게 하냐면...
        processLine(line, line_counter);
        
        // Move to the start of the next line
        ptr = newline + 1;
        line_counter++;
    }

    if (PRINT_LOG){
        chrono::steady_clock::time_point end_time = chrono::steady_clock::now();
        cout << "FileRead Elapsed Time : " << chrono::duration_cast<chrono::milliseconds>(end_time - start_time).count() << " ms\n";
    } 


    // 결과 출력 ===================================================================
    // if (PRINT_LOG) {
    //     cout << "\t\t";
    //     for (int j = 0; j < 8; j++) {
    //         cout << genebarcode_table.mv_barcodes[j] << "\t";
    //     }
    //     cout << endl;

    //     for (int i = 0; i < 10; i++) {

    //         // 이름
    //         cout << genebarcode_table.mv_gene_names[i] << "\t";
    //         // values
    //         for (int j = 0; j < 8; j++) {
    //             cout << genebarcode_table.mvv_table[i][j] << "\t";
    //         }
    //         cout << endl;
    //     }
    //     cout << genebarcode_table.mvv_table.size() << " rows" << endl;
    //     cout << genebarcode_table.mv_barcodes.size() << " cols" << endl;
    // }

}

// 스레드 서브루틴
void processLine(const string& line, int n_line) {
    
    // Processing..
    vector<string> tmpline = split(line, DELIMITER); // 쉼표 단위로 분리
    genebarcode_table.mv_gene_names[n_line] = tmpline[0]; // 이름 넣기
    
    // 돌아가면서
    for (int i=1; i<tmpline.size(); i++) {
        if (tmpline[i] != "0") { // 0 아닌 것들만
            genebarcode_table.mvv_table[i-1].emplace_back(n_line, stoi(tmpline[i]));
        }
    }

    // genebarcode_table.take_row(tmpline, n_line); // 그 줄을 넣음

}

// Pseudobulk 
void pseudobulk() {

    chrono::steady_clock::time_point start_time = chrono::steady_clock::now();

    // psd <- aggregate.Matrix(t(mtx), groupings = groups, fun = "sum")

    // patient 수 * gene의 개수 ( 6 * 19263 ) TABLE 생성
    pseudobulk_table.resize(N_SAMPLE, vector<double>(N_GENE, 0));

    // 사람 수에 따른 벡터를 만든다.
    vector<vector<int>> baridx_list(N_SAMPLE, vector<int>());

    // barcode name 벡터를 돌리면서 해당하는 인덱스들을 넣는다.
    for (int idx_barcode = 0; idx_barcode < N_BARCODE; idx_barcode++) {
        BarcodeType compBarcode(genebarcode_table.mv_barcodes[idx_barcode], 0); // 지금의 바코드명을 비교용 바코드 타입에 넣는다.
        int who = barname_table[binarySearch_idx(barname_table, compBarcode)].m_idx; // 표 상 위치를 binary search로 찾아, 해당하는 col 인덱스를 가져온다.
        baridx_list[who].emplace_back(idx_barcode);
    }


    // 이제 사람 수만큼 독립적으로 멀티스레딩. 뮤텍스 칠 필요가 없음.
    for (int n_sample=0; n_sample < N_SAMPLE; n_sample++) {
        namecc_table[n_sample].m_barcount = baridx_list[n_sample].size(); // 속한 바코드 개수를 저장해둠
        vec_bulk_thread.emplace_back(pseudobulk_process_col, baridx_list[n_sample], n_sample, N_GENE); // 스레드에 할당
    }

    for (thread& t : vec_bulk_thread) t.join();
    vec_bulk_thread.clear();


    if (PRINT_LOG){
        chrono::steady_clock::time_point end_time = chrono::steady_clock::now();
        cout << "Pseudobulk Elapsed Time : " << chrono::duration_cast<chrono::milliseconds>(end_time - start_time).count() << " ms\n";
    } 

}

// 스레드 서브루틴 (pseudobulk)
void pseudobulk_process_col(const vector<int>& idx_list, int n_sample, const int& n_gene) {

    // vector<int> idx_list :   현재 연산해야 할 barcode들의 queue 
    // int n_sample         :   해당하는 샘플
    // int n_gene           :   총 gene의 개수

    // 처리해야 할 바코드 인덱스 큐를 돌리면서
    for (const int& col : idx_list) {

        // cell별로
        for (const GeneVal& geneval : genebarcode_table.mvv_table[col]) {

            // 결과 테이블에다가 누산
            pseudobulk_table[n_sample][geneval.m_geneidx] += geneval.m_val;

        }
            

    }

}

// R 결과값과 consistency 체크
int pseudobulk_validation(string& path, bool printlog) {
    
    ifstream file(path);
    
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << path << endl;
        return -1;
    }

    string line;
    int line_count = 0;
    int mismatch_count = 0;
    getline(file, line); // 첫 줄은 버린다.

    while (getline(file, line)) {
        
        vector<string> tmpline = split(line, DELIMITER); // 쉼표 단위로 분리 (csv 호환)
        
        for (int i=1; i<tmpline.size(); i++) {
            if ( pseudobulk_table[i-1][line_count] != stoi(tmpline[i]) ) {

                if (printlog) {
                    printf(
                        "Mismatch ! (%s, %s)\t계산값 : %f\t검증값 : %s\n",
                        genebarcode_table.mv_gene_names[line_count].c_str(),
                        namecc_table[i-1].m_name.c_str(),
                        pseudobulk_table[i-1][line_count],
                        tmpline[i].c_str());
                }
                mismatch_count++;
            }
        }
        line_count++;
        // if (line_count > 20) break;
    }
    file.close(); // 닫기

    // 결과에 따른 리턴
    if (mismatch_count) {
        cout << "Mismatch Count : " << mismatch_count << endl;
        return -1;
    } 
    
    return 1; // 정상일 경우
}

// Pseudobulk printing
void print_pseudobulk() {
    
    for (int name = 0; name < 6; name++) {
        cout << namecc_table[name].m_name << "\t";
    }
    cout << endl;

    for (int i = 0; i < 6; i++) {
        cout << genebarcode_table.mv_gene_names[i] << "\t";
        for (int j = 0; j < 6; j++) {
            cout << pseudobulk_table[j][i] << "\t";
        }
        cout << endl;
    }
    
}

void logcpm() {
    
    // cpm > log transformation
    // psd = log2(edgeR::cpm(psd) + 1)

    for (vector<double>& vec : pseudobulk_table) {
        double sum = accumulate(vec.begin(), vec.end(), 0);
        for (double& val : vec) {
            if (val != 0) {
                val = log2((val * 1'000'000 / sum) + 1);
            }
        }
    }

}

void make_lfc() {

    // case_m = psd[, pheno == "case"] % > % rowMeans
    // control_m = psd[, pheno == "control"] % > % rowMeans
    // lfc = log2(case_m / control_m)


    // 평균을 구하기 위한 case/control 별 전체 개수
    int n_case_tot = 0;
    int n_ctrl_tot = 0;

    // 개수 집계
    for (const NameType& namecc : namecc_table) {
        switch (namecc.m_casecontrol) {
            case CASE :
                n_case_tot++; break;
            case CONTROL :
                n_ctrl_tot++; break;
        }
    }
    
    // lfc table 크기 설정
    lfc_table.resize(N_GENE, 0);

    // gene별로
    for (int gene = 0; gene < N_GENE; gene++) {

        int case_over = 0;
        int ctrl_over = 0;

        double sum_ca    = 0;
        double sum_ctrl  = 0;

        // 사람마다 돌아가면서
        for (int person = 0; person < N_SAMPLE; person++) {

            // 자기 인덱스의 Name이 case, control이냐에 따라서 각자 값에다가 누산
            switch (namecc_table[person].m_casecontrol) {

                case CASE:
                    if (pseudobulk_table[person][gene] >= CPM_COUNT_THRESHOLD) case_over++;
                    sum_ca += pseudobulk_table[person][gene]; break;

                case CONTROL:
                    if (pseudobulk_table[person][gene] >= CPM_COUNT_THRESHOLD) ctrl_over++;
                    sum_ctrl += pseudobulk_table[person][gene]; break;

            }
        }

        // 넘은 값들 설정
        lfc_table[gene].m_case_n_over_thres = case_over;
        lfc_table[gene].m_ctrl_n_over_thres = ctrl_over;

        // case/control의 평균값을 받아 lfc값을 구해주는 메소드
        lfc_table[gene].process_lfc_val(sum_ca / n_case_tot, sum_ctrl / n_ctrl_tot);

    }
}

void print_lfc() {

    for (int i = 0; i < lfc_table.size(); i++) {
        cout << genebarcode_table.mv_gene_names[i] << "\t" << lfc_table[i].m_val << "\n";
    }

}

void save_lfc(string path) {

    ofstream file(path, ios::binary | ios::out);

    // 반복하면서
    for (int i = 0; i < lfc_table.size(); i++) {
        string line = genebarcode_table.mv_gene_names[i] + DELIMITER + to_string(lfc_table[i].m_val) + "\n";
        file.write(line.c_str(), line.size());
    }   
    
    file.close();

    cout << "LFC Table Saved !\n" << endl;

}



void pseudobulk_shuffle() {

    // chrono::steady_clock::time_point start_time = chrono::steady_clock::now();

    // 슈도벌크 테이블을 새로 만든다.
    // 이 때, 스레드로 할당하는 인덱스 위치를  start와 end로 주면 된다. end-start가 그 사람의 개수가 되도록.
    // pseudobulk_table.resize(namecc_table.size(), vector<double>(genebarcode_table.mv_gene_names.size(), 0));
    for (auto &vec : pseudobulk_table) fill(vec.begin(), vec.end(), 0); // 초기화
    
    // size init time : 0 ms

    // 인덱스 리스트를 셔플한다.
    random_device rd; mt19937 mersenne(rd());
    vector<int> shuffled_idx_list(N_BARCODE, 0);
    vector<vector<int>> shuffle_slice(N_SAMPLE, vector<int>());
    for (int i=0; i<N_BARCODE; i++) shuffled_idx_list[i] = i;
    shuffle(shuffled_idx_list.begin(), shuffled_idx_list.end(), mersenne);

    // shuffle time : 0 ms

    // 그걸 개수에 맞게 분배한다.
    int start = 0;
    int end = 0;
    for (int n_sample=0; n_sample < N_SAMPLE; n_sample++) {
        
        // 인덱싱할 상-하한선 결정
        start = end;
        end += namecc_table[n_sample].m_barcount;
        
        // 복사해가기
        for (int i=start; i<end;i++) shuffle_slice[n_sample].emplace_back(shuffled_idx_list[i]);

    }

    // 사람별 인덱스 리스트를 크기순으로 정렬한다.
    for (int n_sample=0; n_sample < N_SAMPLE; n_sample++) {        
        sort(shuffle_slice[n_sample].begin(), shuffle_slice[n_sample].end());
    }

    // index sort time : 3 ms

    // 이제 셔플된 인덱스 리스트를 할당하면서 멀티스레딩.
    for (int n_sample=0; n_sample < N_SAMPLE; n_sample++) {
        vec_bulk_thread.emplace_back(pseudobulk_process_col, shuffle_slice[n_sample], n_sample, N_GENE); // 스레드에 할당
    }

    // alloc time : 0 ms

    for (thread& t : vec_bulk_thread) t.join();
    vec_bulk_thread.clear();

    // thread calc time : 1771 ms

    // if (PRINT_LOG){
    //     chrono::steady_clock::time_point end_time = chrono::steady_clock::now();
    //     cout << "Pseudobulk Elapsed Time : " << chrono::duration_cast<chrono::milliseconds>(end_time - start_time).count() << " ms\n";
    // }
    
}

void make_lfc_tmp() {

    // 평균을 구하기 위한 case/control 별 전체 개수
    int n_case_tot = 0;
    int n_ctrl_tot = 0;

    // 개수 집계
    for (const NameType& namecc : namecc_table) {
        switch (namecc.m_casecontrol) {
            case CASE :
                n_case_tot++; break;
            case CONTROL :
                n_ctrl_tot++; break;
        }
    }

    // lfc table 초기화
    lfc_table_tmp.resize(N_GENE, 0);
    fill(lfc_table_tmp.begin(), lfc_table_tmp.end(), 0);

    // gene별로
    for (int gene = 0; gene < N_GENE; gene++) {

        double sum_ca    = 0;
        double sum_ctrl  = 0;

        // 사람마다 돌아가면서
        for (int person = 0; person < N_SAMPLE; person++) {

            // 자기 인덱스의 Name이 case, control이냐에 따라서 각자 값에다가 누산
            switch (namecc_table[person].m_casecontrol) {

                case CASE:
                    sum_ca += pseudobulk_table[person][gene]; break;

                case CONTROL:
                    sum_ctrl += pseudobulk_table[person][gene]; break;

            }
        }

        // 평균을 구한 후
        sum_ca /= n_case_tot;
        sum_ctrl /= n_ctrl_tot;

        // 해당 값을 저장
        lfc_table_tmp[gene] = sum_ca - sum_ctrl;

    }

}

void process_permutation() {
    
    // LFC랑 LFC_tmp를 비교한다.
    for (int gene = 0; gene < N_GENE; gene++) { // gene마다 돌아가면서..

        int compval = lfc_table_tmp[gene];

        // 자기보다 큰 것의 개수를 센다. (one-tailed)
        if (lfc_table[gene].m_val < compval) lfc_table[gene].m_n_over_1t++;
        
        // (절댓값 기준) 자기보다 더 큰 것의 개수를 센다. (two-tailed)
        if (abs(lfc_table[gene].m_val) < abs(compval)) lfc_table[gene].m_n_over_2t++;

    }

}


// 결과를 출력한다.....
void print_result_pval(string path) {


    ofstream file(path, ios::binary | ios::out);

    // 출력 포맷은
    // Gene name, case mean, control mean, case에서 threshold 넘은 개수, ctrl에서 넘은 개수, lfc value, pvalue \n
    string colline = "Gene_Name,Case_Mean,Control_Mean,overthreshold_number_in_case,overthreshold_number_in_control,lfc_value,p-value_one_tailed,p-value_two_tailed\n";
    file.write(colline.c_str(), colline.size());


    for (int i_gene = 0; i_gene < N_GENE; i_gene++) { // gene별로 반복하면서
        
        double pval_1t = 0;
        double pval_2t = ( ((double)lfc_table[i_gene].m_n_over_2t) + 1) / N_PERMUTATION;

        // 상위 50%면, 
        // 자신의 rank / Permutation 횟수 = ( 나보다 컸던 숫자 + 1 ) / Permutation 횟수
        if (lfc_table[i_gene].m_n_over_1t < N_PERMUTATION/2) {
            pval_1t = ((double) lfc_table[i_gene].get_rank()) / N_PERMUTATION;
        }
        
        // 하위 50%면, rank를 뒤집어서,
        // 자신의 rank / Permutation 횟수 = (Permutation 횟수 - 나보다 컸던 숫자) / Permutation 횟수
        else pval_1t = 1 - ( ( (double) lfc_table[i_gene].m_n_over_1t ) / N_PERMUTATION);


        // 라인 만들어 출력 WARNING : DELIMITER를 합칠 때, "NA"와 같이 직접적인 string과 합치면 힙을 읽는 에러가 발생함. output 출력 이상해짐
        // string myline = to_string(63) + DELIMITER;  // 작동
        // string myline = "NA" + DELIMITER;           // 오작동
        string line;
        line = genebarcode_table.mv_gene_names[i_gene] + DELIMITER
                + to_string(lfc_table[i_gene].m_mean_case) + DELIMITER
                + to_string(lfc_table[i_gene].m_mean_ctrl) + DELIMITER
                + to_string(lfc_table[i_gene].m_case_n_over_thres) + DELIMITER
                + to_string(lfc_table[i_gene].m_ctrl_n_over_thres) + DELIMITER
                + to_string(lfc_table[i_gene].m_val) + DELIMITER
                + to_string(pval_1t) + DELIMITER
                + to_string(pval_2t)
                + "\n";

        file.write(line.c_str(), line.size());
        
    }   
    
    file.close();
    
    printf("P-value Table Saved !\n");

}


// ============================================================================================================================
// ============================================================================================================================

// Main Function

// ============================================================================================================================
// ============================================================================================================================

// command : ./main ./d_NameCC.txt ./d_BarcodeName.txt ./d_GeneBarcode.txt ./result_p_val.txt 100 1

int main(int argc, char* argv[]) {

    if (input_checker(argc, argv) < 0) return 1;

    // signal 등록
    signal(SIGINT, SignalHandler);

    // 필요한 정보 불러오기
    load_data_from_txt(PARAM_LIST[1], PARAM_LIST[2], PARAM_LIST[3]);

    
    // 기본 1회 ===========================================================
    
    // 변수 세팅
    N_SAMPLE = namecc_table.size();
    N_BARCODE = genebarcode_table.mv_barcodes.size();
    N_GENE = genebarcode_table.mv_gene_names.size();
    
    pseudobulk();
    // string path = "./data/psd_check2.txt";
    // if (pseudobulk_validation(path, false) < 0){
    //     cerr << "Error. Terminating..." << endl; exit(0);
    // }
    // print_pseudobulk();

    // 이후 처리
    logcpm();
    cout << "\nlogcpm Done!\n";
    // print_pseudobulk();

    make_lfc();
    cout << "\nLFC Table Done!\n";
    // print_lfc();
    // path = "lfc.txt";
    // save_lfc(path);


    // 나머지 횟수를 permutation =========================================
    cout << "\nPermutation Starting...\n";

    // 시간 측정
    thread t_eta(process_eta);
    time_permutation_start = chrono::steady_clock::now();

    for (I_PERMUTATION=0; I_PERMUTATION<N_PERMUTATION-1; I_PERMUTATION++) {
        
        pseudobulk_shuffle(); // 153 ms
        logcpm(); // 15 ms
        make_lfc_tmp(); // 4 ms
        process_permutation(); // 1 ms
    
    }

    t_eta.join();
    cout << "Permutation End. Elapsed Time : " << chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now() - time_permutation_start).count() << " ms\n";

    // 결과 프린트
    print_result_pval(PARAM_LIST[4]);

    printf("Programme Done! Terminating...\n");
}