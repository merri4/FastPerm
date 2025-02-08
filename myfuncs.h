#pragma once

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <time.h>
#include <chrono>

using namespace std;

// ==========================================================================================
// case = 0. control = 1
// ==========================================================================================
enum CaseControl { CASE, CONTROL };


// ==========================================================================================
// NameType
// ==========================================================================================
class NameType {

public :

    // attrs
    string m_name;
    short m_casecontrol;
    int m_barcount;

    // methods
    NameType(string name, CaseControl casecontrol) : m_name(name), m_casecontrol(casecontrol), m_barcount(0) {}
    ~NameType();

};


// ==========================================================================================
// 바코드 (barcode, idx)
// ==========================================================================================

class BarcodeType {

public :

    // attr
    string m_barcode;
    short m_idx;

    // methods
    BarcodeType(string barcode, short idx) : m_barcode(barcode), m_idx(idx) {}
    ~BarcodeType();
    bool operator>(BarcodeType& barcodetype);
    bool operator<(BarcodeType& barcodetype);
    bool operator>=(BarcodeType& barcodetype);
    bool operator<=(BarcodeType& barcodetype);
    bool operator==(BarcodeType& barcodetype);
    bool operator!=(BarcodeType& barcodetype);

};


// ==========================================================================================
// GeneVal
// ==========================================================================================
class GeneVal {
public :

    int m_geneidx;
    int m_val;

    // Methods
    GeneVal() : m_geneidx(0), m_val(0) {};
    GeneVal(int idx, int val) : m_geneidx(idx), m_val(val) {};
    ~GeneVal();

};


// ==========================================================================================
// GeneBarcodeTable
// ==========================================================================================
class GeneBarcodeTable {

public :

    // Attributes
    vector<string> mv_gene_names; // gene 이름을 저장 (행)
    vector<string> mv_barcodes; // 바코드 인덱스를 저장 (열)
    vector<vector<GeneVal>> mvv_table;

    // Methods
    GeneBarcodeTable();
    ~GeneBarcodeTable();
    // void take_row(vector<string>& numarr, int n_line);

};


// ==========================================================================================
// LFC Value
// ==========================================================================================

class LfcVal {
public :

    // Attr
    double m_val;               // 실제 자기 값.
    
    double m_mean_case;         // case mean
    double m_mean_ctrl;         // control mean
    
    int m_case_n_over_thres;    // case에서 threshold 넘은 개수
    int m_ctrl_n_over_thres;    // control에서 threshold 넘은 개수

    int m_n_over_1t;            // 자기보다 더 큰 것의 개수. 곧 랭크가 됨. (전체 N회 - 결측값 개수) 중에서 자기보다 위가 3개 있었으면 4등.
    int m_n_over_2t;            // (절댓값 기준; two-tailed) 자기보다 더 큰 것의 개수.

    // Methods
    LfcVal() : m_val(0), m_mean_case(0), m_mean_ctrl(0), m_case_n_over_thres(0), m_ctrl_n_over_thres(0), m_n_over_1t(0), m_n_over_2t(0) {};
    LfcVal(double val) : m_val(val), m_mean_case(0), m_mean_ctrl(0), m_case_n_over_thres(0), m_ctrl_n_over_thres(0), m_n_over_1t(0), m_n_over_2t(0) {};
    ~LfcVal();

    int get_rank(); // 랭크를 구해줌
    void process_lfc_val(double m_case, double m_ctrl); // 두개의 평균값을 받아 lfc value를 계산 후 저장


    LfcVal& operator=(const double val); // 숫자 대입
    bool operator>(LfcVal& other);
    bool operator<(LfcVal& other);
    bool operator>=(LfcVal& other);
    bool operator<=(LfcVal& other);
    bool operator==(LfcVal& other);
    bool operator!=(LfcVal& other);

};

// ==========================================================================================
// 함수들
// ==========================================================================================

// 스트링을 입력받아서, namecc_table을 참고해 idx값을 리턴받는다.
int get_index(const vector<NameType>& table, const string& name);

// 특정 Delimiter로 분리
vector<string> split(const string& str, const char& Delimiter);

// 이진탐색
int binarySearch_idx(vector<BarcodeType>& vec, BarcodeType& data);