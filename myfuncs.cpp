#include "myfuncs.h"

// ==========================================================================================
// NameType Methods ==================================================
// ==========================================================================================

NameType::~NameType() {}



// ==========================================================================================
// BarcodeType Methods ==================================================
// ==========================================================================================

BarcodeType::~BarcodeType() {}

bool BarcodeType::operator>(BarcodeType& barcodetype) {
    return (this->m_barcode > barcodetype.m_barcode);
}

bool BarcodeType::operator<(BarcodeType& barcodetype) {
    return (this->m_barcode < barcodetype.m_barcode);
}

bool BarcodeType::operator<=(BarcodeType& barcodetype) {
    return (this->m_barcode <= barcodetype.m_barcode);
}

bool BarcodeType::operator>=(BarcodeType& barcodetype) {
    return (this->m_barcode >= barcodetype.m_barcode);
}

bool BarcodeType::operator==(BarcodeType& barcodetype) {
    return (this->m_barcode == barcodetype.m_barcode);
}

bool BarcodeType::operator!=(BarcodeType& barcodetype) {
    return (this->m_barcode != barcodetype.m_barcode);
}


// ==========================================================================================
// GeneBarcodeTable Methods ==================================================
// ==========================================================================================

GeneBarcodeTable::GeneBarcodeTable() {
}
GeneBarcodeTable::~GeneBarcodeTable() {
}

// void GeneBarcodeTable::take_row(vector<string>& numarr, int n_line) {

//     const size_t size = numarr.size(); // todo : 여기 매번 계산하지 않고 상수 크기로 받아와도됨
//     this->mv_gene_names[n_line] = numarr[0]; // todo : gene name이 변할 수도 있음.

//     vector<int> tmparr(size-1);
//     for (int i = 0; i < size - 1; i++) {
//         tmparr[i] = stoi(numarr[i + 1]); // str을 int로 변환하여 넣고
//     }
//     this->mvv_table[n_line] = tmparr;
// }


GeneVal::~GeneVal(){}


// ==========================================================================================
// LFC Value
// ==========================================================================================

LfcVal::~LfcVal(){
}

int LfcVal::get_rank() {
    return this->m_n_over_1t+1; // 랭크를 구함
} 


void LfcVal::process_lfc_val(double m_case, double m_ctrl){
    this->m_mean_case = m_case;
    this->m_mean_ctrl = m_ctrl;
    this->m_val = m_case - m_ctrl;
}


LfcVal& LfcVal::operator=(const double val) { // 숫자 대입
    this->m_val = val;
    return *this;  // Return a reference to myself.
}

bool LfcVal::operator>(LfcVal& other) {
    return (this->m_val > other.m_val);
}

bool LfcVal::operator<(LfcVal& other) {
    return (this->m_val < other.m_val);
}

bool LfcVal::operator>=(LfcVal& other) {
    return (this->m_val >= other.m_val);
}

bool LfcVal::operator<=(LfcVal& other) {
    return (this->m_val <= other.m_val);
}

bool LfcVal::operator==(LfcVal& other) {
    return (this->m_val == other.m_val);
}

bool LfcVal::operator!=(LfcVal& other) {
    return (this->m_val != other.m_val);
}


// ==========================================================================================
// 기타 함수들 ===============================================================================
// ==========================================================================================
int get_index(const vector<NameType>& table, const string& name) {
    for (int i = 0; i < table.size(); i++) {
        if (name == table[i].m_name) {
            return i;
        }
    }
    return -1;
}

vector<string> split(const string& str, const char& Delimiter) {

    istringstream iss(str); // istringstream에 str을 담는다.
    string buffer;          // 구분자를 기준으로 절삭된 문자열이 담겨지는 버퍼
    vector<string> result;  // to return

    // istringstream은 istream을 상속받으므로 getline을 사용할 수 있다.
    while (getline(iss, buffer, Delimiter)) {
        if (buffer != "") result.emplace_back(buffer); // 절삭된 문자열을 vector에 저장
    }

    return result;
}


int binarySearch_idx(vector<BarcodeType>& vec, BarcodeType& data) {

    // 초기 인덱스 세팅
    int lo = 0;
    int hi = vec.size() - 1;
    int mid;

    // 볼 게 남아있는 동안 반복
    while (hi >= lo) {
        mid = (hi + lo) / 2; // mid 재계산
        if (vec[mid] == data) return mid; // 찾으면 즉시 종료
        else if (vec[mid] < data) lo = mid + 1; // 하한선 올림
        else hi = mid - 1; // 상한선 내림
    }

    return -1; // 못찾았으면 -1 던짐.
}
