

## What is FastPerm?
FastPerm은 single-cell RNA-seq data에 대해 비교군/대조군 간 통계 검정을 수행하는 프로그램입니다. scRNA-seq data를 pseudobulk하고, non-parametric test인 permutation test(순열 검정법)을 사용합니다. FastPerm을 활용해 DE analysis를 수행하거나 GSEA(Gene Set Enrichment Analysis)등을 수행할 수 있습니다. 특히, pseudobulk시에 샘플의 수가 적어져 parametric test를 사용하기 어려운 상황에서 사용할 수 있습니다.


### Main features
- single-cell method들(wilcoxon, MAST)에 비해 낮은 alpha error
- bulk method들(DESeq2, edgdR)에 비해 낮은 beta error
- [Sparse Matrix](https://en.wikipedia.org/wiki/Sparse_matrix) -> [LIL (List of Lists)](https://en.wikipedia.org/wiki/Sparse_matrix#List_of_lists_(LIL)) conversion을 통한 연산 효율 향상
- Multithreading을 활용한 연산 효율 향상
- Runtime동안 고정적인 memory usage


### Development Version Info
- OS : Ubuntu 22.04.1 LTS
- Compiler : g++ (Ubuntu 11.3.0-1ubuntu1~22.04) 11.3.0



&nbsp;
## How to build
git clone

```
> git clone git@github.com:merri4/FastPerm.git
```

해당 디렉토리 이동 후, make main 입력

```
> make main
```

다음과 같은 로그가 나오면 정상입니다.

```
g++    -c -o main.o main.cpp
g++    -c -o myfuncs.o myfuncs.cpp
g++ -o main main.o myfuncs.o -lpthread
```

&nbsp;
## How to use
필수적으로 넣어줘야 하는 파라미터는 5가지입니다.
3가지는 input 파일, 1가지는 output 경로, 나머지 하나는 permutation 횟수입니다.

- `NameCC table file_name` : 샘플의 이름과 case / control 여부가 담긴 txt 파일 
- `BarcodeName table file_name` : cell의 barcode와 각 cell이 어느 샘플에 속하는지를 나타낸 txt 파일
- `GeneBarcode file_name` : 행이 gene이고, 열이 cell인 count matrix txt 파일
- `Result save path` : 결과 파일이 저장될 경로 및 파일명
- `Permutation Number` : permutation할 횟수
- `cpm count threshold` : pseudobulk 후 cpm (count per million) 과정에서 저발현하는 유전자들을 검정에서 제외하기 위해 사용하는 threshold입니다. (0, 100,000] 사이의 정수로만 입력되어야 합니다. 별도로 입력하지 않을 경우 1로 들어갑니다.

다시 말해 `NameCC table`과 `BarcodeName Table`은 metadata이고, `GeneBarcode`는 RNA count matrix입니다.

test/ 안에 3가지 input 파일의 예시가 있으니 포맷을 확인하시기 바랍니다.

예를 들어 이런 커맨드를 사용하면, 

```
./main perm_sampledata.txt perm_coldata.txt perm_counts.txt result.csv 1000
```

- `perm_sampledata.txt`를 `NameCC table`로,
- `perm_coldadta.txt`을 `BarcodeName table`로,
- `perm_coutns.txt`를 `GeneBarcode table`로,
- `result.csv`를 output 경로로,
- `permutation 횟수`를 `1000회`로

지정하고 프로그램이 실행됩니다. 

&nbsp;
## Troubleshooting
오류가 났을 경우 체크해봐야 할 주요 사항은 다음과 같습니다.
- `NameCC table` 의 구분자는 쉼표, case와 control은 반드시 **소문자로만** 작성되어야 합니다.
- `BarcodeName table`에서의 cell에 대응하는 sample 정보는 `NameCC table`의 샘플 이름들과 반드시 일치해야 합니다. `NameCC table`에서 나오지 않은 샘플명이 `BarcodeName table`에서 나오면 프로그램이 오작동합니다. 
- `GeneBarcode table`의 첫 줄은 cell의 이름(barcode)를 나타냅니다. `BarcodeName table`에서의 cell 이름과 개수가 정확히 일치해야 합니다 (순서가 일치해야 할 필요는 없음).


&nbsp;
## Reading logs
```
// 인풋에 오류가 없을 경우 이 메시지가 출력됩니다. 받은 파라미터들을 보여줍니다.
Input Taken.
Permutation Number :    1000
Cpm count threshold :   1
Start Processing...


// 받은 샘플의 정보를 보여줍니다. case는 0, control은 1로 표시됩니다.
// (인덱스) : (샘플명) : (case / control) 여부 형태로 출력됩니다.
./data/perm_sampledata.txt Opened. Reading file..  
0 : s1 : 0
1 : s2 : 0
2 : s3 : 0
3 : s4 : 0
4 : s5 : 0
5 : s6 : 1
6 : s7 : 1
7 : s8 : 1
8 : s9 : 1
9 : s10 : 1

// 받은 cell 이름과 sample 소속 정보를 보여줍니다. 
// (Cell 이름) : (소속 sample의 인덱스) 로 표시됩니다.
// 예를 들어 C1_1은 0번 인덱스이므로, s1 샘플에 속하는 cell이라는 의미입니다.
./data/perm_coldata.txt Opened. Reading file..
C1_1 : 0        (s1)
C1_10 : 0       (s1)
C1_100 : 0      (s1)
C1_101 : 1      (s2)
C1_102 : 1      (s2)
1000 rows

// 가장 시간이 오래 걸리는 지점입니다. 일반적인 크기일 경우 2-3분, 길게는 10분까지 소요될 수 있습니다.
./data/perm_counts.txt Opened. Reading file..
FileRead Elapsed Time : 5296 ms

// 실제 count matrix에 대한 pseudobulk 수행 시간을 나타냅니다.
Pseudobulk Elapsed Time : 8 ms

// pseudobulk 후 logcpm normalization 과정이 정상적으로 수행되었음을 보여줍니다.
logcpm Done!

// 원본 데이터에 대한 lfc table을 정상적으로 생성했음을 보여줍니다.
LFC Table Done!

// permutation 시작됩니다. 10%마다 로그가 출력되며, 남은 예상 소요 시간이 초 단위로 출력됩니다.
Permutation Starting...
10 % Done       Permutation : 100 / 1000        ETA : 17 seconds
20 % Done       Permutation : 200 / 1000        ETA : 15 seconds
30 % Done       Permutation : 300 / 1000        ETA : 13 seconds
40 % Done       Permutation : 400 / 1000        ETA : 11 seconds
50 % Done       Permutation : 500 / 1000        ETA : 9 seconds
60 % Done       Permutation : 600 / 1000        ETA : 7 seconds
70 % Done       Permutation : 700 / 1000        ETA : 5 seconds
80 % Done       Permutation : 800 / 1000        ETA : 3 seconds
90 % Done       Permutation : 900 / 1000        ETA : 1 seconds
Permutation End. Elapsed Time : 18511 ms // 총 걸린 시간 표시

// 결과가 정상적으로 저장되었음을 보여줍니다.
P-value Table Saved !

// 프로그램 실행이 완료되었습니다.
Programme Done! Terminating...
```



&nbsp;
## Reading output result
```
Gene,Case_Mean,Control_Mean,pct_case,pct_ctrl,lfc,p_val_1t,p_val_2t
gene1,8.066215,7.748651,5,5,0.317564,0.100000,0.200000
gene2,5.580030,5.030147,5,5,0.549883,0.200000,0.700000
gene3,inf,inf,0,0,NA,NA,NA
...
```

|  Gene  | case mean | control mean | pct_case | pct_ctrl |   lfc    | p_val_1t | p_val_2t |
| :----: | :-------: | :----------: | :------: | :------: | :------: | :------: | :------: |
| gene 1 | 8.066215  |   7.748651   |    5     |    5     | 0.317564 | 0.100000 | 0.200000 |

- `Gene` : gene 이름
- `case mean` : case 측의 mean
- `control mean` : control 측의 mean
- `pct_case` : threshold를 넘은 case 측 샘플의 개수
- `pct_control` : threshold를 넘은 control 측 샘플의 개수
- `lfc` : case와 control의 log2-fold change. (NA라고 표기되었다면 pct_case나 pct_control 중 어느 하나라도 0인 경우입니다. 따라서 pvalue도 구할 수 없기에 NA로 표기됩니다.)
- `p_val_1t` : 1 tailed p-value. lfc의 rank를 나열할 때 양수/음수측을 구분하여 정렬해 계산한 p-value
- `p_val_2t` : 2 tailed p-value. lfc의 rank를 나열할 때 절댓값을 기준으로 정렬해 계산한 p-value

일반적으로 one-tailed보다 two-tailed가 좀 더 보수적인 결과를 보여줍니다. alpha error와 FDR control을 고려하여, 기준이 더 엄격한 **2-tailed p-value를 사용하길 권장**합니다.





&nbsp;
## Development History
- binary search로 열별 누산.
- ifstream 파일 여는 코드 더 간단하게 수정 (while(file.good()) 대신 while(getline())으로)
- 동기식 처리에서 비동기 처리로 전환 (parse_gene_barcode_table, pseudobulk). 비동기 처리 시 pseudobulk 후 mismatch가 발생하여, 뮤텍스 사용
- load_data_from_text 데이터 불러오는 하나의 래퍼함수로 묶음, 파일 실행 시 경로를 arg로 받도록 수정
- split 함수의 개선
- mmap을 사용하여 30배정도 빠른 txt reading
- 세마포어를 사용해 스레드 개수를 변화시키면서 성능을 테스트해보았는데, 스레드가 3개일 때 가장 빠르게 작동하지만, 어떤 경우에도 순차처리보다 빠르지는 않았음. 순차처리를 하는 걸 기본으로 하되, sparse table를 파싱하는 개선된 방식을 사용해 테스트해보기.
- sparse table을 LIL(List of Lists)로 파싱하여 뮤텍스가 필요 없는 새로운 병렬처리 도입(사람별로 해당하는 바코드의 인덱스를 받아 순회).
- LFC table 만들 때 log2의 NaN또는 INF값의 처리 >> permutation시 rank 매길 수 있도록 아예 LfcVal 구조체를 만듬. 또, NA값은 rank에서 제외하므로 그 개수도 세도록 함. 결측값들은 제외
- 전역변수를 정리하고 하드코딩을 많이 제거. 코드 가독성 높임.
- lfc table 출력 가능하게 만듬
- 사용하는 size 관련 변수들 모두 전역으로 뺌
- 완료된 gene별 lfc, p-value table을 txt로 export함
- permutation 수행, 파라미터로 넣을 수 있음
- output 출력 변경에 따른 파라미터, 클래스 수정
- 변수 정리 및 전역변수화, two-tail 결과 출력, ETA 10% 단위로 출력하도록 수정
