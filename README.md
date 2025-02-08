### Version Info
OS : Ubuntu 22.04.1 LTS
Compiler : g++ (Ubuntu 11.3.0-1ubuntu1~22.04) 11.3.0

### Structure
구조 설명
- 클래스. (2차원 벡터)
- 어떻게 데이터가 파싱되는지
- pseudobulk의 과정

### History
- binary search로 열별 누산.
- 필리아 테이블 (예시 데이터)로 pseudobulk 검증했을 때 문제 없음. >> 데이터 문제로 밝혀짐
- ifstream 파일 여는 코드 더 간단하게 수정 (while(file.good()) 대신 while(getline())으로)
- 동기식 처리에서 비동기 처리로 전환 (parse_gene_barcode_table, pseudobulk). 비동기 처리 시 pseudobulk 후 mismatch가 발생하여, 뮤텍스 사용. 뮤텍스를 사용하였더니 성능이 저하됨...??
- load_data_from_text 데이터 불러오는 하나의 래퍼함수로 묶음, 파일 실행 시 경로를 arg로 받도록.
- split 함수의 개선
- mmap을 사용하여 30배정도 빠른 txt reading
- 세마포어를 사용해 스레드 개수를 변화시키면서 성능을 테스트해보았는데, 스레드가 3개일 때 가장 빠르게 작동하지만, 어떤 경우에도 순차처리보다 빠르지는 않았음. 순차처리를 하는 걸 기본으로 하되, sparse table를 파싱하는 개선된 방식을 사용해 테스트해보기.
- 새로운 병렬처리 도입 (사람별로 해당하는 바코드의 인덱스를 받아 순차처리). 대강 1400 ms 나옴
- LFC table 만들 때 log2의 NaN또는 INF값의 처리 >> permutation시 rank 매길 수 있도록 아예 LfcVal 구조체를 만듬. 또, NA값은 rank에서 제외하므로 그 개수도 세도록 함. 결측값들은 제외
- 전역변수를 정리하고 하드코딩을 많이 제거. 코드 가독성 높임.
- lfc table 출력 가능하게 만듬
- 사용하는 size 관련 변수들 모두 전역으로 뺌
- 완료된 gene별 lfc, p-value table을 txt로 export함
- permutation 수행, 파라미터로 넣을 수 있음
- output 출력 변경에 따른 파라미터, 클래스 수정
- 변수 정리 및 전역변수화, two-tail 결과 출력, ETA 10% 단위로 출력하도록 수정

### Features
- 행렬곱을 대체하는 Binary Search를 통한 빠른 pseudobulk + 열(바코드)별 멀티스레딩
- nmap을 사용한 빠른 txt file reading + 멀티스레딩
- log계산 시 0인 것 생략
- signal handler
- sparse matrix 특징 이용한 적은 메모리 소요 및 pseudobulk 시간 단축


### TODOLIST
- [ ] API 추가
- [ ] 쓰지 않는 print 관련 함수들 제거
- [ ] usage 출력 포맷 다른 것들 참고