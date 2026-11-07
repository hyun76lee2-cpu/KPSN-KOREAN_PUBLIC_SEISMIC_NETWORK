# IMM 관측 지도

참여자들이 1초 간격 IMM 관측값을 올리고, 누구나 지도에서 볼 수 있는 사이트예요.
올리는 값: IMM_(0.2,0.5), IMM_(0.5,1.5), IMM_(1.5,7.5), 대표 IMM(0.1~10초 최댓값)과 대표 주기.
개인 관측 기록 공유용이며 지진 예보·경보가 아니에요.

## 파일

| 파일 | 역할 |
|---|---|
| `index.html` | 지도, 관측소 목록, 폰 센서·CSV로 참여하기 |
| `imm-core.js` | IMM 계산 (IMM 진도계와 같은 방식) |
| `firebase-config.js` | 내 Firebase 프로젝트 값 (직접 채움) |
| `database.rules.json` | 데이터베이스 보안 규칙 |
| `pc/imm_share.example.json` | `imm_live.py --share`용 설정 예시 |

## 1. Firebase 준비 (무료 Spark 요금제)

1. https://console.firebase.google.com 에서 프로젝트를 만들어요 (Google 애널리틱스는 꺼도 돼요).
2. **빌드 > Authentication > 시작하기 > 로그인 방법**에서 **익명**을 사용 설정해요.
3. **빌드 > Realtime Database > 데이터베이스 만들기**. 위치는 한국에서 가까운 **asia-southeast1(싱가포르)**, 보안 규칙은 "잠금 모드"로 시작해요.
4. Realtime Database의 **규칙** 탭에 `database.rules.json` 내용을 통째로 붙여넣고 게시해요.
5. **프로젝트 설정(톱니바퀴) > 내 앱 > 웹 앱 추가(</>)**. 앱 이름만 넣고 등록하면 나오는 `firebaseConfig` 값을 `firebase-config.js`에 옮겨 적어요. `databaseURL`은 Realtime Database 화면 맨 위에 나오는 주소예요.

> Firestore가 아니라 **Realtime Database**를 써야 해요. Firestore는 쓰기 횟수로 과금돼서 1초마다 쓰면 무료 한도를 금방 넘어요.

## 2. GitHub Pages로 올리기

1. GitHub에 새 공개 저장소를 만들고 이 폴더의 파일을 모두 올려요 (`pc` 폴더는 안 올려도 돼요).
2. 저장소의 **Settings > Pages**에서 Source를 "Deploy from a branch", 브랜치는 `main`, 폴더는 `/ (root)`로 저장해요.
3. 1~2분 뒤 `https://<아이디>.github.io/<저장소이름>/` 에서 열려요.
4. Firebase **Authentication > 설정 > 승인된 도메인**에 `<아이디>.github.io`를 추가해요.

## 3. 참여하기

- **폰**: 사이트를 폰 브라우저로 열고 "내 관측소로 참여하기"에서 이름과 위치를 정한 뒤 참여 시작. 폰은 바닥에 테이프로 고정하고, 처음 한 번 조용할 때 잡음 측정(30초)을 해요. 아이폰은 시작할 때 동작 센서 권한을 물어봐요.
- **CSV**: 입력을 "CSV 파형 재생"으로 바꾸고 `시간, X, Y, Z` 파일을 골라요.
- **ADXL355 (컴퓨터)**: `pc/imm_share.example.json`을 `imm_share.json`으로 복사해 값을 채우고
  `python imm_live.py --serial COM5 --share imm_share.json`
  처음 실행할 때 익명 계정이 만들어지고 `imm_share_token.json`에 저장돼서, 다음에도 같은 관측소로 이어져요.

## 데이터와 개인정보

- 공개되는 것: 관측소 이름, 위치(소수 둘째 자리, 약 1 km 단위), 센서 종류, 1초 관측값.
- 30초 넘게 갱신이 없으면 지도에서 회색, 하루가 지나면 숨겨요. 데이터 자체는 Firebase 콘솔에서 지울 수 있어요.
- 각자 자기 관측소 경로에만 쓸 수 있어요 (익명 로그인 uid 기준 보안 규칙).

## 용량 참고 (무료 한도: 저장 1 GB, 내려받기 월 10 GB)

관측소 하나가 1초에 약 100바이트를 올리고, 지도를 보는 사람은 모든 관측소 값을 계속 받아요.
관측소 10곳을 한 사람이 하루 종일 보면 하루 약 100 MB라, 시험 규모에서는 충분하지만
참여자와 보는 사람이 많아지면 갱신 간격을 늘리거나 유료 요금제를 검토해야 해요.
