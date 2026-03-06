# GitHub 업로드 방법 (SCT27)

## 전체 순서

1. GitHub 계정 준비
2. 로컬에서 git 초기화
3. GitHub에 repository 생성
4. push

---

## Step 1: GitHub 계정 및 설정

아직 없으면 https://github.com 에서 계정 만들기.

로컬 git 설정 (처음 한 번만):
```bash
git config --global user.name "Your Name"
git config --global user.email "your@email.com"
```

---

## Step 2: 로컬 git 초기화

```bash
cd sct27_github/          # 이 패키지 폴더로 이동

git init                  # git 저장소 초기화
git add .                 # 모든 파일 스테이징
git status                # 스테이징된 파일 확인 (선택)
git commit -m "Initial release: SCT27 NLTE code (PRE 2024)"
```

---

## Step 3: GitHub에 Repository 생성

### 웹 브라우저에서:
1. https://github.com → 로그인
2. 오른쪽 상단 **+** → **New repository**
3. 입력:
   - **Repository name:** `sct27`  (또는 `sct27-nlte-eedf` 등)
   - **Description:** `NLTE plasma code (SCFLY+ZELDA) for XFEL experiments. Cho et al. PRE 2024.`
   - **Public** 선택 (논문 재현 목적이므로 공개 권장)
   - **⚠️ "Add a README file" 체크 해제** (우리가 이미 만들었음)
   - **⚠️ ".gitignore" 체크 해제** (이미 있음)
   - **⚠️ "Choose a license" 체크 해제** (이미 있음)
4. **Create repository** 클릭

---

## Step 4: 로컬에서 GitHub로 연결 & Push

GitHub 페이지에서 복사한 URL을 사용:

```bash
# HTTPS 방식 (간단, 비밀번호 또는 token 필요):
git remote add origin https://github.com/YOUR_USERNAME/sct27.git

# 또는 SSH 방식 (key 설정 완료된 경우 더 편리):
git remote add origin git@github.com:YOUR_USERNAME/sct27.git

# 업로드:
git branch -M main
git push -u origin main
```

### HTTPS 인증 (token):
GitHub는 2021년부터 비밀번호 대신 Personal Access Token을 사용합니다.
- GitHub → Settings → Developer settings → Personal access tokens → Tokens (classic)
- **Generate new token** → `repo` 권한 체크 → 생성
- push할 때 비밀번호 대신 이 token 입력

---

## Step 5: 이후 코드 수정 시 업데이트

```bash
# 파일 수정 후:
git add 수정된파일.f90       # 특정 파일만
# 또는
git add .                    # 전체

git commit -m "Fix: description of what you changed"
git push
```

---

## 권장 Repository 설정 (GitHub 웹에서)

Push 완료 후 GitHub 웹에서:

1. **About** (오른쪽 ⚙️):
   - Description: `NLTE plasma code coupling SCFLY population kinetics with ZELDA Boltzmann EEDF solver`
   - Website: DOI 링크 `https://doi.org/10.1103/PhysRevE.109.045207`
   - Topics: `plasma-physics`, `xfel`, `boltzmann-equation`, `nlte`, `fortran`, `eedf`

2. **Releases** (오른쪽 패널):
   - **Create a new release** → tag: `v1.0.0`, title: `PRE 2024 release`
   - 이렇게 하면 논문에 `Code available at: github.com/username/sct27 (v1.0.0)` 로 인용 가능

3. **README** 확인: GitHub에서 자동으로 렌더링되어 보임

---

## Zenodo 연동 (DOI 부여, 선택 사항)

학술 코드는 Zenodo를 통해 DOI를 받으면 논문에서 인용하기 좋습니다:

1. https://zenodo.org → GitHub 계정으로 로그인
2. **GitHub** 탭 → 저장소 토글 ON
3. GitHub에서 Release 생성하면 Zenodo가 자동으로 DOI 부여

---

## 전체 파일 구조 확인

```
sct27/
├── README.md               ← GitHub 메인 페이지에 표시됨
├── LICENSE
├── Makefile                ← 수정된 빌드 파일
├── .gitignore
├── src/                    ← 모든 소스 코드
│   ├── scmn_v3.f           (수정: ee_scale 추가)
│   ├── zelda.f90           (수정: 3개 버그 수정 + ee_scale 적용)
│   ├── zcoeff_elastic.f90  (수정: tgev 버그 수정)
│   ├── zchk_stuff.f90      (수정: ee_scale 변수 추가)
│   └── ... (기타 소스 파일들)
├── input/
│   └── neon_2keV/          ← PRE 2024 Neon 케이스 입력 파일들
├── scripts/
│   ├── run_fig1b.sh        ← Fig.1(b) 재현 스크립트
│   └── run_fig1c.sh        ← Fig.1(c) 재현 스크립트
└── docs/
    └── physics_notes.md    ← 물리/수치 방법 설명
```
