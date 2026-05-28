# enumerate_li_configs_Final_optimized_motif_v5_mtp — 사용 매뉴얼

> **대상 파일:** `enumerate_li_configs_Final_optimized_motif_v5_mtp.py`  
> **버전:** v5 (MTP 에너지 백엔드 + 모티프 시딩 GA + v5 최적화)  
> **기반 논문:** Kim et al., *ACS Nano* — "Revealing Li Staging Reactions in Graphite via a GA Coupled with MLIP"

---

## 목차

1. [개요](#1-개요)
2. [요구 사항](#2-요구-사항)
3. [입력 파일](#3-입력-파일)
4. [출력 디렉토리 구조](#4-출력-디렉토리-구조)
5. [전체 옵션 레퍼런스](#5-전체-옵션-레퍼런스)
6. [사용 예시](#6-사용-예시)
7. [워크플로우 가이드](#7-워크플로우-가이드)
8. [성능 튜닝 가이드](#8-성능-튜닝-가이드)
9. [출력 파일 설명](#9-출력-파일-설명)
10. [자주 묻는 질문](#10-자주-묻는-질문)

---

## 1. 개요

흑연 삽입 화합물(GIC, Graphite Intercalation Compound)에서 Li-vacancy 배열을 열거하고, 각 SOC(State of Charge) 레벨에서 가장 안정한 구조를 탐색하는 스크립트입니다.

### 동작 모드

| 모드 | 설명 | 사용 시기 |
|------|------|----------|
| **전수 열거 (full enumeration)** | 대칭 유일 대표 구조를 모두 생성 | 소형 슈퍼셀, Li 수가 적을 때 |
| **GA (Genetic Algorithm)** | 유전 알고리즘으로 광대한 공간 탐색 | 대형 슈퍼셀, Li 수가 많을 때 |

### 에너지 백엔드

| 백엔드 | 설명 | 기본값 |
|--------|------|--------|
| **`mtp-relax`** | `mlp relax`로 MLIP 구조 완화 후 에너지 | ✅ 기본 |
| **`ewald`** | Ewald 전기정적 에너지 (빠르나 근사) | 대안 |

---

## 2. 요구 사항

```bash
# Python 3.8+
pip install pymatgen numpy matplotlib tqdm

# MLIP 실행파일 (mtp-relax 백엔드 사용 시)
# mlp 바이너리가 PATH에 있어야 함
mlp --version

# 필수 파일
# - 구조 파일: LiC6.cif 또는 POSCAR
# - (mtp-relax 시) mlip.ini: MLIP 포텐셜 설정 파일
```

---

## 3. 입력 파일

### 구조 파일 (`--structure`)
- LiC6 또는 LiC12 등 GIC 단위셀 CIF 또는 POSCAR
- Li와 C 원소가 모두 포함되어야 함
- 슈퍼셀 확장 전 단위셀 형태로 제공

### mlip.ini (`--mlip-ini`)
- MLIP 포텐셜 파라미터 파일
- `mlp relax` 실행에 필요
- **mtp-relax 백엔드 사용 시 필수**

---

## 4. 출력 디렉토리 구조

```
li_configs/                              ← --output (기본: ./li_configs)
│
├── enumeration_summary.json             ← 전체 실행 요약 (SOC별 최안정 구조, 통계)
│
├── SOC_0.5000_nLi24/                    ← SOC=0.5, n_Li=24 인 경우
│   ├── top_stable_cif/
│   │   ├── rank001_AAAA_nLi24_SOC0.5000_MTP_E-123.456789eV.cif   ← 최안정 구조
│   │   ├── rank002_AAAA_nLi24_SOC0.5000_MTP_E-123.400000eV.cif
│   │   └── ...  (최대 20개)
│   ├── ga_convergence_nLi24.png         ← GA 수렴 그래프 (--ga 시)
│   │
│   ├── mtp_cfg/                         ← MTP 입출력 CFG 파일
│   │   ├── IN_ga_gen000.cfg             ← 각 세대 입력 구조 (--mtp-write-aggregate 시)
│   │   ├── IN_ga_gen001.cfg
│   │   ├── ...
│   │   └── relaxed_SOC_0.5000_nLi24_all.cfg  ← 완화된 구조 전체 (SOC 완료 후 통합)
│   │
│   └── mtp_logs/
│       └── relax_ga_gen000.log          ← mlp relax 실행 로그 (SOC별 단일 파일)
│
└── SOC_0.4167_nLi20/
    └── ...
```

> **주의:** MTP 임시 작업 파일(`mtp_work/`)은 기본적으로 실행 후 자동 삭제됩니다.  
> `--mtp-keep-workdirs` 옵션으로 보존 가능.

---

## 5. 전체 옵션 레퍼런스

> 🔵 = **기본값 있음 (생략 가능)**  
> 🔴 = **필수 또는 기본값 없음**  
> ✅ = **기본 활성화된 플래그**

---

### 5.1 기본 설정

| 옵션 | 타입 | 기본값 | 설명 |
|------|------|--------|------|
| `--structure` / `-s` 🔴 | str | — | 단위셀 CIF 또는 POSCAR 파일 경로 |
| `--output` / `-o` 🔵 | str | `./li_configs` | 출력 디렉토리 루트 |
| `--supercell NX NY NZ` 🔵 | int×3 | `2 2 2` | 슈퍼셀 확장 비율 (a, b, c 방향) |

---

### 5.2 열거 모드

| 옵션 | 타입 | 기본값 | 설명 |
|------|------|--------|------|
| `--no-symmetry` 🔵 | flag | off | 대칭 비고려 전수 열거 (C(N,k) 전체) |
| `--symprec` 🔵 | float | `0.01` Å | SpacegroupAnalyzer 대칭 정밀도 |
| `--match-tol` 🔵 | float | `0.01` | Li 사이트 매칭 분수좌표 허용오차 |
| `--soc-levels N [N ...]` 🔵 | int+ | 모든 0…N_Li | 처리할 Li 수 목록 (예: `--soc-levels 12 16 20 24`) |

---

### 5.3 GA (유전 알고리즘) 설정

| 옵션 | 타입 | 기본값 | 설명 |
|------|------|--------|------|
| `--ga` 🔵 | flag | off | GA 활성화 (미지정 시 전수 열거) |
| `--ga-pop` 🔵 | int | `120` | 세대당 population 크기 |
| `--ga-gens` 🔵 | int | `60` | 총 세대 수 |
| `--ga-threshold-ev` 🔵 | float | `0.5` eV | 생존자 선택 에너지 임계값 |
| `--ga-seed` 🔵 | int | `42` | 난수 시드 (재현성) |
| `--stacking-mutation-prob` 🔵 | float | `0.15` | 스태킹 돌연변이 비율 (새 개체 중) |
| `--ga-min-li-li-same-layer` 🔵 | float | `4.00` Å | 같은 Li 레이어 내 최소 Li-Li 거리 |

---

### 5.4 스태킹 시퀀스 설정

| 옵션 | 타입 | 기본값 | 설명 |
|------|------|--------|------|
| `--stacking-sequences SEQ [SEQ ...]` 🔵 | str+ | 단층(AAAA 등) | 탐색할 C 레이어 스태킹 시퀀스. `ALL`로 모든 A/B 조합 사용 |
| `--stacking-modes AA\|AB` 🔵 | str+ | — | 레거시 단축 옵션 (AA=전층A, AB=ABAB...) |
| `--ab-shift DA DB` 🔵 | float×2 | 자동감지 | B 레이어 분수좌표 면내 변위 (미지정 시 자동) |
| `--carbon-layer-tol` 🔵 | float | `0.05` | 탄소 레이어 그룹핑 분수 z 허용오차 |

#### SOC 기반 스태킹 정책

| 옵션 | 타입 | 기본값 | 설명 |
|------|------|--------|------|
| `--soc-stacking-policy` 🔵 | choice | `conservative` | SOC에 따른 스태킹 공간 축소 정책 |
| `--aa-only-soc-threshold` 🔵 | float | 정책에 따름 | AA-only 적용 SOC 임계값 직접 지정 |

**`--soc-stacking-policy` 선택지:**

| 값 | AA-only 시작 SOC | 근거 |
|----|-----------------|------|
| `none` | 없음 (모든 SOC에 전체 시퀀스) | 정책 없음 |
| `conservative` ✅ | SOC ≥ 0.50 | 보수적 마진 |
| `literature` | SOC ≥ 0.40 | Kim et al. 논문 Fig. 4c 기준 |

---

### 5.5 Li 사이트 확장 (GA 전용)

| 옵션 | 타입 | 기본값 | 설명 |
|------|------|--------|------|
| `--ga-expand-li-sites` ✅ | flag | **on** | GA Li 후보 사이트를 흑연 헥사곤 중심으로 확장 |
| `--no-ga-expand-li-sites` | flag | — | 입력 구조의 Li 사이트만 사용 |
| `--hollow-site-tol` 🔵 | float | `0.35` Å | 빈 사이트 중복 제거 허용 거리 |
| `--cc-bond-cutoff` 🔵 | float | `1.85` Å | C-C 결합 판별 기준거리 (헥사곤 감지용) |

---

### 5.6 모티프 시딩 (Motif-Seeded GA)

고SOC에서 물리적으로 타당한 초기 population을 제공하여 수렴 속도를 향상시킵니다.

| 옵션 | 타입 | 기본값 | 설명 |
|------|------|--------|------|
| `--ga-motif-seeding` 🔵 | choice | `high-soc` | 모티프 시딩 모드 |
| `--motif-seed-soc-threshold` 🔵 | float | `0.50` | high-soc 모드 발동 SOC 임계값 |
| `--motif-rh-fraction` 🔵 | float | `0.70` | RH-like/균일 분포 시드 비율 |
| `--motif-facing-fraction` 🔵 | float | `0.20` | Li-Li facing pair 시드 비율 |
| `--motif-random-fraction` 🔵 | float | `0.10` | 무작위 시드 비율 (세 비율의 합이 1이 되도록 정규화) |
| `--facing-mutation-prob` 🔵 | float | `0.20` | facing 돌연변이 비율 (새 개체 중) |
| `--ga-uniform-vacancy-mode` 🔵 | choice | `seed` | 균일 vacancy 분포 강제 방식 |
| `--uniform-vacancy-tolerance` 🔵 | int | `1` | 레이어간 vacancy 수 최대 허용 편차 |

**`--ga-motif-seeding` 선택지:**

| 값 | 설명 |
|----|------|
| `none` | 순수 무작위 초기 population |
| `high-soc` ✅ | SOC ≥ threshold 시 모티프 시딩 활성화 |
| `always` | 모든 SOC에서 모티프 시딩 |

**`--ga-uniform-vacancy-mode` 선택지:**

| 값 | 설명 |
|----|------|
| `none` | 균일 vacancy 미강제 |
| `seed` ✅ | 모티프 시드 생성 시에만 균일 vacancy 선호 |
| `hard` | 모든 GA 개체에 균일 vacancy 강제 (엄격) |

---

### 5.7 중복 제거 (Deduplication)

| 옵션 | 타입 | 기본값 | 설명 |
|------|------|--------|------|
| `--dedup-mode` 🔵 | choice | `symmetry` | 중복 구조 필터링 방식 |
| `--dedup-ltol` 🔵 | float | `0.2` | StructureMatcher 격자 길이 허용오차 |
| `--dedup-stol` 🔵 | float | `0.3` | StructureMatcher 사이트 허용오차 |
| `--dedup-angle-tol` 🔵 | float | `5.0` °| StructureMatcher 각도 허용오차 |

**`--dedup-mode` 선택지:**

| 값 | 설명 |
|----|------|
| `symmetry` ✅ | pymatgen StructureMatcher로 대칭 동치 제거 |
| `exact` | 정확히 동일한 분수좌표 집합만 제거 |
| `none` | 중복 제거 비활성화 |

---

### 5.8 에너지 백엔드 (MTP / Ewald)

| 옵션 | 타입 | 기본값 | 설명 |
|------|------|--------|------|
| `--energy-backend` 🔵 | choice | `mtp-relax` | 에너지 평가 방법 |
| `--mlp-bin` 🔵 | str | `mlp` | MLIP 실행파일 경로 또는 이름 |
| `--mlip-ini` 🔵 | str | `mlip.ini` | `mlp relax`에 사용할 ini 파일 |
| `--mtp-force-tolerance` 🔵 | float | 없음 | `--force-tolerance` 값 (미지정 시 mlp 기본값) |
| `--mtp-workers` 🔵 | int | `1` | 동시 `mlp relax` 프로세스 수 |
| `--mtp-batch-size` 🔵 | int | `0` | CFG 청크 크기 (0=SOC/세대 전체를 하나로) |
| `--mtp-keep-workdirs` 🔵 | flag | off | 임시 CFG/로그 보존 (디버깅용) |
| `--mtp-failed-energy` 🔵 | float | `1e100` | 완화 실패 시 패널티 에너지 |
| `--mtp-extra-args` 🔵 | str* | 없음 | `mlp relax`에 추가할 인자들 |

---

### 5.9 v5 I/O 최적화 옵션

GA에서 세대마다 대량으로 생성되던 파일을 줄이는 옵션들입니다.

| 옵션 | 기본값 | 설명 |
|------|--------|------|
| `--mtp-write-aggregate` ✅ | **on** | 세대별 `IN_*.cfg` 집계 파일 쓰기 |
| `--no-mtp-write-aggregate` | — | `IN_*.cfg` 생략 (대폭 I/O 절감) |
| `--mtp-single-log-per-soc` ✅ | **on** | SOC별 단일 로그 파일 (chunk별 분리 대신) |
| `--no-mtp-single-log-per-soc` | — | chunk별 개별 로그 (v4 동작) |
| `--mtp-flush-relaxed-per-soc` ✅ | **on** | SOC 완료 후 relaxed CFG 통합 1회 쓰기 |
| `--no-mtp-flush-relaxed-per-soc` | — | 세대마다 즉시 쓰기 (v4 동작) |

> **I/O 절감 효과 예시** (60 세대, 10 SOC 레벨):  
> - `--no-mtp-write-aggregate`: `IN_*.cfg` 파일 600개 → 0개  
> - `--mtp-single-log-per-soc`: 로그 파일 600개 → 10개  
> - `--mtp-flush-relaxed-per-soc`: `relaxed_*.cfg` 600개 → 10개

---

### 5.10 v5 GA 고급 옵션

#### 조기 종료 (Early Stopping)

| 옵션 | 타입 | 기본값 | 설명 |
|------|------|--------|------|
| `--ga-early-stop-gens` 🔵 | int | `0` (꺼짐) | N 세대 연속 개선 없으면 GA 조기 종료 |

#### 적응적 돌연변이 (Adaptive Mutation)

| 옵션 | 타입 | 기본값 | 설명 |
|------|------|--------|------|
| `--ga-adaptive-mutation` ✅ | flag | **on** | 정체 시 random/fullwipe 비율 자동 증가 |
| `--no-ga-adaptive-mutation` | flag | — | 적응적 돌연변이 비활성화 |
| `--ga-adaptive-after-gens` 🔵 | int | `5` | 몇 세대 정체 후 adaptive 발동할지 |
| `--ga-fullwipe-fraction` 🔵 | float | `0.05` | fullwipe(전체 Li 재생성) 비율 (정체 시) |

#### 다양성 필터 (Diversity Filter)

| 옵션 | 타입 | 기본값 | 설명 |
|------|------|--------|------|
| `--ga-min-symdiff` 🔵 | int | `0` (꺼짐) | 새 개체와 현 세대 멤버 간 최소 Hamming 거리 |

> `--ga-min-symdiff 2`: 기존 개체와 Li 사이트 2개 이상 달라야 population 추가  
> 조기 수렴 방지에 효과적. 값이 클수록 탐색 다양성↑ 속도↓

#### Elite 섭동 (Elite Perturbation)

| 옵션 | 타입 | 기본값 | 설명 |
|------|------|--------|------|
| `--ga-elite-perturb-fraction` 🔵 | float | `0.05` | 최우수 구조 1-2 사이트 교체 비율 (PSO small_perturb) |

---

### 5.11 출력 제어

| 옵션 | 타입 | 기본값 | 설명 |
|------|------|--------|------|
| `--no-poscar` 🔵 | flag | off | POSCAR 파일 생략 (전수 열거 모드) |
| `--sig-figs` 🔵 | int | `8` | POSCAR 좌표 유효 자릿수 |

---

## 6. 사용 예시

### 예시 1: 최소 실행 (Ewald 에너지, 전수 열거)

```bash
python enumerate_li_configs_Final_optimized_motif_v5_mtp.py \
    --structure LiC6.cif \
    --energy-backend ewald
```
- 2×2×2 슈퍼셀, AA 스태킹, 모든 SOC 레벨 열거
- 대칭 유일 구조만 선택, Ewald 에너지 계산
- 빠른 탐색에 적합

---

### 예시 2: 특정 SOC만 처리 (전수 열거)

```bash
python enumerate_li_configs_Final_optimized_motif_v5_mtp.py \
    --structure LiC6.cif \
    --energy-backend ewald \
    --soc-levels 12 16 20 24 \
    --output ./results_ewald
```

---

### 예시 3: GA 기본 실행 (Ewald 에너지)

```bash
python enumerate_li_configs_Final_optimized_motif_v5_mtp.py \
    --structure LiC6.cif \
    --energy-backend ewald \
    --ga \
    --ga-pop 120 \
    --ga-gens 60
```
- population 120, 60 세대 GA
- 모티프 시딩: SOC ≥ 0.50에서 자동 활성 (`high-soc` 기본값)
- 스태킹: AA only (SOC ≥ 0.50 보수적 정책 기본값)

---

### 예시 4: GA + 모든 스태킹 + 논문 기반 정책

```bash
python enumerate_li_configs_Final_optimized_motif_v5_mtp.py \
    --structure LiC6.cif \
    --energy-backend ewald \
    --ga \
    --stacking-sequences ALL \
    --soc-stacking-policy literature \
    --ga-motif-seeding always \
    --output ./results_ga_all_stacking
```
- 모든 A/B 스태킹 시퀀스 탐색
- SOC ≥ 0.40부터 AA-only (논문 Kim et al. 기준)
- 모든 SOC에서 모티프 시딩 활성

---

### 예시 5: MTP 에너지 기본 실행 ⭐ (권장)

```bash
python enumerate_li_configs_Final_optimized_motif_v5_mtp.py \
    --structure LiC6.cif \
    --mlip-ini /path/to/mlip.ini \
    --ga \
    --stacking-sequences ALL \
    --soc-stacking-policy literature \
    --ga-pop 120 \
    --ga-gens 60 \
    --output ./results_mtp
```
- MTP 에너지 백엔드 (기본값)
- v5 I/O 최적화 자동 적용 (기본값)
- 논문 기반 스태킹 정책

---

### 예시 6: MTP + 병렬 평가 (멀티 chunk)

```bash
python enumerate_li_configs_Final_optimized_motif_v5_mtp.py \
    --structure LiC6.cif \
    --mlip-ini mlip.ini \
    --ga \
    --mtp-workers 4 \
    --mtp-batch-size 30 \
    --output ./results_mtp_parallel
```
- 한 세대 120개 구조를 30개씩 4개 청크로 분할
- 4개 `mlp relax` 프로세스 병렬 실행
- **주의:** `mlp` 자체가 이미 멀티코어 사용 시 workers=1이 더 효율적

---

### 예시 7: v5 고급 GA + I/O 절감 ⭐ (대규모 계산 권장)

```bash
python enumerate_li_configs_Final_optimized_motif_v5_mtp.py \
    --structure LiC6.cif \
    --mlip-ini mlip.ini \
    --ga \
    --stacking-sequences ALL \
    --soc-stacking-policy literature \
    --ga-pop 150 \
    --ga-gens 80 \
    --ga-early-stop-gens 20 \
    --ga-min-symdiff 2 \
    --ga-fullwipe-fraction 0.08 \
    --ga-elite-perturb-fraction 0.08 \
    --no-mtp-write-aggregate \
    --output ./results_v5_optimized
```
- 조기 종료: 20 세대 개선 없으면 중단
- Hamming distance ≥ 2 다양성 필터
- fullwipe 8%, elite perturbation 8%
- `IN_*.cfg` 집계 파일 생략 (파일 수 대폭 감소)

---

### 예시 8: 특정 슈퍼셀 + 특정 스태킹 + 특정 SOC

```bash
python enumerate_li_configs_Final_optimized_motif_v5_mtp.py \
    --structure LiC6.cif \
    --supercell 3 3 2 \
    --stacking-sequences AAAA ABAA ABAB \
    --soc-levels 18 24 30 36 \
    --energy-backend ewald \
    --output ./results_3x3x2
```

---

### 예시 9: 완전한 디버그 실행 (파일 전부 보존)

```bash
python enumerate_li_configs_Final_optimized_motif_v5_mtp.py \
    --structure LiC6.cif \
    --mlip-ini mlip.ini \
    --ga \
    --ga-gens 10 \
    --ga-pop 30 \
    --mtp-keep-workdirs \
    --mtp-write-aggregate \
    --no-mtp-single-log-per-soc \
    --no-mtp-flush-relaxed-per-soc \
    --soc-levels 12 \
    --output ./debug_run
```
- 임시 CFG 파일 보존
- 세대별 aggregate 입력 파일 쓰기
- chunk별 로그 파일 분리
- 단일 SOC(n_Li=12)만 실행

---

### 예시 10: v4 완전 호환 모드 (v5 신기능 모두 끄기)

```bash
python enumerate_li_configs_Final_optimized_motif_v5_mtp.py \
    --structure LiC6.cif \
    --mlip-ini mlip.ini \
    --ga \
    --no-mtp-flush-relaxed-per-soc \
    --no-mtp-single-log-per-soc \
    --no-ga-adaptive-mutation \
    --ga-early-stop-gens 0 \
    --ga-min-symdiff 0 \
    --output ./results_v4_compat
```

---

## 7. 워크플로우 가이드

### 탐색 규모별 추천 설정

#### 소형 탐색 (테스트, 빠른 확인)
```bash
python enumerate_li_configs_Final_optimized_motif_v5_mtp.py \
    --structure LiC6.cif \
    --energy-backend ewald \
    --soc-levels 12 24 \
    --output ./test_run
```

#### 중형 탐색 (단일 스태킹, GA)
```bash
python enumerate_li_configs_Final_optimized_motif_v5_mtp.py \
    --structure LiC6.cif \
    --energy-backend ewald \
    --ga --ga-pop 100 --ga-gens 50 \
    --soc-stacking-policy conservative \
    --output ./medium_run
```

#### 대형 탐색 (모든 스태킹, MTP, 논문 재현)
```bash
python enumerate_li_configs_Final_optimized_motif_v5_mtp.py \
    --structure LiC6.cif \
    --mlip-ini mlip.ini \
    --ga --ga-pop 200 --ga-gens 100 \
    --stacking-sequences ALL \
    --soc-stacking-policy literature \
    --ga-motif-seeding always \
    --ga-early-stop-gens 25 \
    --ga-min-symdiff 2 \
    --no-mtp-write-aggregate \
    --output ./full_run
```

---

### SOC-의존 스태킹 정책 선택 가이드

```
SOC 구간    | 물리적 상태         | 권장 정책
-----------|--------------------|-----------
SOC < 0.20 | C-C dominant      | none (AB 허용)
SOC < 0.40 | Li-C competitive  | none 또는 conservative
SOC ≥ 0.40 | Li-Li dominant    | literature (AA-only)
SOC ≥ 0.50 | AA 완전 지배       | conservative 또는 literature
```

---

## 8. 성능 튜닝 가이드

### GA 수렴 속도 향상

| 상황 | 권장 조치 |
|------|----------|
| 수렴이 너무 느림 | `--ga-early-stop-gens 15` 설정 |
| 조기 수렴 의심 | `--ga-min-symdiff 2` 또는 `3` 설정 |
| 국소 최소에 갇힘 | `--ga-fullwipe-fraction 0.10` 증가 |
| 미세 탐색 강화 | `--ga-elite-perturb-fraction 0.10` 증가 |
| 고SOC 수렴 불량 | `--ga-motif-seeding always` + `--motif-rh-fraction 0.80` |

### MTP I/O 최적화

| 상황 | 권장 조치 |
|------|----------|
| 디스크 공간 부족 | `--no-mtp-write-aggregate` (파일 수 ~60배 감소) |
| 로그 파일 과다 | `--mtp-single-log-per-soc` (기본값, 이미 활성) |
| mlp 자체가 멀티코어 | `--mtp-workers 1` (기본값, 충돌 방지) |
| mlp 싱글코어 + 서버 멀티코어 | `--mtp-workers 4 --mtp-batch-size 30` |

### 메모리 관리

- `--ga-pop`이 너무 크면 메모리 증가 → 256 이하 권장
- MTP `relaxed_structure_cache`가 누적됨 → 매우 많은 SOC 레벨 실행 시 주의
- `--soc-levels`로 필요한 SOC만 지정하여 실행 범위 제한

---

## 9. 출력 파일 설명

### `enumeration_summary.json`

전체 실행의 JSON 요약 파일. 주요 필드:

```json
{
  "structure_file": "LiC6.cif",
  "supercell": [2, 2, 2],
  "ga_enabled": true,
  "energy_backend": "mtp-relax",
  "soc_stacking_policy": "literature",
  "stable_by_soc": [
    {
      "n_Li": 24,
      "SOC": 0.5,
      "method": "GA",
      "best_stacking_sequence": "AAAA",
      "best_mtp_relaxed_energy_eV": -123.456789,
      "best_occupied_Li_sites": [0, 3, 6, ...]
    }
  ],
  "soc_levels": [
    {
      "n_Li": 24,
      "SOC": 0.5,
      "n_ga_survivors": 45,
      "ga_motif_seeding_active": true,
      "n_ga_motif_rh_like_seeds": 84,
      "n_ga_motif_facing_seeds": 24,
      "ga_early_stop_gen": 35,
      "n_ga_elite_perturb_candidates": 312,
      "n_ga_fullwipe_candidates": 58,
      "top_stable": [ ... ]
    }
  ]
}
```

### CIF 파일 (`top_stable_cif/`)

- MTP 완화된 구조 (가능할 경우)
- 파일명: `rank001_AAAA_nLi24_SOC0.5000_MTP_E-123.456789eV.cif`
- 최대 20개 (symmetry-unique, 에너지 오름차순)

### GA 수렴 그래프 (`ga_convergence_nLi{N}.png`)

세대별 각 개체의 에너지 분포를 시각화. 색상 구분:
- 🟢 초록: Survivor (에너지 임계 이내)
- 🔵 파랑: Random
- 🔴 빨강: Vacancy swap 돌연변이
- 🟠 주황: Li layer swap 돌연변이
- 🟣 보라: Stacking flip 돌연변이
- 🟤 갈색: Facing pair 돌연변이
- 🩵 시안: RH-like 시드
- 💜 마젠타: Facing 시드
- 🟩 진초록: Elite perturbation (v5 신규)
- ⬜ 회색: Fullwipe restart (v5 신규)

---

## 10. 자주 묻는 질문

**Q. GA와 전수 열거 중 무엇을 써야 하나요?**  
A. Li 수가 많아 C(N,k) > 10⁶이면 GA 권장. `--ga` 플래그를 추가하면 됩니다.

**Q. `mlip.ini` 파일은 어떻게 준비하나요?**  
A. 기존 MTP 포텐셜 파일(`.mtp`)과 함께 `mlip.ini`를 작성해야 합니다. MLIP 소프트웨어 문서를 참고하세요.

**Q. `--mtp-workers`를 늘려도 빠르지 않습니다.**  
A. `mlp relax`가 이미 내부적으로 멀티코어를 사용하는 경우 workers=1이 최적입니다. `mlp relax`가 싱글코어이고 서버에 여유 코어가 있을 때만 workers > 1이 효과적입니다.

**Q. GA가 너무 일찍 수렴합니다 (에너지가 초반에 멈춤).**  
A. `--ga-min-symdiff 2`로 다양성 필터를 추가하고, `--ga-fullwipe-fraction 0.10`으로 fullwipe 비율을 높이세요.

**Q. 특정 SOC에서 `NoValidGACandidateError`가 발생합니다.**  
A. `--ga-min-li-li-same-layer` 값을 낮추거나 (예: 3.5), `--no-ga-expand-li-sites`로 확장 사이트를 끄거나, 해당 SOC를 `--soc-levels`에서 제외하세요.

**Q. 논문 결과를 재현하려면?**  
A. `--soc-stacking-policy literature --ga-motif-seeding always --stacking-sequences ALL` 조합을 사용하세요. 랜덤 시드는 `--ga-seed 42` (기본값)입니다.

**Q. 디스크 사용량이 너무 큽니다.**  
A. `--no-mtp-write-aggregate`를 추가하면 세대별 `IN_*.cfg` 파일이 생략되어 디스크 사용량이 크게 줄어듭니다. `--mtp-flush-relaxed-per-soc`(기본 활성)도 완화 CFG 파일을 SOC당 1개로 통합합니다.

---

*마지막 수정: 2026-05-17 | 스크립트 버전: v5_mtp*
