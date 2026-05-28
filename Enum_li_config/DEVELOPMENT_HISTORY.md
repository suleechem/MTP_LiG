# enumerate_li_configs — 개발 히스토리

> **대상 코드:** `enumerate_li_configs_Final_optimized_motif_v5_mtp.py` 및 전신 버전들  
> **목적:** 흑연 삽입 화합물(GIC)에서 Li-vacancy 배열을 열거하고 가장 안정한 구조를 탐색  
> **기반 논문:** Kim et al., *ACS Nano* — "Revealing Li Staging Reactions in Graphite via a GA Coupled with MLIP"

---

## 버전 계보

```
enumerate_li_configs_Final_v0a.py               (2,540 lines)
        │  ↓ v0b: JSON 인코더, 스태킹/dedup 체계 개편
enumerate_li_configs_Final_optimized_v0b.py     (2,554 lines)
        │  ↓ v1: SOC 정책 + 모티프 시딩 GA 신설
enumerate_li_configs_Final_optimized_motif_v1_Ewald.py    (3,029 lines)
        │  ↓ v2: GA 내부 루프 알고리즘 최적화 (O(N²) → O(N))
enumerate_li_configs_Final_optimized_motif_v2_Ewald.py    (3,052 lines)
        │  ↓ v3: 멀티프로세싱 병렬 에너지 평가
enumerate_li_configs_Final_optimized_motif_v3_Ewald_parallel.py  (3,201 lines)
        │  ↓ v4: 에너지 백엔드를 MTP(MLIP)로 교체
enumerate_li_configs_Final_optimized_motif_v4_mtp.py      (3,569 lines)
        │  ↓ v5: I/O 최적화 + GA 고급 연산자
enumerate_li_configs_Final_optimized_motif_v5_mtp.py      (3,866 lines)
```

---

## v0a — 초기 버전: 핵심 기능 확립

**파일:** `enumerate_li_configs_Final_v0a.py` (2,540 lines)

### 배경 및 목적
GIC 슈퍼셀에서 SOC(State of Charge)별 Li-vacancy 배열을 열거하는 코드의 최초 완성본.
전수 열거와 GA 탐색 두 가지 모드를 모두 갖춘 기반 프레임워크.

### 핵심 기능

#### 구조 처리
- `load_and_build_supercell()` — CIF/POSCAR 읽기, a×b×c 슈퍼셀 생성
- `extract_sublattices()` — Li / C 부분격자 분리
- `detect_layers_by_z()`, `detect_li_layers()` — c축 방향 층 자동 감지
- `auto_detect_b_layer_shift()` — AB 스태킹 B 레이어 면내 변위 자동 계산 (Li 사이트 → 헥사곤 중심 활용)
- `make_carbon_stacking_fracs()` — 스태킹 시퀀스(AAAA, ABAB 등)별 탄소 분수좌표 생성
- `expand_li_sites_with_hollow_centers()` — 흑연 헥사곤 중심 빈 사이트 자동 감지 및 GA Li 후보 풀 확장
- `build_same_layer_compatibility_matrix()` — 같은 층 내 Li-Li 최소 거리 제약 행렬 사전 계산

#### 열거 엔진
- `enumerate_symmetric()` — 대칭 유일 대표 구조 열거 (canonical form)
- `enumerate_all()` — C(N,k) 전수 열거
- `get_li_site_permutations()` — SpacegroupAnalyzer로 Li 사이트 대칭 순열 추출
- `generate_all_stacking_sequences()` — 첫 층 A 고정, 나머지 A/B 모든 조합 생성 (`ALL`)

#### 에너지 평가
- `compute_ewald_energy()` — Li=+1, C=−n_Li/n_C 산화수로 Ewald 전기정적 에너지 계산
- `EwaldPrecomputed` 클래스 — 스태킹 시퀀스별 Ewald 행렬 사전 계산  
  → 기존 O((N_Li+N_C)²)에서 설정 당 O(n_li²)로 대폭 가속

#### GA (유전 알고리즘)
- `GAOptimizer` 클래스 — 돌연변이-선택 기반 GA (교차 없음)
  - 연산자: `_vacancy_swap`, `_layer_swap`, `_stacking_flip`
  - `_compatible_add_sites()` — 쌍 제약 만족하는 추가 가능 사이트 반환
  - `_enumerate_valid_occupied_pool()` — 하드코어 제약 하 유효 구성 역추적 열거
  - 에너지 캐시, 세대별 survivor 선택
  - 수렴 플롯 (`matplotlib`)

#### 중복 제거 및 출력
- `select_top_symmetry_unique_candidates()` — 에너지 정렬 → 좌표 해시 → StructureMatcher 순으로 상위 20개 symmetry-unique 구조 선택
- `dedupe_stacking_candidates_by_structure()` — 전체 중복 제거
- CIF, POSCAR, JSON 요약 출력

---

## v0b — 안정화: JSON 인코더 + 스태킹/dedup 체계 정리

**파일:** `enumerate_li_configs_Final_optimized_v0b.py` (2,554 lines, +14)

### 변경 배경
v0a 실행 결과 JSON 저장 시 `numpy` 타입(`np.integer`, `np.floating`, `np.bool_`, `np.ndarray`)이 직렬화 불가 오류 발생.
스태킹-aware 중복 제거 로직과 JSON 요약 필드 구조 개편.

### 주요 변경

| 구분 | 내용 |
|------|------|
| **`_NumpyEncoder`** 신설 | numpy 스칼라/배열 타입을 JSON 직렬화 가능하게 처리. `json.dump(..., cls=_NumpyEncoder)` |
| **`perms_by_stacking`** | 스태킹 시퀀스별 Li 대칭 순열을 딕셔너리로 관리 (단일 `perms` → `perms_by_stacking[seq]`) |
| **`dedup_stats`** | 중복 제거 통계(`n_exact_duplicates_removed`, `n_symmetry_duplicates_skipped_for_top20`)를 별도 dict로 정리 |
| **`stacking_policy_info`** | 스태킹 정책 정보를 `soc_summary`에 포함 |
| **`n_candidate_key_duplicates_removed`** | GA 세대 내 중복 키 통계 집계 |

---

## v1 — 물리 기반 스태킹 정책 + 모티프 시딩 GA 신설

**파일:** `enumerate_li_configs_Final_optimized_motif_v1_Ewald.py` (3,029 lines, +475)

### 변경 배경
논문(Kim et al., *ACS Nano*)의 핵심 물리 결과 반영:

| SOC 구간 | 지배 상호작용 | 안정 스태킹 |
|---------|-------------|-----------|
| < 20% | C–C dominant | AB 선호 |
| 20–40% | Li–C competitive | AA/AB 경쟁 |
| ≥ 40% | Li–Li dominant | **AA only** |

고SOC에서 GA 초기 population이 무작위 생성으로 인해 물리적으로 비현실적인 구성에 집중되는 문제 → 모티프 시딩으로 해결.

### 새로운 함수

| 함수 | 역할 |
|------|------|
| `get_aa_only_soc_threshold()` | 정책(`none/conservative/literature`)에서 SOC 임계값 반환 |
| `active_stacking_sequences_for_soc()` | 현재 SOC에서 허용되는 스태킹 시퀀스 필터링 |
| `normalize_fraction_triplet()` | 세 비율 합이 1이 되도록 정규화 |
| `motif_seeding_active_for_soc()` | `none/high-soc/always` 모드에서 모티프 시딩 활성 여부 판단 |

### SOC 기반 스태킹 정책 (`--soc-stacking-policy`)

```
none         → 모든 SOC에서 모든 스태킹 시퀀스 사용
conservative → SOC ≥ 0.50부터 AA-only
literature   → SOC ≥ 0.40부터 AA-only (논문 Kim et al. 기준)
--aa-only-soc-threshold → 직접 임계값 지정
```

### 모티프 시딩 GA (`--ga-motif-seeding`)

GAOptimizer에 3종 시드 타입 추가:

| 시드 타입 | 생성 방법 | 인코딩하는 물리 |
|----------|----------|----------------|
| **RH-like / alternating** (`_random_balanced_occupied`) | 층별 균등 분배 | 고SOC 균일 vacancy 분포 (RH 상) |
| **Facing-enhanced** (`_facing_enhanced_occupied`) | 동일 xy 열 우선 배치 | Li-Li 층간 facing pair 선호 |
| **Random** | 순수 무작위 | 탐색 다양성 유지 |

```
--ga-motif-seeding  none / high-soc / always
--motif-seed-soc-threshold  0.50
--motif-rh-fraction         0.70
--motif-facing-fraction     0.20
--motif-random-fraction     0.10
--facing-mutation-prob      0.20
--ga-uniform-vacancy-mode   none / seed / hard
--uniform-vacancy-tolerance 1
```

새 GA 연산자: **facing mutation** — 기존 survivor에서 facing 구성으로 변환.

새 GA 통계 필드:
- `n_motif_rh_seeded`, `n_motif_facing_seeded`, `n_motif_random_seeded`
- `n_facing_mutations_generated`, `n_uniform_vacancy_rejections`
- `_xy_columns_cache` — `_sites_by_xy_column()` 결과 지연 캐시

---

## v2 — GA 내부 루프 알고리즘 최적화

**파일:** `enumerate_li_configs_Final_optimized_motif_v2_Ewald.py` (3,052 lines, +23)

### 변경 배경
v1에서 고SOC 실행 시 GA 속도 병목 분석 결과:
- `_random_balanced_occupied()` 내부 층별 체크가 O(N²) list comprehension
- `_facing_enhanced_occupied()`, `_random_occupied()` 등 모든 greedy 루프에서 `_is_valid_occupied()` (O(k²) 행렬 검사)를 매 사이트 추가마다 호출

### 핵심 최적화

| 위치 | 변경 전 | 변경 후 | 효과 |
|------|--------|--------|------|
| `_random_balanced_occupied()` | `len([x for x in chosen if x in layer_sites])` | `layer_count` 정수 카운터 | O(N²) → O(1) per site |
| `_random_occupied()` 내부 | `tuple(sorted()) + _is_valid_occupied()` | `_can_add_site(site, chosen)` | O(k²) → O(k) |
| `_facing_enhanced_occupied()` | `_is_valid_occupied()` | `_can_add_site()` | O(k²) → O(k) |
| `_layer_swap()` 대상 사이트 필터 | `_is_valid_occupied()` | `_can_add_site()` | O(k²) → O(k) |
| `_vacancy_swap()` | `_compatible_add_sites()` 후 `_is_valid_occupied()` 재확인 | 두 번째 호출 제거 | 중복 검사 제거 |
| `_compatible_add_sites()` vacant_mask | Python for-loop | `np.fromiter()` 벡터화 | 벡터화 |
| `normalize_fraction_triplet()` | `soc_summary`에서 3회 호출 | 1회 호출 후 `dict(zip())` | 중복 제거 |

**신설 헬퍼:**
```python
def _can_add_site(self, site: int, existing: List[int]) -> bool:
    """O(k) 단일 사이트 호환성 검사. compatible_pair_matrix[site, occ] 행만 검사."""
    if self.compatible_pair_matrix is None or not existing:
        return True
    occ = np.asarray(existing, dtype=int)
    return bool(np.all(self.compatible_pair_matrix[site, occ]))
```

**버그 수정:**
- `_facing_enhanced_occupied()` 내부 fallback 루프에서 `n_motif_random_seeded` 카운터 잘못 증가 → 제거
- `optimizer` 변수 스코프: `if 'optimizer' in locals()` → `optimizer: Optional[GAOptimizer] = None` 명시적 선언

---

## v3 — 멀티프로세싱 병렬 에너지 평가

**파일:** `enumerate_li_configs_Final_optimized_motif_v3_Ewald_parallel.py` (3,201 lines, +149)

### 변경 배경
전수 열거 경로에서 수천~수만 개 구조의 Ewald 에너지를 직렬로 계산하는 병목.
대형 슈퍼셀(3×3×2 등)에서 SOC별 수만 개 후보 처리 시 시간 과다.

### 새로운 함수

| 함수 | 역할 |
|------|------|
| `resolve_worker_count(requested)` | `--workers` 값 정규화 (0 또는 음수 → 전체 CPU 코어) |
| `_init_parallel_energy_worker(ep_dict)` | multiprocessing worker 초기화 (EwaldPrecomputed 딕셔너리 전달) |
| `_parallel_candidate_energy(candidate)` | 단일 후보 에너지 계산 (worker 프로세스 내) |
| `evaluate_candidates_precomputed(...)` | 직렬/병렬 모드 분기, `Pool.imap`으로 청크 단위 평가 |

### 신설 CLI 인수

```
--workers / -j         병렬 프로세스 수 (기본 1, 0=자동)
--parallel-chunksize   Pool.imap 청크 크기 (기본 64)
```

### GAOptimizer 병렬화

GA 내 `_evaluate_population()` 메서드 추가:
- population 일괄 에너지 평가를 `evaluate_candidates_precomputed()`로 위임
- `--workers > 1`이고 `EwaldPrecomputed` 사용 시 세대 내 병렬 평가 적용

### 아키텍처 패턴

```
multiprocessing.Pool(processes=n_workers,
                     initializer=_init_parallel_energy_worker,
                     initargs=(ewald_dict,))
pool.imap(_parallel_candidate_energy, candidates, chunksize)
```

**주의:** `EwaldPrecomputed` 객체가 각 worker 프로세스에 전역으로 공유 (pickle 오버헤드 최소화).

---

## v4 — MLIP/MTP 에너지 백엔드 도입

**파일:** `enumerate_li_configs_Final_optimized_motif_v4_mtp.py` (3,569 lines, +368)

### 변경 배경
Ewald 전기정적 에너지는 Li-Li 거리 제약 외 국소 구조 효과를 포착 못함.
MLIP(Machine Learning Interatomic Potential)로 `mlp relax`를 통해 실제 구조 완화 에너지 계산 필요.
PSO 코드(`PSO_opt_parallel_var.py`)의 병렬 평가 패턴 참조.

### 새로운 함수

| 함수 | 역할 |
|------|------|
| `write_mtp_cfg_structure()` | pymatgen Structure → MLIP CFG 형식 직렬화 |
| `split_cfg_blocks_from_text()` | CFG 텍스트를 `BEGIN_CFG`/`END_CFG` 블록으로 분리 |
| `split_cfg_blocks_file()` | CFG 파일에서 블록 파싱 |
| `parse_mtp_cfg_block()` | 완화된 CFG 블록 → Structure + energy 파싱 |
| `_candidate_feature()` | 후보 식별 메타데이터 문자열 생성 |
| `write_soc_top_cifs_mtp()` | MTP 완화 구조 우선 사용하는 CIF 출력 |

### `MTPRelaxationEvaluator` 클래스 (신설)

```
역할: GA와 전수 열거 양 경로의 에너지 평가를 `mlp relax` subprocess로 대체

핵심 동작:
  1. 후보 구조들을 IN_{label}.cfg에 CFG 형식으로 기록
  2. energy_cache에 없는 구조만 골라 to-relax_{ci}.cfg로 청크 분할
  3. 각 청크에 `mlp relax mlip.ini --cfg-filename=... --save-relaxed=...` 실행
  4. --mtp-workers > 1이면 ThreadPoolExecutor로 청크 병렬 실행
  5. 완화된 CFG 파싱 → energy_cache, relaxed_structure_cache, relaxed_block_cache 업데이트
  6. 모든 청크 결과를 relaxed_{label}.cfg로 통합 기록

3-tier 캐시:
  energy_cache           : (stacking, occupied) → float
  relaxed_structure_cache: (stacking, occupied) → Structure
  relaxed_block_cache    : (stacking, occupied) → str (원본 CFG 블록)
```

### 새로운 CLI 인수

```
--energy-backend      mtp-relax / ewald   (기본: mtp-relax)
--mlp-bin             mlp 실행파일 경로
--mlip-ini            mlip.ini 파일
--mtp-force-tolerance --force-tolerance 값
--mtp-workers         동시 mlp relax 프로세스 수
--mtp-batch-size      CFG 청크 크기 (0=전체)
--mtp-keep-workdirs   임시 CFG 파일 보존
--mtp-failed-energy   완화 실패 시 패널티 에너지
--mtp-extra-args      mlp relax 추가 인수
```

### 에너지 평가 경로 분기

```
--energy-backend ewald   → EwaldPrecomputed (v3까지 동작)
--energy-backend mtp-relax → MTPRelaxationEvaluator
                             └ GA:        _energies_for_population() → evaluate_candidates()
                             └ 전수 열거: evaluate_candidates() 직접 호출
```

### 파일 출력 구조 (세대별)

```
SOC_{soc}_nLi{n}/
  mtp_cfg/
    IN_ga_gen000.cfg        ← 세대 전체 입력 (집계)
    relaxed_ga_gen000.cfg   ← 세대 전체 완화 결과
  mtp_logs/
    relax_ga_gen000_chunk0000.txt  ← 청크별 mlp 로그
  mtp_work/ga_gen000/       ← 임시 청크 파일 (실행 후 삭제)
```

---

## v5 — I/O 최적화 + GA 고급 연산자 (PSO 패턴 적용)

**파일:** `enumerate_li_configs_Final_optimized_motif_v5_mtp.py` (3,866 lines, +297)

### 변경 배경
v4에서 확인된 두 가지 문제:

1. **I/O 과다**: 60 세대 × 10 SOC = 600개 `IN_*.cfg` + 600개 `relaxed_*.cfg` + 600개 로그 파일 생성
2. **GA 조기 수렴**: 순수 돌연변이-선택만으로는 정체(stagnation) 극복 어려움

PSO 코드(`PSO_opt_parallel_var.py`)의 다양성 유지 전략 참조:
- `small_perturb()`: elite 구조 1-2 사이트 교체
- `fullwipe`: 전체 Li 재생성
- `min_symdiff`: Hamming 거리 기반 다양성 필터

### I/O 최적화 (`MTPRelaxationEvaluator` 개선)

| 기능 | v4 동작 | v5 기본 동작 | 절감 효과 |
|------|--------|------------|---------|
| 집계 입력 파일 | 매 세대 `IN_*.cfg` 기록 | 선택적 (기본 on, `--no-mtp-write-aggregate`로 끄기) | 파일 수 ×60 감소 |
| 완화 결과 파일 | 세대마다 즉시 기록 | SOC 완료 후 1회 통합 기록 | 파일 수 ×60 감소 |
| 로그 파일 | 청크×세대별 분리 | SOC별 단일 파일 | 파일 수 ×60 감소 |

**신설 메서드:**
```python
def flush_relaxed_combined(self, soc_dir, label):
    """SOC 완료 후 _pending_relaxed_blocks 버퍼를 relaxed_{label}_all.cfg로 1회 기록."""
```

**인메모리 버퍼:**
```python
self._pending_relaxed_blocks: List[str] = []   # 세대별 완화 블록 누적
self._pending_relaxed_label: Optional[str] = None
```

### GA 고급 연산자

#### 1. 조기 종료 (Early Stopping)
```
--ga-early-stop-gens N   N 세대 연속 개선 없으면 GA 중단 (0=꺼짐)
```
- `n_stagnant_gens` 카운터로 정체 세대 추적
- 조기 종료 시 마지막 세대 survivor를 `all_survivors`에 추가 후 종료

#### 2. 적응적 돌연변이 (Adaptive Mutation)
```
--ga-adaptive-mutation         (기본 on)
--ga-adaptive-after-gens 5     5 세대 정체 후 발동
--ga-fullwipe-fraction 0.05    fullwipe 비율
```
- 정체 깊이에 따라 선형 증폭:
  ```
  effective_frac_random = min(0.60, base + 0.20 × ramp)
  effective_fullwipe    = min(0.20, base + 0.10 × ramp)
  ```

#### 3. Fullwipe 연산자
```
--ga-fullwipe-fraction 0.05
```
PSO의 `fullwipe` 전략에서 차용. 정체된 population에 완전히 새로운 구성 주입 → 국소 최솟값 탈출.

#### 4. Elite 섭동 (Elite Perturbation)
```
--ga-elite-perturb-fraction 0.05
```
PSO의 `small_perturb()`에서 차용.
```python
def _elite_perturb(self, candidate, n_swaps=1):
    """최우수 구조에서 1-2개 사이트 교체 → 미세 국소 탐색."""
```

#### 5. 다양성 필터 (min_symdiff)
```
--ga-min-symdiff 2
```
PSO의 `min_symdiff` 패턴에서 차용. 새 개체와 현 세대 모든 멤버 간 Hamming 거리가 임계값 미만이면 거부.
```python
def _hamming_distance(self, a, b):
    """Li 사이트 집합의 대칭 차이 크기의 절반 = 실제로 다른 사이트 수."""
    return len(set(a).symmetric_difference(set(b))) // 2
```

#### 6. `_add()` 클로저 (코드 정리)
GA 메인 루프 내 `_append_unique_candidate()` + occupancy 추적을 단일 클로저로 캡슐화:
```python
def _add(candidate, ctype):
    added = self._append_unique_candidate(..., existing_occupancies=new_occupancies)
    if added and new_occupancies is not None:
        new_occupancies.append(candidate[1])
    return added
```

### 수렴 플롯 추가 타입
- `elite_perturb` (진초록)
- `fullwipe` (회색)

### 새로운 CLI 인수 (총 14개 추가)

**I/O 제어 (3쌍):**
```
--mtp-write-aggregate / --no-mtp-write-aggregate
--mtp-single-log-per-soc / --no-mtp-single-log-per-soc
--mtp-flush-relaxed-per-soc / --no-mtp-flush-relaxed-per-soc
```

**GA 고급 (6개):**
```
--ga-early-stop-gens
--ga-adaptive-mutation / --no-ga-adaptive-mutation
--ga-adaptive-after-gens
--ga-fullwipe-fraction
--ga-min-symdiff
--ga-elite-perturb-fraction
```

---

## 전체 기능 누적 요약

| 기능 범주 | v0a | v0b | v1 | v2 | v3 | v4 | v5 |
|----------|:---:|:---:|:--:|:--:|:--:|:--:|:--:|
| **열거 엔진** (대칭/전수) | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ |
| **스태킹 시퀀스** (AA/AB/ALL) | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ |
| **GA** (vacancy_swap, layer_swap, stacking_flip) | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ |
| **EwaldPrecomputed** (O(n²) 고속) | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ |
| **Li 사이트 확장** (hollow centers) | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ |
| **StructureMatcher dedup** | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ |
| **_NumpyEncoder** (JSON) | — | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ |
| **perms_by_stacking** (스태킹별 대칭) | — | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ |
| **SOC 스태킹 정책** (none/conservative/literature) | — | — | ✅ | ✅ | ✅ | ✅ | ✅ |
| **모티프 시딩 GA** (RH-like, facing) | — | — | ✅ | ✅ | ✅ | ✅ | ✅ |
| **facing mutation** 연산자 | — | — | ✅ | ✅ | ✅ | ✅ | ✅ |
| **균일 vacancy 강제** (seed/hard) | — | — | ✅ | ✅ | ✅ | ✅ | ✅ |
| **`_can_add_site()`** O(k) 최적화 | — | — | — | ✅ | ✅ | ✅ | ✅ |
| **`_random_balanced_occupied()`** O(N²)→O(N) | — | — | — | ✅ | ✅ | ✅ | ✅ |
| **`_xy_columns_cache`** 지연 캐시 | — | — | — | ✅ | ✅ | ✅ | ✅ |
| **멀티프로세싱** (Pool.imap, Ewald 병렬) | — | — | — | — | ✅ | ✅ | ✅ |
| **MTP/MLIP 에너지 백엔드** | — | — | — | — | — | ✅ | ✅ |
| **MTPRelaxationEvaluator** (CFG I/O, cache) | — | — | — | — | — | ✅ | ✅ |
| **ThreadPoolExecutor** (청크 병렬) | — | — | — | — | — | ✅ | ✅ |
| **조기 종료** (early stopping) | — | — | — | — | — | — | ✅ |
| **적응적 돌연변이** (adaptive mutation) | — | — | — | — | — | — | ✅ |
| **Fullwipe 연산자** | — | — | — | — | — | — | ✅ |
| **Elite 섭동** (small_perturb) | — | — | — | — | — | — | ✅ |
| **min_symdiff 다양성 필터** | — | — | — | — | — | — | ✅ |
| **I/O 지연 플러시** (SOC별 통합 기록) | — | — | — | — | — | — | ✅ |

---

## 알고리즘 복잡도 개선 이력

| 연산 | v0a | v2 이후 |
|------|-----|--------|
| Ewald 에너지 (전체 구조) | O((N_Li+N_C)²) per config | O(n_li²) via EwaldPrecomputed |
| 사이트 추가 호환성 검사 | O(k²) (전체 submatrix) | O(k) via `_can_add_site` |
| `_random_balanced_occupied` 층별 카운트 | O(N²) list comprehension | O(1) counter |
| `_compatible_add_sites` vacant mask | O(N) Python loop | O(N) numpy 벡터화 |
| 에너지 평가 (전수 열거) | 직렬 | v3: multiprocessing Pool / v4: ThreadPoolExecutor (청크) |

---

## 외부 참조 및 영향

| 참조 소스 | 도입 버전 | 도입된 개념 |
|----------|---------|-----------|
| Kim et al., *ACS Nano* 논문 | v1 | SOC 임계값(0.40/0.50), AA-only 체제, RH-like/facing 모티프 |
| PSO_opt_parallel_var.py | v4, v5 | ThreadPoolExecutor 청크 병렬, small_perturb, fullwipe, min_symdiff |

---

*작성일: 2026-05-17*
