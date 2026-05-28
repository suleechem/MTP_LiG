# `enumerate_li_configs_Final_optimized_motif.py` 사용 인스트럭션

## 1. 코드 목적

이 스크립트는 graphite intercalation compound 계열 구조에서 Li/vacancy configuration과 carbon-layer stacking sequence를 함께 고려하여, 각 SOC(State of Charge)별로 Ewald electrostatic energy 기준의 안정 구조를 탐색하는 도구입니다.

최종 motif 버전은 다음 기능을 포함합니다.

- CIF 또는 POSCAR 구조 파일 입력
- supercell 생성
- SOC별 Li/vacancy configuration 생성
- 기본 symmetry-aware enumeration
- `--no-symmetry` brute-force enumeration
- `--ga` 지정 시 genetic algorithm 방식 탐색
- carbon stacking sequence 고려
  - 예: `AA`, `AB`
  - 4층 예: `AAAA`, `ABAA`, `AABA`, `ABAB`
- `--ga` 없이도 `--stacking-sequences`를 full enumeration에 적용
- AB stacking을 graphene hexagon center와 carbon vertex registry 기준으로 정의
- GA에서 입력 Li site뿐 아니라 graphene hollow site를 Li candidate site로 확장
- 같은 Li layer 내 Li-Li 최소 거리 constraint 적용
- high SOC에서 valid configuration이 없으면 해당 SOC만 skip
- SOC-dependent stacking pruning
  - high SOC에서 all-A stacking만 사용하도록 search space 축소
- high-SOC motif-biased GA
  - RH-like / Alternating template seed
  - facing-pair enhanced seed
  - random seed
  - facing-pair mutation 강화
  - same-layer vacancy uniformity bias 또는 hard constraint
- 중복 구조 제거
  - hollow site 확장 단계: 좌표 기준 Li site 중복 제거
  - GA candidate 생성 단계: candidate key 및 coordinate hash 중복 제거
  - 최종 CIF 선택 단계: StructureMatcher 기반 symmetry-equivalence 중복 제거
- 각 SOC별 symmetry-unique top 20 CIF 저장
- 각 SOC별 최안정 stacking sequence, Li site, Ewald energy 출력
- 최종 `enumeration_summary.json` 생성

---

## 2. 실행 환경

필수 package:

```bash
pip install numpy pymatgen
```

선택 package:

```bash
pip install matplotlib tqdm
```

- `matplotlib`: GA convergence plot 저장
- `tqdm`: 긴 enumeration에서 progress bar 표시
- Python 3.8 호환 방식입니다.
- `argparse.BooleanOptionalAction`을 사용하지 않습니다.

---

## 3. 기본 실행

### 3.1 기본 symmetry-aware enumeration

```bash
python enumerate_li_configs_Final_optimized_motif.py \
  --structure LiC6.cif
```

동작:

- `--ga`가 없으므로 full enumeration 사용
- 기본적으로 symmetry-aware representative만 평가
- stacking sequence는 all-A sequence만 사용
  - carbon layer 2개: `AA`
  - carbon layer 4개: `AAAA`
- 각 SOC별 Ewald energy 기준 top 20 CIF 저장

---

### 3.2 no-symmetry brute-force enumeration

```bash
python enumerate_li_configs_Final_optimized_motif.py \
  --structure LiC6.cif \
  --no-symmetry
```

동작:

- symmetry reduction 없이 모든 `C(N,k)` Li 조합 평가
- stacking sequence를 여러 개 지정하면 계산량이 stacking sequence 수만큼 증가합니다.

---

### 3.3 GA 실행

```bash
python enumerate_li_configs_Final_optimized_motif.py \
  --structure LiC6.cif \
  --ga
```

동작:

- 모든 SOC level에 GA 적용
- GA는 자동으로 켜지지 않습니다.
- 반드시 `--ga`를 명시해야 합니다.
- GA 사용 시 hollow Li site 확장이 기본적으로 켜집니다.
- 같은 Li layer 내 Li-Li 최소 거리 constraint가 적용됩니다.
- 기본 설정에서는 SOC 50% 이상에서 high-SOC motif-biased GA가 작동합니다.

---

### 3.4 특정 SOC만 계산

```bash
python enumerate_li_configs_Final_optimized_motif.py \
  --structure LiC6.cif \
  --ga \
  --soc-levels 2 3 4 5 6
```

`--soc-levels`의 숫자는 SOC fraction이 아니라 점유 Li 개수입니다.

예를 들어 입력 구조 기준 Li site가 16개이고 `--soc-levels 10`이면:

```text
SOC = 10 / 16 = 0.625
```

GA에서 hollow site 확장이 켜져 있으면 실제 탐색에 쓰이는 Li candidate site 수는 입력 Li site 수보다 클 수 있습니다. 그러나 SOC 표기는 입력 구조의 기준 Li 수를 따라 출력됩니다.

---

## 4. Stacking sequence

### 4.1 기본 stacking

기본 stacking은 모든 carbon layer가 A registry인 구조입니다.

```text
2 carbon layers → AA
4 carbon layers → AAAA
```

---

### 4.2 모든 A/B stacking sequence 생성

```bash
python enumerate_li_configs_Final_optimized_motif.py \
  --structure LiC6.cif \
  --stacking-sequences ALL
```

동작:

- carbon layer 수 자동 검출
- 첫 번째 layer는 A로 고정
- 나머지 layer는 A/B 조합을 모두 생성

4 carbon layer 예:

```text
AAAA
AAAB
AABA
AABB
ABAA
ABAB
ABBA
ABBB
```

첫 layer를 A로 고정하는 이유는 전체 in-plane translation에 해당하는 중복을 줄이기 위해서입니다.

---

### 4.3 특정 stacking sequence만 지정

```bash
python enumerate_li_configs_Final_optimized_motif.py \
  --structure LiC6.cif \
  --stacking-sequences AAAA ABAA AABA ABAB
```

지정한 sequence만 고려합니다.

---

### 4.4 GA와 stacking sequence 동시 사용

```bash
python enumerate_li_configs_Final_optimized_motif.py \
  --structure LiC6.cif \
  --ga \
  --stacking-sequences ALL
```

GA candidate는 다음 형태입니다.

```python
candidate = (stacking_sequence, occupied_Li_sites)
```

예:

```python
("AA", (0, 3, 7))
("AB", (1, 16, 18))
("ABAB", (0, 3, 7, 10))
```

같은 Li configuration이라도 stacking sequence가 다르면 별도 candidate로 평가됩니다.

---

## 5. AA/AB stacking의 물리적 정의

이 코드에서 stacking은 단순 fractional coordinate shift가 아니라 graphene registry 기준으로 정의됩니다.

### 5.1 A registry

```text
A layer:
graphene hexagon center들이 reference A layer와 겹치는 registry
```

### 5.2 B registry

```text
B layer:
shift된 layer의 carbon atom, 즉 hexagon vertex가
reference A layer의 graphene hexagon center 위/아래에 놓이는 registry
```

즉, AB stacking은 다음 조건을 만족하도록 생성됩니다.

```text
reference graphene hexagon center
↔ shifted layer carbon vertex
```

---

## 6. AB shift 자동 검출

`--ab-shift`를 지정하지 않으면 B-layer in-plane shift는 자동으로 계산됩니다.

자동 검출 절차:

1. carbon atom을 fractional z-coordinate 기준으로 carbon layer로 분류
2. 첫 carbon layer를 reference A layer로 사용
3. Li site xy projection을 graphene hexagon center 후보로 사용
4. reference carbon layer의 carbon atom을 carbon vertex 후보로 사용
5. carbon vertex가 hexagon center 위치로 이동하도록 하는 가장 짧은 in-plane vector 계산
6. 이 vector를 B-layer shift로 사용

summary JSON에는 다음과 유사한 정보가 기록됩니다.

```json
"stacking_definition": {
  "mode": "auto_ring_center_to_carbon_vertex",
  "definition": {
    "A": "hexagon centers overlap between carbon layers",
    "B": "a shifted-layer carbon vertex is placed above/below a reference-layer hexagon center"
  },
  "b_shift_supercell_frac": [ ... ],
  "center_to_vertex_distance_angstrom": ...
}
```

---

### 6.1 AB shift 수동 지정

필요하면 수동으로 B-layer shift를 지정할 수 있습니다.

```bash
python enumerate_li_configs_Final_optimized_motif.py \
  --structure LiC6.cif \
  --stacking-sequences AB \
  --ab-shift 0.1666667 0.3333333
```

주의:

- `--ab-shift`는 생성된 supercell fractional coordinate 기준입니다.
- 입력 CIF의 cell setting에 따라 올바른 값이 달라질 수 있습니다.
- 일반적으로 자동 검출을 권장합니다.

---

## 7. SOC-dependent stacking pruning

### 7.1 목적

논문 결과에 따르면 high SOC에서는 AB stacking이 energetically unfavorable/inaccessible해지고, lithiation pathway가 AA-driven route로 수렴합니다. 이를 코드에 반영하여 high SOC에서 non-all-A stacking sequence를 제외할 수 있습니다.

즉, 사용자가 `--stacking-sequences ALL`을 지정하더라도, high SOC에서는 다음처럼 줄어듭니다.

```text
4 carbon layers:
  원래: AAAA, AAAB, AABA, AABB, ABAA, ABAB, ABBA, ABBB
  high SOC pruning 후: AAAA only
```

2 carbon layers에서는:

```text
AA only
```

---

### 7.2 옵션

```bash
--soc-stacking-policy none
--soc-stacking-policy conservative
--soc-stacking-policy literature
```

기본값:

```bash
--soc-stacking-policy conservative
```

---

### 7.3 policy별 의미

#### `none`

SOC와 관계없이 사용자가 지정한 stacking sequence를 모두 사용합니다.

```bash
--soc-stacking-policy none
```

사용 상황:

- 모든 stacking possibility를 그대로 열어두고 싶을 때
- 논문 기반 pruning 없이 완전한 탐색 공간을 유지하고 싶을 때

---

#### `conservative`

기본값입니다.

```text
SOC ≥ 0.50:
  all-A stacking only
```

예:

```text
2 layers: AA only
4 layers: AAAA only
```

사용 상황:

- SOC 50% 이상에서 candidate 수를 줄이고 싶지만, SOC 40-50% 영역은 보수적으로 열어두고 싶을 때

---

#### `literature`

논문 결과를 더 적극적으로 반영합니다.

```text
SOC ≥ 0.40:
  all-A stacking only
```

사용 상황:

- 논문에서 제시한 high-SOC AA-dominant regime을 강하게 적용하고 싶을 때
- AB stacking이 40% 이상에서 물리적으로 의미 없는 후보를 많이 만든다고 판단할 때

---

### 7.4 threshold 직접 지정

```bash
--aa-only-soc-threshold 0.50
```

예:

```bash
python enumerate_li_configs_Final_optimized_motif.py \
  --structure LiC6.cif \
  --ga \
  --stacking-sequences ALL \
  --soc-stacking-policy conservative \
  --aa-only-soc-threshold 0.55
```

이 경우 SOC 55% 이상에서 all-A only가 적용됩니다.

---

### 7.5 로그 예

```text
SOC=0.6250 stacking policy: SOC >= 0.500: AA-only stacking constraint active | active_sequences=['AAAA']
GA mode: raw_candidates is the search-space size only; full enumeration is not performed.
```

해석:

- SOC 0.625에서 stacking pruning이 작동
- active stacking sequence는 `AAAA` 하나
- raw candidate 수는 전체 탐색 공간 크기일 뿐이며, GA에서는 전부 계산하지 않음

---

## 8. GA에서 hollow Li site 확장

### 8.1 목적

초기 입력 구조에 정의된 Li site만 쓰면, 같은 Li layer의 비어 있는 graphene hollow site가 탐색에서 제외될 수 있습니다. 최종 버전에서는 GA 사용 시 다음 site들을 Li candidate site로 사용합니다.

```text
입력 구조의 Li site
+ graphene hexagon center hollow site
```

---

### 8.2 기본 동작

GA에서는 hollow site 확장이 기본적으로 켜져 있습니다.

```bash
python enumerate_li_configs_Final_optimized_motif.py \
  --structure LiC6.cif \
  --ga
```

이는 다음과 동일합니다.

```bash
python enumerate_li_configs_Final_optimized_motif.py \
  --structure LiC6.cif \
  --ga \
  --ga-expand-li-sites
```

---

### 8.3 hollow site 확장 끄기

입력 구조에 이미 정의된 Li site만 사용하려면:

```bash
python enumerate_li_configs_Final_optimized_motif.py \
  --structure LiC6.cif \
  --ga \
  --no-ga-expand-li-sites
```

---

### 8.4 hollow site 확장 단계의 중복 제거

hollow site 확장 단계에서는 좌표 기준으로 중복 Li candidate site를 제거합니다.

즉, 다음 경우는 하나의 Li candidate site로 통합됩니다.

```text
- 입력 Li site와 새로 검출된 hollow site가 같은 위치
- 서로 다른 ring 검출 경로에서 같은 hollow center가 반복 생성
- periodic boundary 조건에서 같은 xy 위치로 겹치는 hollow site
```

summary JSON에는 다음과 같은 항목이 기록됩니다.

```json
"ga_li_candidate_site_generation": {
  "enabled": true,
  "n_input_li_sites": 16,
  "n_ga_li_candidate_sites": 48,
  "n_hollow_xy_before_dedupe": 64,
  "n_hollow_xy_duplicates_removed": 16
}
```

---

### 8.5 hollow site 관련 옵션

```bash
--cc-bond-cutoff 1.85
```

C-C bond 판단 cutoff입니다. graphene ring/hollow site 검출에 사용됩니다.

```bash
--hollow-site-tol 0.35
```

hollow site 좌표 중복 제거 tolerance입니다.

---

## 9. 같은 Li layer 내 Li-Li distance constraint

GA에서는 같은 Li layer 안에서 Li-Li 거리가 너무 가까운 configuration을 금지합니다.

기본값:

```bash
--ga-min-li-li-same-layer 4.35
```

의미:

```text
같은 Li layer 내 모든 Li-Li distance ≥ 4.35 Å
```

적용되는 단계:

- initial population 생성
- vacancy-swap mutation
- layer-swap mutation
- facing mutation
- fallback valid occupation pool 생성

constraint 완화 예:

```bash
python enumerate_li_configs_Final_optimized_motif.py \
  --structure LiC6.cif \
  --ga \
  --ga-min-li-li-same-layer 4.0
```

---

## 10. high-SOC motif-biased GA

### 10.1 목적

high SOC에서는 전체 random 탐색보다 물리적으로 유리한 motif를 우선 sampling하는 것이 효율적입니다.

논문 기반 해석을 반영하여 SOC 50% 이상에서는 다음 경향을 GA initial population과 mutation에 반영합니다.

1. AA stacking only
2. RH-like / Alternating layer filling 우선
3. same-layer vacancy가 균일하게 분포하도록 유도 또는 강제
4. interlayer facing pair를 증가시키는 mutation 강화
5. random generation 대신 motif-seeded initial population 사용

기본 high-SOC 초기 population 비율:

```text
SOC ≥ 0.50:
  70% = RH-like / Alternating, layer-balanced seed
  20% = facing-pair enhanced seed
  10% = random seed
```

---

### 10.2 motif seeding 활성 옵션

```bash
--ga-motif-seeding none
--ga-motif-seeding high-soc
--ga-motif-seeding always
```

기본값:

```bash
--ga-motif-seeding high-soc
```

#### `none`

motif seeding을 사용하지 않습니다. 기존 random initial population을 사용합니다.

```bash
--ga-motif-seeding none
```

#### `high-soc`

SOC가 threshold 이상일 때만 motif seeding을 사용합니다.

```bash
--ga-motif-seeding high-soc
```

기본 threshold:

```bash
--motif-seed-soc-threshold 0.50
```

#### `always`

모든 SOC에서 motif seeding을 사용합니다.

```bash
--ga-motif-seeding always
```

일반적으로 low SOC에서는 구조 다양성이 중요하므로 `always`는 신중하게 사용하십시오.

---

### 10.3 motif seeding threshold

```bash
--motif-seed-soc-threshold 0.50
```

의미:

```text
SOC ≥ 0.50에서 motif-seeded initial population 사용
```

예:

```bash
python enumerate_li_configs_Final_optimized_motif.py \
  --structure LiC6.cif \
  --ga \
  --motif-seed-soc-threshold 0.60
```

이 경우 SOC 60% 이상에서 motif seeding이 작동합니다.

---

### 10.4 motif seed 비율

기본값:

```bash
--motif-rh-fraction 0.70
--motif-facing-fraction 0.20
--motif-random-fraction 0.10
```

의미:

```text
initial population:
  70% RH-like / Alternating seed
  20% facing-pair enhanced seed
  10% random seed
```

비율은 반드시 합이 정확히 1.0일 필요는 없습니다. 내부적으로 normalize되거나 target population size에 맞춰 정수 개수로 변환됩니다.

---

### 10.5 RH-like / Alternating seed

RH-like / Alternating seed는 Li가 특정 layer나 특정 region에 과도하게 몰리지 않도록 layer-balanced하게 배치하는 seed입니다.

목표:

```text
- layer별 Li/vacancy 분포 균일화
- high SOC에서 AA-driven Alternating motif를 우선 탐색
- 특정 Li layer에 vacancy가 몰리는 후보 억제
```

이는 hard constraint가 아니라, 기본적으로 초기 population bias입니다.

---

### 10.6 facing-pair enhanced seed

facing-pair enhanced seed는 서로 다른 Li layer 사이에서 xy 위치가 같거나 가까운 Li site pair, 즉 interlayer facing pair를 많이 포함하도록 만든 seed입니다.

목표:

```text
- high SOC에서 Li-Li interlayer interaction을 반영
- SOC가 높아질수록 facing tendency가 중요해지는 경향을 sampling에 반영
```

---

### 10.7 random seed

motif bias가 너무 강해져 탐색 다양성이 줄어드는 것을 방지하기 위해 일부 population은 random으로 유지합니다.

기본값:

```text
10% random
```

---

### 10.8 로그 예

```text
SOC=0.6250 GA motif seeding active:
RH-like/alternating=0.70, facing-enhanced=0.20, random=0.10;
facing_mutation_prob=0.20, uniform_vacancy_mode=seed

GA motif-seeded initial population:
RH-like/alternating=84, facing-enhanced=24, random=12
```

해석:

- SOC 0.625에서 motif seeding 활성
- population 120 기준 RH-like seed 84개, facing seed 24개, random 12개 생성
- 이후 GA mutation에는 facing mutation이 포함됨

---

## 11. same-layer vacancy uniformity

### 11.1 목적

high SOC에서 RH-like / Alternating filling을 우선하려면 vacancy가 특정 Li layer에 몰리는 구조보다, layer별 vacancy가 균일하게 분포한 구조를 우선하는 것이 유리합니다.

이를 위해 다음 옵션을 사용합니다.

```bash
--ga-uniform-vacancy-mode none
--ga-uniform-vacancy-mode seed
--ga-uniform-vacancy-mode hard
```

기본값:

```bash
--ga-uniform-vacancy-mode seed
```

---

### 11.2 vacancy count 정의

각 Li layer마다 vacancy 수를 계산합니다.

```text
vacancy_count(layer) = 해당 Li layer의 candidate site 수 - 해당 layer에 점유된 Li 수
```

균일성 조건:

```text
max(vacancy_count) - min(vacancy_count) ≤ --uniform-vacancy-tolerance
```

기본 tolerance:

```bash
--uniform-vacancy-tolerance 1
```

---

### 11.3 `none`

```bash
--ga-uniform-vacancy-mode none
```

vacancy uniformity를 고려하지 않습니다.

사용 상황:

- motif bias 없이 구조 다양성을 최대한 유지하고 싶을 때
- high SOC에서 vacancy distribution을 제한하고 싶지 않을 때

---

### 11.4 `seed`

```bash
--ga-uniform-vacancy-mode seed
```

기본값입니다.

동작:

- initial motif seed 생성 시 vacancy가 균일하게 분포하도록 유도
- random generation 및 mutation 결과에는 hard filter를 걸지 않음

장점:

- high SOC search space를 물리적으로 유리한 영역으로 bias
- 그래도 GA가 다른 구조를 탐색할 수 있음
- 일반적으로 가장 안전한 기본값

---

### 11.5 `hard`

```bash
--ga-uniform-vacancy-mode hard
```

동작:

- initial population
- random candidate
- vacancy swap
- layer swap
- facing mutation

모든 단계에서 vacancy uniformity 조건을 만족하지 않는 candidate는 invalid로 버립니다.

예:

```text
Li layer별 vacancy 수 = [2, 2, 1, 1]
max-min = 1 → tolerance 1에서 통과

Li layer별 vacancy 수 = [4, 0, 1, 0]
max-min = 4 → tolerance 1에서 탈락
```

장점:

- candidate 수를 크게 줄일 수 있음
- RH-like / layer-balanced search를 강하게 적용

주의:

- 너무 강한 제약일 수 있음
- high SOC에서 valid candidate 수가 매우 적어질 수 있음
- GA initial population 생성이 어려워지거나 특정 SOC가 skip될 수 있음

완화 예:

```bash
--ga-uniform-vacancy-mode hard --uniform-vacancy-tolerance 2
```

---

## 12. facing-pair mutation

### 12.1 목적

high SOC에서는 interlayer Li-Li interaction이 중요해지므로, 서로 다른 Li layer에서 xy 위치가 같은 Li site pair, 즉 facing pair를 증가시키는 mutation을 강화합니다.

옵션:

```bash
--facing-mutation-prob 0.20
```

기본값:

```text
0.20
```

의미:

```text
새 generation을 만들 때 일부 mutation을 facing-pair enhanced mutation으로 배정
```

값을 높이면 facing-rich candidate가 더 많이 생성됩니다.

예:

```bash
python enumerate_li_configs_Final_optimized_motif.py \
  --structure LiC6.cif \
  --ga \
  --facing-mutation-prob 0.35
```

주의:

- 너무 크게 설정하면 diversity가 줄어들 수 있습니다.
- RH-like / Alternating seed와 random seed를 함께 유지하는 것이 좋습니다.

---

## 13. high SOC motif 옵션 권장값

### 13.1 기본 권장

```bash
--ga-motif-seeding high-soc
--motif-seed-soc-threshold 0.50
--motif-rh-fraction 0.70
--motif-facing-fraction 0.20
--motif-random-fraction 0.10
--ga-uniform-vacancy-mode seed
--facing-mutation-prob 0.20
```

일반적으로 이 조합이 가장 안전합니다.

---

### 13.2 더 강한 pruning

```bash
--soc-stacking-policy literature
--ga-uniform-vacancy-mode hard
--uniform-vacancy-tolerance 1
```

효과:

- SOC 40% 이상에서 all-A only
- vacancy uniformity를 hard constraint로 적용
- candidate 수를 크게 줄임

주의:

- 일부 SOC에서 valid candidate가 부족할 수 있습니다.

---

### 13.3 diversity 보존

```bash
--ga-uniform-vacancy-mode seed
--motif-random-fraction 0.20
--facing-mutation-prob 0.10
```

효과:

- motif bias는 유지하되 random exploration을 더 확보

---

## 14. GA 알고리즘

현재 GA는 crossover 없는 mutation-selection evolutionary search입니다.

### 14.1 Candidate representation

```python
candidate = (stacking_sequence, occupied_Li_sites)
```

예:

```python
("AB", (1, 2, 5, 19, 21, 22))
```

---

### 14.2 Energy cache

같은 candidate는 재평가하지 않습니다.

```python
energy_cache[(stacking_sequence, occupied_Li_sites)]
```

같은 Li configuration이라도 stacking sequence가 다르면 다른 candidate입니다.

```python
energy_cache[("AA", occupied)]
energy_cache[("AB", occupied)]
```

---

### 14.3 Survivor 기준

각 generation에서 survivor는 generation-local minimum을 기준으로 결정됩니다.

```text
E(candidate) ≤ E_min_generation + ga_threshold_ev
```

기본값:

```bash
--ga-threshold-ev 0.5
```

예:

```text
E_min_generation = -20.0 eV
ga_threshold_ev = 0.5 eV

survivor 조건:
E(candidate) ≤ -19.5 eV
```

Survivor는 다음 generation에 유지되고 mutation parent pool로 사용됩니다.

---

### 14.4 Mutation 종류

#### vacancy swap

Li 하나를 제거하고 다른 valid vacant site 하나를 채웁니다.

```text
("AB", (1, 2, 5, 19, 21, 22))
→
("AB", (1, 2, 6, 19, 21, 22))
```

#### layer swap

Li 하나를 다른 Li layer의 valid vacant site로 이동합니다.

#### stacking mutation

Li configuration은 유지하고 stacking sequence만 변경합니다.

```text
("AA", occupied) → ("AB", occupied)
```

4층 예:

```text
("AAAA", occupied) → ("ABAB", occupied)
```

stacking mutation 비율:

```bash
--stacking-mutation-prob 0.15
```

#### facing mutation

high SOC에서 interlayer facing pair를 증가시키는 방향으로 Li를 재배치합니다.

옵션:

```bash
--facing-mutation-prob 0.20
```

---

## 15. GA 중복 제거 전략

최종 버전에서는 GA 효율과 최종 CIF의 구조적 독립성을 위해 여러 단계의 중복 제거를 사용합니다.

### 15.1 1단계: hollow site 확장 단계

좌표 기준으로 중복 Li candidate site를 제거합니다.

summary 항목:

```json
"n_hollow_xy_before_dedupe": 64,
"n_hollow_xy_duplicates_removed": 16
```

---

### 15.2 2단계: GA candidate 생성 단계

GA population 생성 및 mutation 과정에서 다음 중복을 제거합니다.

#### candidate key 중복

동일한 candidate key:

```python
(stacking_sequence, occupied_Li_sites)
```

가 같은 generation 또는 candidate 생성 과정에서 반복되면 제거합니다.

#### coordinate hash 중복

candidate key가 달라도 실제 생성되는 좌표가 같으면 빠른 coordinate hash로 중복을 제거합니다.

summary 항목:

```json
"n_ga_candidate_key_duplicates_removed": 12,
"n_ga_candidate_coordinate_duplicates_removed": 5
```

---

### 15.3 3단계: GA energy cache

이미 평가된 candidate는 다시 Ewald energy를 계산하지 않습니다.

```python
if candidate in energy_cache:
    return energy_cache[candidate]
```

summary 항목:

```json
"n_evaluated_candidates": 360
```

---

### 15.4 4단계: GA 완료 후 energy 정렬

GA가 끝나면 실제 평가된 모든 candidate를 Ewald energy 오름차순으로 정렬합니다.

```text
lowest energy → highest energy
```

---

### 15.5 5단계: 최종 CIF 선택 단계의 symmetry-equivalence deduplication

최종 CIF 저장 전에는 `pymatgen.analysis.structure_matcher.StructureMatcher`를 사용하여 symmetry-equivalent 중복을 제거합니다.

중복으로 보는 경우:

```text
- periodic translation
- origin shift
- rotation/reflection
- site permutation
- fixed supercell 내 symmetry-equivalent structure
```

처리 방식:

```text
energy 오름차순 후보 list
  ↓
exact coordinate hash 중복 먼저 제거
  ↓
StructureMatcher로 symmetry-equivalence 검사
  ↓
unique structure이면 top list에 추가
  ↓
symmetry-unique top 20이 모이면 즉시 중단
```

summary 항목:

```json
"n_exact_duplicates_removed": 10,
"n_symmetry_duplicates_skipped_for_top20": 35,
"n_symmetry_unique_cifs_written": 20
```

---

## 16. 중복 제거 옵션

기본값은 symmetry-equivalence 기준입니다.

```bash
--dedup-mode symmetry
```

가능한 값:

```bash
--dedup-mode symmetry
--dedup-mode exact
--dedup-mode none
```

### 16.1 `symmetry`

StructureMatcher를 사용해 symmetry-equivalent 구조를 제거합니다. 기본값입니다.

### 16.2 `exact`

좌표 hash가 같은 구조만 제거합니다. 더 빠르지만 translation/rotation/reflection equivalent 구조는 남을 수 있습니다.

### 16.3 `none`

최종 CIF 선택 단계의 중복 제거를 비활성화합니다. 일반적으로 권장하지 않습니다.

---

### 16.4 StructureMatcher tolerance 옵션

```bash
--dedup-ltol 0.2
--dedup-stol 0.3
--dedup-angle-tol 5.0
```

- `--dedup-ltol`: lattice length tolerance
- `--dedup-stol`: site tolerance
- `--dedup-angle-tol`: lattice angle tolerance

---

## 17. Full enumeration에서 stacking sequence 적용

`--ga`가 없어도 `--stacking-sequences`는 적용됩니다.

예:

```bash
python enumerate_li_configs_Final_optimized_motif.py \
  --structure LiC6.cif \
  --stacking-sequences AA AB
```

동작:

```text
각 SOC:
  각 stacking sequence:
    Li/vacancy configuration enumeration
    Ewald energy 계산
  모든 stacking 후보를 energy 기준 통합 정렬
  최종 CIF 선택 단계에서 dedup
  symmetry-unique top 20 CIF 저장
```

symmetry-aware mode에서는 stacking sequence마다 symmetry operation을 따로 계산합니다. AB-type stacking은 AA보다 symmetry가 낮을 수 있기 때문입니다.

---

## 18. Energy 계산

각 candidate에 대해 Ewald electrostatic energy를 계산합니다.

산화수 설정:

```text
Li = +1
C  = -n_Li / n_C
```

`n_Li=0`이면 energy는 `0.0 eV`로 처리합니다.

같은 SOC에서는 C-C contribution이 일정하므로, energy 차이는 주로 다음을 반영합니다.

```text
Li arrangement
carbon stacking registry
Li/hollow site position
```

최종 optimized 버전에서는 stacking sequence별 Ewald matrix를 precompute하여, 각 configuration의 Ewald energy 평가를 빠르게 수행합니다.

---

## 19. 출력 파일 구조

기본 output directory:

```bash
./li_configs
```

변경:

```bash
--output my_output_dir
```

예상 구조:

```text
li_configs/
  SOC_0.7500_nLi12/
    top_stable_cif/
      rank001_AA_nLi12_SOC0.7500_E-25.458eV.cif
      rank002_AA_nLi12_SOC0.7500_E-25.100eV.cif
      ...
    ga_convergence_nLi12.png
  enumeration_summary.json
```

---

## 20. CIF 저장 규칙

각 SOC별로 최대 20개 구조만 CIF로 저장합니다.

```text
symmetry-unique 후보 수 < 20:
  전부 CIF 저장

symmetry-unique 후보 수 ≥ 20:
  Ewald energy 기준 top 20 CIF 저장
```

GA에서는 top 20이 GA survivor만에서 선택되는 것이 아니라, 실제 평가된 전체 candidate에서 energy 기준으로 선택됩니다.

---

## 21. 로그 해석

예:

```text
SOC 11/17 START | n_Li=10/16, SOC=0.6250, C(48,k)=6,540,715,896, stacking_sequences=1, raw_candidates=6,540,715,896
```

의미:

- 입력 기준 SOC는 `10/16 = 0.625`
- hollow 확장 후 GA Li candidate site가 48개
- Li 조합 수는 `C(48,10)=6,540,715,896`
- SOC-dependent stacking pruning으로 active stacking sequence가 1개
- raw candidate 수는 전체 탐색 공간 크기
- GA에서는 이 수를 전부 계산하지 않음

GA 완료 로그:

```text
GA complete: 720 unique candidates evaluated (95 survivors), global best E = -35.1234 eV
```

의미:

- raw candidate 전체 중 GA가 720개를 실제 평가
- survivor criterion을 통과한 candidate는 95개
- 최저 Ewald energy는 `-35.1234 eV`

최종 selection 로그:

```text
SOC=0.6250: final selection from 720 GA-evaluated candidates (mode=symmetry):
exact duplicates removed=12,
symmetry duplicates skipped while selecting top 20=35,
symmetry-unique CIF candidates=20.
```

의미:

- GA 평가 후보 720개를 energy 오름차순으로 정렬
- exact duplicate 12개 제거
- StructureMatcher 기준 symmetry-equivalent 중복 35개 skip
- 최종적으로 symmetry-unique CIF 후보 20개 확보

---

## 22. `enumeration_summary.json`

계산 종료 후 output directory에 저장됩니다.

```text
li_configs/enumeration_summary.json
```

주요 상위 항목:

```json
{
  "structure_file": "...",
  "supercell": [2, 2, 2],
  "space_group": "...",
  "ga_enabled": true,
  "n_carbon_layers": 2,
  "stacking_sequences": ["AA", "AB"],
  "soc_stacking_policy": "conservative",
  "aa_only_soc_threshold": 0.50,
  "ga_motif_seeding": "high-soc",
  "motif_seed_soc_threshold": 0.50,
  "ga_motif_seed_fractions": {
    "rh_like": 0.70,
    "facing": 0.20,
    "random": 0.10
  },
  "ga_uniform_vacancy_mode": "seed",
  "facing_mutation_prob": 0.20,
  "max_cif_per_soc": 20,
  "duplicate_structure_filter": "symmetry",
  "soc_levels": [ ... ],
  "stable_by_soc": [ ... ]
}
```

---

### 22.1 SOC별 summary 항목

예:

```json
{
  "n_Li": 10,
  "SOC": 0.625,
  "n_raw_combinations_per_stacking": 6540715896,
  "n_stacking_sequences": 1,
  "n_raw_stacking_candidates": 6540715896,
  "active_stacking_sequences": ["AA"],
  "stacking_policy_applied": true,
  "motif_seeding_active": true,
  "method": "GA",
  "best_ewald_energy_eV": -35.123456,
  "best_stacking_sequence": "AA",
  "best_occupied_Li_sites": [0, 4, 8, 12, 16, 20, 24, 28, 32, 36],
  "n_ga_survivors": 95,
  "n_evaluated_candidates": 720,
  "n_exact_duplicates_removed": 12,
  "n_symmetry_duplicates_skipped_for_top20": 35,
  "n_symmetry_unique_cifs_written": 20,
  "n_ga_candidate_key_duplicates_removed": 8,
  "n_ga_candidate_coordinate_duplicates_removed": 4,
  "n_motif_seed_rh_like": 84,
  "n_motif_seed_facing": 24,
  "n_motif_seed_random": 12,
  "top_stable": [...]
}
```

skip된 SOC 예:

```json
{
  "n_Li": 15,
  "SOC": 0.9375,
  "n_unique_or_found": 0,
  "method": "GA",
  "best_ewald_energy_eV": NaN,
  "best_stacking_sequence": null,
  "best_occupied_Li_sites": [],
  "n_ga_survivors": 0,
  "n_evaluated_candidates": 0,
  "top_stable": []
}
```

---

### 22.2 `stable_by_soc`

각 SOC별 최안정 구조만 모은 요약입니다.

```json
"stable_by_soc": [
  {
    "n_Li": 10,
    "SOC": 0.625,
    "method": "GA",
    "best_stacking_sequence": "AA",
    "best_ewald_energy_eV": -35.123456,
    "best_occupied_Li_sites": [0, 4, 8, 12, 16, 20, 24, 28, 32, 36]
  }
]
```

---

## 23. 주요 옵션 정리

### 입력/출력

| 옵션 | 기본값 | 설명 |
|---|---:|---|
| `--structure`, `-s` | 필수 | 입력 CIF 또는 POSCAR |
| `--output`, `-o` | `./li_configs` | 출력 directory |
| `--supercell NX NY NZ` | `2 2 2` | supercell expansion |
| `--soc-levels N ...` | all | 계산할 점유 Li 개수 |
| `--no-poscar` | false | POSCAR 출력 skip |
| `--sig-figs` | `8` | POSCAR coordinate significant figures |

### Enumeration

| 옵션 | 기본값 | 설명 |
|---|---:|---|
| `--no-symmetry` | false | symmetry reduction 없이 모든 조합 평가 |
| `--symprec` | `0.01` | SpacegroupAnalyzer tolerance |
| `--match-tol` | `0.01` | Li-site matching tolerance |

### GA

| 옵션 | 기본값 | 설명 |
|---|---:|---|
| `--ga` | false | GA 활성화 |
| `--ga-pop` | `120` | GA population size |
| `--ga-gens` | `60` | GA generation 수 |
| `--ga-threshold-ev` | `0.5` | survivor threshold |
| `--ga-seed` | `42` | random seed |
| `--stacking-mutation-prob` | `0.15` | stacking mutation 비율 |
| `--facing-mutation-prob` | `0.20` | facing-pair mutation 비율 |

### SOC stacking policy

| 옵션 | 기본값 | 설명 |
|---|---:|---|
| `--soc-stacking-policy` | `conservative` | `none`, `conservative`, `literature` |
| `--aa-only-soc-threshold` | `0.50` | all-A only 적용 SOC threshold |

### Motif-biased GA

| 옵션 | 기본값 | 설명 |
|---|---:|---|
| `--ga-motif-seeding` | `high-soc` | `none`, `high-soc`, `always` |
| `--motif-seed-soc-threshold` | `0.50` | motif seeding 시작 SOC |
| `--motif-rh-fraction` | `0.70` | RH-like / Alternating seed 비율 |
| `--motif-facing-fraction` | `0.20` | facing-enhanced seed 비율 |
| `--motif-random-fraction` | `0.10` | random seed 비율 |
| `--ga-uniform-vacancy-mode` | `seed` | `none`, `seed`, `hard` |
| `--uniform-vacancy-tolerance` | `1` | layer별 vacancy 수 차이 허용값 |

### Stacking

| 옵션 | 기본값 | 설명 |
|---|---:|---|
| `--stacking-sequences` | all-A | 직접 sequence 지정 또는 `ALL` |
| `--stacking-modes` | None | legacy shortcut: `AA`, `AB` |
| `--ab-shift DA DB` | auto | B-layer shift 수동 지정 |
| `--carbon-layer-tol` | `0.05` | carbon layer grouping tolerance |

### GA Li site 확장 및 Li-Li constraint

| 옵션 | 기본값 | 설명 |
|---|---:|---|
| `--ga-expand-li-sites` | true | GA에서 hollow site 확장 |
| `--no-ga-expand-li-sites` | false | 입력 Li site만 사용 |
| `--ga-min-li-li-same-layer` | `4.35` Å | 같은 Li layer 내 최소 Li-Li 거리 |
| `--hollow-site-tol` | `0.35` Å | hollow site 중복 제거 tolerance |
| `--cc-bond-cutoff` | `1.85` Å | C-C bond cutoff |

### 중복 제거

| 옵션 | 기본값 | 설명 |
|---|---:|---|
| `--dedup-mode` | `symmetry` | `symmetry`, `exact`, `none` |
| `--dedup-ltol` | `0.2` | StructureMatcher lattice length tolerance |
| `--dedup-stol` | `0.3` | StructureMatcher site tolerance |
| `--dedup-angle-tol` | `5.0` | StructureMatcher angle tolerance |

---

## 24. 추천 실행 예시

### 24.1 기본 추천: GA + stacking ALL + high-SOC motif bias

```bash
python enumerate_li_configs_Final_optimized_motif.py \
  --structure LiC6.cif \
  --ga \
  --stacking-sequences ALL
```

기본적으로 다음이 적용됩니다.

```text
SOC ≥ 0.50:
  all-A stacking only
  motif-seeded initial population
  RH-like 70%, facing 20%, random 10%
  uniform vacancy mode = seed
```

---

### 24.2 논문 기준을 더 강하게 반영

```bash
python enumerate_li_configs_Final_optimized_motif.py \
  --structure LiC6.cif \
  --ga \
  --stacking-sequences ALL \
  --soc-stacking-policy literature
```

효과:

```text
SOC ≥ 0.40:
  all-A stacking only
```

---

### 24.3 high SOC search space를 더 강하게 줄이기

```bash
python enumerate_li_configs_Final_optimized_motif.py \
  --structure LiC6.cif \
  --ga \
  --stacking-sequences ALL \
  --soc-stacking-policy literature \
  --ga-uniform-vacancy-mode hard \
  --uniform-vacancy-tolerance 1
```

---

### 24.4 high SOC에서 valid candidate가 부족할 때 완화

```bash
python enumerate_li_configs_Final_optimized_motif.py \
  --structure LiC6.cif \
  --ga \
  --stacking-sequences ALL \
  --ga-uniform-vacancy-mode seed \
  --uniform-vacancy-tolerance 2 \
  --ga-min-li-li-same-layer 4.0
```

---

### 24.5 motif seeding 끄기

```bash
python enumerate_li_configs_Final_optimized_motif.py \
  --structure LiC6.cif \
  --ga \
  --stacking-sequences ALL \
  --ga-motif-seeding none
```

---

### 24.6 기존처럼 모든 SOC에서 모든 stacking 유지

```bash
python enumerate_li_configs_Final_optimized_motif.py \
  --structure LiC6.cif \
  --ga \
  --stacking-sequences ALL \
  --soc-stacking-policy none
```

---

### 24.7 GA 없이 stacking-aware full enumeration

```bash
python enumerate_li_configs_Final_optimized_motif.py \
  --structure LiC6.cif \
  --stacking-sequences AA AB
```

---

### 24.8 입력 Li site만 사용

```bash
python enumerate_li_configs_Final_optimized_motif.py \
  --structure LiC6.cif \
  --ga \
  --stacking-sequences ALL \
  --no-ga-expand-li-sites
```

---

## 25. 문제 해결

### 25.1 `BooleanOptionalAction` 오류

오류:

```text
AttributeError: module 'argparse' has no attribute 'BooleanOptionalAction'
```

해결:

- 이 최종 버전은 Python 3.8 호환 방식입니다.
- `--ga-expand-li-sites`와 `--no-ga-expand-li-sites`는 별도 flag로 구현되어 있습니다.

---

### 25.2 high SOC에서 valid Li configuration 없음

로그:

```text
skipped GA search because no valid Li configuration satisfies d_Li-Li(same layer) >= 4.350 Å
```

해결:

```bash
--ga-min-li-li-same-layer 4.0
```

또는:

```bash
--ga-uniform-vacancy-mode seed
```

또는 hard mode 사용 중이면:

```bash
--uniform-vacancy-tolerance 2
```

---

### 25.3 symmetry dedup이 너무 느림

최종 top 20 selection에서 StructureMatcher를 사용하므로 후보가 많고 구조가 복잡하면 시간이 걸릴 수 있습니다.

해결 옵션:

```bash
--dedup-mode exact
```

또는:

```bash
--dedup-mode none
```

단, `none`은 중복 CIF가 저장될 수 있으므로 권장하지 않습니다.

---

### 25.4 AB stacking이 의도와 다름

자동 AB shift가 입력 CIF convention과 맞지 않을 수 있습니다. 이 경우 수동 shift를 지정하십시오.

```bash
--ab-shift DA DB
```

`DA`, `DB`는 supercell fractional coordinate 기준입니다.

---

### 25.5 raw_candidates가 너무 큼

예:

```text
raw_candidates=52,325,727,168
```

해석:

- 이것은 전체 search-space 크기입니다.
- GA에서는 전체를 계산하지 않습니다.
- 실제 Ewald 계산 수는 다음 로그에서 확인합니다.

```text
GA complete: XXX unique candidates evaluated
```

candidate 수를 줄이려면:

```bash
--soc-stacking-policy literature
--ga-uniform-vacancy-mode hard
--uniform-vacancy-tolerance 1
```

또는 motif seeding은 유지하되 diversity를 줄입니다.

```bash
--motif-random-fraction 0.05
--facing-mutation-prob 0.30
```

---

## 26. 전체 workflow 요약

```text
입력 structure 읽기
  ↓
supercell 생성
  ↓
Li / C sublattice 분리
  ↓
carbon layer, Li layer 검출
  ↓
stacking sequence 생성
  ↓
SOC-dependent stacking policy 준비
  ↓
AB shift 자동 검출 또는 수동 적용
  ↓
stacking sequence별 carbon coordinate 생성
  ↓
GA 사용 시:
    hollow Li site 확장
    hollow site 좌표 중복 제거
    same-layer Li-Li compatibility matrix 생성
  ↓
각 SOC loop
    ↓
    SOC-dependent active stacking sequence 결정
    ↓
    high-SOC motif seeding 활성 여부 결정
    ↓
    GA 또는 full enumeration
    ↓
    candidate key / coordinate hash 중복 억제
    ↓
    Ewald energy 계산 또는 precomputed Ewald matrix 평가
    ↓
    평가된 후보 전체를 energy 오름차순 정렬
    ↓
    최종 CIF 선택:
      exact coordinate hash dedup
      StructureMatcher symmetry-equivalence dedup
      symmetry-unique top 20 확보 시 중단
    ↓
    SOC별 최안정 stacking / energy / Li sites 출력
  ↓
SOC별 안정 구조 요약 출력
  ↓
enumeration_summary.json 저장
```

---

## 27. 해석상 주의

- `raw_candidates`는 전체 탐색 공간 크기입니다.
- `--ga`에서는 `raw_candidates` 전체를 계산하지 않습니다.
- `n_evaluated_candidates`는 GA가 실제 Ewald energy를 계산한 candidate 수입니다.
- `n_symmetry_unique_cifs_written`는 최종 CIF로 저장된 symmetry-unique 구조 수입니다.
- `n_ga_survivors`는 GA selection criterion을 통과한 survivor 수입니다.
- survivor 수와 CIF 저장 수는 다릅니다.
- CIF top 20은 survivor subset이 아니라, 평가된 전체 candidate에서 energy 기준으로 고른 symmetry-unique top 20입니다.
- `occupied_Li_sites`는 hollow site 확장 후의 GA Li candidate site index일 수 있습니다.
- `--ga-uniform-vacancy-mode hard`는 search space를 크게 줄이지만, 너무 강하면 valid candidate가 부족해질 수 있습니다.
- `--soc-stacking-policy literature`는 논문 기반 pruning을 적극적으로 적용하므로, 검증 목적의 완전 탐색에서는 `none` 또는 `conservative`를 사용하십시오.
