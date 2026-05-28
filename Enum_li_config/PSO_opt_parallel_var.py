import os
import sys
import uuid
import fcntl
import time
import random
import traceback
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from concurrent.futures import ThreadPoolExecutor, as_completed  # ← Thread 병렬 (외부 프로세스 호출에 적합)
from collections import defaultdict
from functools import reduce

from ase import Atoms
from ase.atom import Atom
from ase.io import read, write
from ase.build import stack
from scipy.spatial.distance import pdist
from scipy.spatial import cKDTree
from sklearn.cluster import KMeans
import networkx as nx

# ----------------------------
# 유틸: 안전한 append (파일락)
# ----------------------------
def _safe_append(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "a") as f:
        fcntl.flock(f, fcntl.LOCK_EX)
        f.write(text)
        fcntl.flock(f, fcntl.LOCK_UN)

# ----------------------------
# 분석/생성 관련 함수군 (원본과 동일/보강)
# ----------------------------
def analyze_li_layer_distances(structure_file, file_format='cif'):
    atoms = read(structure_file, format=file_format)
    li_positions = np.array([atom.position for atom in atoms if atom.symbol == "Li"])

    dists = pdist(li_positions)
    r_min = np.min(dists[dists > 0]) * 0.95 if len(dists[dists > 0]) > 0 else 1.5

    z_coords = li_positions[:, 2].reshape(-1, 1)
    unique_z = np.unique(np.round(z_coords.flatten(), 2))
    k = len(unique_z)

    if k < 2:
        z_layer_tol = 0.5  # fallback
    else:
        kmeans = KMeans(n_clusters=k, n_init=10).fit(z_coords)
        centers = sorted([c[0] for c in kmeans.cluster_centers_])
        z_layer_tol = 0.5 * np.min(np.diff(centers))

    r_min = 1.0  ## 모든 hex_center 를 사용함. (valid site = all site)
    return r_min, z_layer_tol

def estimate_layer_spacing(atoms):
    z_coords = sorted([atom.position[2] for atom in atoms if atom.symbol == 'C'])
    z_unique = np.unique(np.round(z_coords, 3))
    if len(z_unique) < 2:
        raise ValueError("Carbon layer가 2개 이상 필요합니다.")
    diffs = np.diff(z_unique)
    spacing = np.mean(diffs)
    return spacing

def generate_graphite_stack(stack_type, layers=4, spacing=3.35, Ncell=2, struc_type=100):
    a = 2.46
    a1 = np.array([a, 0])
    a2 = np.array([a/2, a * np.sqrt(3) / 2])

    if struc_type == 100:
        shift_dict = {'A': np.array([0.0, 0.0]), 'B': (1.0 / (3.0 * Ncell)) * a1, 'C': (2.0 / (3.0 * Ncell)) * a1}
    elif struc_type == 75:
        shift_dict = {'A': np.array([0.0, 0.0]), 'B':  (1.0 / (6.0 * Ncell)) * a1 - (1.0 / (6.0 * Ncell)) * a2,
                      'C': -(1.0 / (6.0 * Ncell)) * a1 + (1.0 / (6.0 * Ncell)) * a2}
    else:
        shift_dict = {'A': np.array([0.0, 0.0]), 'B': (1.0 / (3.0 * Ncell)) * a1, 'C': (2.0 / (3.0 * Ncell)) * a1}

    base_layer = Atoms('C2', positions=[(0, 0, 0), (1.23, 0, 0)], cell=[a, a, spacing * layers], pbc=True)
    stacked = base_layer.copy()
    for i, layer in enumerate(stack_type[1:], 1):
        shifted = base_layer.copy()
        shift = shift_dict.get(layer, np.array([0.0, 0.0]))
        shifted.translate([*shift, spacing * i])
        stacked += shifted
    return stacked

def generate_graphite_stack_from_input(base_c_atoms, stack_type, Ncell, struc_type):
    spacing = estimate_layer_spacing(base_c_atoms)
    cell = base_c_atoms.get_cell()
    a1_vec, a2_vec = cell[0], cell[1]

    if struc_type == 100:
        shift_dict = {'A': np.array([0.0, 0.0, 0.0]), 'B': (1.0 / (3.0 * Ncell)) * a1_vec, 'C': (2.0 / (3.0 * Ncell)) * a1_vec}
    elif struc_type == 75:
        shift_dict = {'A': np.array([0.0, 0.0, 0.0]),
                      'B':  (1.0 / (6.0 * Ncell)) * a1_vec - (1.0 / (6.0 * Ncell)) * a2_vec,
                      'C': -(1.0 / (6.0 * Ncell)) * a1_vec + (1.0 / (6.0 * Ncell)) * a2_vec}
    else:
        shift_dict = {'A': np.array([0.0, 0.0, 0.0]), 'B': (1.0 / (3.0 * Ncell)) * a1_vec, 'C': (2.0 / (3.0 * Ncell)) * a1_vec}

    z_coords = np.array([atom.position[2] for atom in base_c_atoms if atom.symbol == "C"])
    z_min = np.min(z_coords)
    layer_tol = 0.1
    base_layer_atoms = [atom for atom in base_c_atoms if atom.symbol == "C" and abs(atom.position[2] - z_min) < layer_tol]
    base_layer = Atoms(base_layer_atoms, cell=cell, pbc=True)

    stacked = base_layer.copy()
    for i, layer_label in enumerate(stack_type[1:], 1):
        shift_xy = shift_dict.get(layer_label, np.array([0.0, 0.0]))
        shift_vec = np.array([shift_xy[0], shift_xy[1], i * spacing])
        shifted = base_layer.copy()
        shifted.translate(shift_vec)
        stacked += shifted
    return stacked

def add_lithium_atoms(base_atoms, li_positions, valid_sites, snap=True):
    new_atoms = base_atoms.copy()
    tree = cKDTree(valid_sites)
    for pos in li_positions:
        if snap:
            _, idx = tree.query(pos)
            snapped_pos = valid_sites[idx]
            new_atoms.append(Atom('Li', position=snapped_pos))
        else:
            key = tuple(np.round(pos, 4))
            valid_set = set(tuple(np.round(site, 4)) for site in valid_sites)
            if key not in valid_set:
                raise ValueError(f"{key} not in valid Li sites")
            new_atoms.append(Atom('Li', position=pos))
    return new_atoms

def group_valid_sites_by_z(valid_sites, z_tol=0.2):
    z_coords = np.array([site[2] for site in valid_sites])
    bins = np.unique(np.round(z_coords / z_tol).astype(int))
    layers = [[] for _ in bins]
    for site in valid_sites:
        z_bin = int(round(site[2] / z_tol))
        idx = np.where(bins == z_bin)[0][0]
        layers[idx].append(site)
    return layers

def minimum_image_distance(p1, p2, cell):
    delta = p1 - p2
    for i in range(3):
        delta[i] -= round(delta[i] / cell[i, i]) * cell[i, i]
    return np.linalg.norm(delta)

def select_li_sites_layerwise(valid_sites, n_li, r_min, max_per_layer, z_tol=0.2, max_tries=10000, cell=None, success_counter=None):
    print("[DEBUG][select_li_sites_layerwise] called", flush=True)
    for attempt in range(max_tries):
        layers = group_valid_sites_by_z(valid_sites, z_tol)
        k = len(layers)
        valid_counts = [len(layer) for layer in layers]
        max_counts = [min(n, max_per_layer) for n in valid_counts]
        total_possible = sum(max_counts)
        if n_li > total_possible:
            continue

        for _ in range(5000):
            counts = [0] * k
            remaining = n_li
            while remaining > 0:
                candidates = [i for i in range(k) if counts[i] < max_counts[i]]
                if not candidates:
                    break
                i = np.random.choice(candidates)
                counts[i] += 1
                remaining -= 1
            if sum(counts) == n_li:
                break
        else:
            continue

        selected = []
        valid = True

        if all(c <= 1 for c in counts):
            for layer, count in zip(layers, counts):
                if count == 1:
                    layer_array = np.array(layer)
                    np.random.shuffle(layer_array)
                    selected.append(layer_array[0])
            if success_counter is not None:
                success_counter["count"] += 1
            print(f"[SUCCESS] Trivial case (all count <= 1), skipping r_min test", flush=True)
            return np.array(selected)

        for layer, count in zip(layers, counts):
            if len(layer) < count:
                valid = False
                break
            layer_array = np.array(layer)
            np.random.shuffle(layer_array)
            if count == 1:
                selected.append(layer_array[0])
                continue
            first_try_limit = 50
            found = False
            for _ in range(first_try_limit):
                np.random.shuffle(layer_array)
                selection = [layer_array[0]]
                for site in layer_array:
                    if any(np.allclose(site, s) for s in selection):
                        continue
                    is_far_enough = True
                    for s in selection:
                        d = minimum_image_distance(site, s, cell)
                        if d < r_min:
                            is_far_enough = False
                            break
                    if is_far_enough:
                        selection.append(site)
                    if len(selection) == count:
                        break
                if len(selection) == count:
                    selected.extend(selection)
                    found = True
                    break
            if not found:
                valid = False
                break

        if valid:
            if success_counter is not None:
                success_counter["count"] += 1
            print(f"[SUCCESS] Valid configuration found ", flush=True)
            return np.array(selected)
    raise RuntimeError("요구 조건을 만족하는 Li 배치를 찾지 못했습니다 (max_tries 초과)")

def get_hex_center_sites(graphite_atoms, bond_cutoff=1.6, z_tol=0.2):
    positions = np.array([atom.position for atom in graphite_atoms if atom.symbol == 'C'])
    cell = graphite_atoms.get_cell()
    print(cell, flush=True)
    a1, a2 = cell[0][:2], cell[1][:2]
    c_length = cell[2, 2]
    z_coords = np.round(positions[:, 2], 3)
    unique_layers = np.unique(z_coords)

    a1_3d = np.array([*a1, 0.0])
    a2_3d = np.array([*a2, 0.0])
    a3_3d = np.array([0.0, 0.0, c_length])
    cell_matrix = np.array([a1_3d, a2_3d, a3_3d]).T
    inv_cell = np.linalg.inv(cell_matrix)

    frac_coords = []
    estimated_spacing = None
    for i, z in enumerate(unique_layers):
        if i < len(unique_layers) - 1:
            z_next = unique_layers[i + 1]
            z_mid = 0.5 * (z + z_next)
            estimated_spacing = z_next - z
        else:
            z_mid = z + 0.5 * estimated_spacing if estimated_spacing else z + 1.7

        layer_indices = np.where(abs(z_coords - z) < z_tol)[0]
        if len(layer_indices) < 6:
            continue

        layer_pos = positions[layer_indices, :2]
        full_pos = []
        image_map = []

        for dx in range(-2, 3):
            for dy in range(-2, 3):
                shift = dx * a1 + dy * a2
                for idx, pos in enumerate(layer_pos):
                    full_pos.append(pos + shift)
                    image_map.append(idx)

        full_pos = np.array(full_pos)
        tree = cKDTree(full_pos)
        pairs = tree.query_pairs(r=bond_cutoff)

        G = nx.Graph()
        for idx in range(len(full_pos)):
            G.add_node(idx)
        for i_, j_ in pairs:
            G.add_edge(i_, j_)

        rings = [cycle for cycle in nx.cycle_basis(G) if len(cycle) == 6]
        for ring in rings:
            coords = full_pos[ring]
            cx, cy = np.mean(coords, axis=0)
            center_cart = np.array([cx, cy, z_mid])
            center_frac = np.dot(inv_cell, center_cart) % 1.0
            if abs(center_frac[2] - 1.0) < z_tol or abs(center_frac[2]) < z_tol:
                center_frac[2] = 0.0
            center_frac = np.round(center_frac, 6)
            center_frac = np.array([0.0 if np.isclose(f, 1.0) else f for f in center_frac])
            frac_coords.append(tuple(center_frac))

    unique_frac = np.unique(np.array(frac_coords), axis=0)
    unique_cart = [np.dot(cell_matrix, f) for f in unique_frac]

    print(f"\n[INFO] Number of unique hex centers: {len(unique_cart)}")
    for i, center in enumerate(unique_cart):
        x, y, z = center
        print(f"  {i+1:3d}: x = {x:.4f}, y = {y:.4f}, z = {z:.4f}")
    return np.array(unique_cart)

def lithium_index_bounds(valid_sites, n_li):
    return [(0, len(valid_sites) - 1)] * n_li

def write_mtp_cfg(atoms, file_handle, energy=None, swarmindex=None):
    file_handle.write("BEGIN_CFG\n")
    file_handle.write(f"Size\n  {len(atoms)}\n")
    file_handle.write("Supercell\n")
    cell = atoms.get_cell()
    for row in cell:
        file_handle.write("  {:.6f} {:.6f} {:.6f}\n".format(*row))
    file_handle.write("AtomData:  id type       cartes_x      cartes_y      cartes_z\n")

    li_atoms = [atom for atom in atoms if atom.symbol == "Li"]
    c_atoms = [atom for atom in atoms if atom.symbol != "Li"]
    sorted_atoms = li_atoms + c_atoms

    for new_id, atom in enumerate(sorted_atoms, start=1):
        atom_type = 0 if atom.symbol == 'Li' else 1
        x, y, z = atom.position
        file_handle.write(f"{new_id:>8} {atom_type:>4} {x:12.6f} {y:12.6f} {z:12.6f}\n")

    if energy is not None:
        file_handle.write("Energy\n  {:.8f}\n".format(energy))

    if swarmindex is not None:
        file_handle.write(f"Feature   Structure_{swarmindex}    Generation\n")
    else:
        file_handle.write("Feature   Structure    Generation\n")
    file_handle.write("END_CFG\n")

def generate_li_variants(
    initial_indices,
    valid_sites,
    target_n_li,
    *,
    # ▼ 새 옵션들 (기본값은 다양성↑ 쪽으로)
    max_variants=None,          # 전체 최대 개수 제한(없으면 생성되는 대로)
    n_random=10,                # 순수 랜덤 제거 기반 추가 생성 개수
    n_layerwise=10,             # 층별(partial) 제거 변형 개수
    n_alt_odd=10,               # 홀수층 우선 제거 변형 개수
    n_fullwipe=10,              # 아래층부터 통째로 지우는 변형 개수
    n_fullwipe_odd=10,          # 홀수층 통째 제거 변형 개수
    n_elite_perturb=20,         # “좋은” 패턴(초기/초기트림)에서 1–2자리만 치환
    # ▼ 새로 추가한 variant 타입
    n_facing=4,                 # facing configuration 변형 개수
    n_pair=4,                   # pair configuration 변형 개수
    min_symdiff=2,              # 기존 variant들과의 대칭차(해밍거리) 최소값
    always_perturb_if_equal=True,  # 초기 개수가 target과 같아도 파생을 더 뽑을지
    rng=None
):
    """
    initial_indices(인덱스 리스트)에서 target_n_li로 줄이는 다양한 variant를 만든다.
    - 기존 5가지 전략 + elite perturbation 유지
    - 추가: facing configuration, pair configuration
    """
    import numpy as np
    from collections import defaultdict
    from scipy.spatial import cKDTree

    if rng is None:
        rng = np.random.default_rng()

    def group_by_layer(indices):
        z_dict = defaultdict(list)
        for idx in indices:
            z = round(float(valid_sites[idx][2]), 3)
            z_dict[z].append(idx)
        # z 오름차순(작은 z → 큰 z). "맨 위"는 가장 큰 z.
        return dict(sorted(z_dict.items(), key=lambda kv: kv[0]))

    def all_indices_grouped_by_layer():
        """valid_sites 전체를 z로 그룹화"""
        z_dict = defaultdict(list)
        for idx in range(len(valid_sites)):
            z = round(float(valid_sites[idx][2]), 3)
            z_dict[z].append(idx)
        return dict(sorted(z_dict.items(), key=lambda kv: kv[0]))

    def layer_xy_maps(layer_indices):
        """
        한 층의 인덱스들을 XY 키로 맵핑.
        - 정확 매칭: (round(x,3), round(y,3))
        - 근사 매칭 대비를 위해 KDTree도 같이 반환
        """
        xs = []
        ys = []
        key_to_idx = {}
        for idx in layer_indices:
            x, y = float(valid_sites[idx][0]), float(valid_sites[idx][1])
            key = (round(x, 3), round(y, 3))
            # 같은 XY키에 여러 사이트가 있을 수 있지만 보통 1개로 가정
            # 여러 개면 첫 번째만 사용
            if key not in key_to_idx:
                key_to_idx[key] = idx
            xs.append(x)
            ys.append(y)
        pts = np.column_stack([xs, ys]) if layer_indices else np.zeros((0, 2))
        kdt = cKDTree(pts) if len(pts) > 0 else None
        return key_to_idx, kdt, pts

    def match_xy_in_layer(key, layer_map, layer_kdt, layer_pts, tol=5e-3):
        """
        주어진 XY key를 그 층에서 찾음.
        - 1) 정확키가 있으면 그 인덱스
        - 2) 없으면 KDTree로 가장 가까운 점(거리 tol 이내) 근사 매칭
        - 3) 없으면 None
        """
        if key in layer_map:
            return layer_map[key]
        if layer_kdt is None or layer_pts.shape[0] == 0:
            return None
        q = np.array([key[0], key[1]])
        dist, j = layer_kdt.query(q)
        if np.isfinite(dist) and dist <= tol:
            # j번째 포인트에 해당하는 실제 인덱스 찾기
            # layer_pts[j] == XY 좌표인 인덱스를 역추적
            xj, yj = layer_pts[j]
            k2 = (round(float(xj), 3), round(float(yj), 3))
            return layer_map.get(k2, None)
        return None

    def remove_random(indices, k):
        return rng.choice(indices, size=len(indices) - k, replace=False)

    def symdiff_size(a, b):
        # 해셋 차이 크기 (가까운 중복 억제용)
        return len(set(a) ^ set(b))

    def add_variant(cand):
        nonlocal variants
        key = tuple(sorted(cand))
        # 정확 중복 제거
        if key in keys:
            return False
        # 근사 중복 제거(해밍거리)
        if min_symdiff is not None and min_symdiff > 0:
            for v in variants:
                if symdiff_size(v, cand) < min_symdiff:
                    return False
        variants.append(np.array(cand, dtype=int))
        keys.add(key)
        return True

    variants = []
    keys = set()

    initial_indices = list(map(int, initial_indices))
    initial_n_li = len(initial_indices)
    if initial_n_li < target_n_li:
        raise ValueError(f"[ERROR] Initial Li count ({initial_n_li}) is less than target ({target_n_li}).")

    # --- 베이스: 초기 패턴을 섞고 target 길이로 트림한 것 1개 ---
    idx0 = initial_indices.copy()
    rng.shuffle(idx0)
    base_trim = np.array(idx0[:target_n_li], dtype=int)
    add_variant(base_trim)

    # 초기 개수 == target 인 경우: 옵션에 따라 섭동만/종료
    if initial_n_li == target_n_li and not always_perturb_if_equal:
        # print(f("[INFO] Generated {len(variants)} Li site variants from POSCAR-based initial_indices."))
        return variants

    # 공통 준비물
    k_remove = initial_n_li - target_n_li
    grouped_initial = group_by_layer(initial_indices)
    z_layers_initial = list(grouped_initial.keys())

    # ---- 기존 전략 1: Random removal ----
    for _ in range(n_random):
        trimmed = remove_random(initial_indices, k_remove)
        add_variant(sorted(trimmed[:target_n_li]))
        if max_variants and len(variants) >= max_variants:
            break

    # ---- 기존 전략 2: Layer-wise partial removal (bottom-up) ----
    for _ in range(n_layerwise):
        temp = initial_indices.copy()
        rng.shuffle(temp)
        removed = 0
        for z in z_layers_initial:
            if removed >= k_remove: break
            layer = grouped_initial[z]
            n = min(k_remove - removed, len(layer))
            selected = rng.choice(layer, n, replace=False).tolist()
            for s in selected:
                if s in temp:
                    temp.remove(s)
            removed += n
        if len(temp) >= target_n_li:
            rng.shuffle(temp)
            add_variant(sorted(temp[:target_n_li]))
        else:
            print(f"[WARN] Variant2: Skipped variant due to insufficient sites: len(temp)={len(temp)}")
        if max_variants and len(variants) >= max_variants:
            break

    # ---- 기존 전략 3: Alternating layer removal (odd layers) ----
    for _ in range(n_alt_odd):
        temp = initial_indices.copy()
        rng.shuffle(temp)
        removed = 0
        for z in z_layers_initial[1::2]:
            if removed >= k_remove: break
            layer = grouped_initial[z]
            n = min(k_remove - removed, len(layer))
            selected = rng.choice(layer, n, replace=False).tolist()
            for s in selected:
                if s in temp:
                    temp.remove(s)
            removed += n
        if len(temp) >= target_n_li:
            rng.shuffle(temp)
            add_variant(sorted(temp[:target_n_li]))
        else:
            print(f"[WARN] Variant3: Skipped variant due to insufficient sites: len(temp)={len(temp)}")
        if max_variants and len(variants) >= max_variants:
            break

    # ---- 기존 전략 4: Full layer wipe (bottom-up) ----
    for _ in range(n_fullwipe):
        temp = initial_indices.copy()
        rng.shuffle(temp)
        removed = 0
        for z in z_layers_initial:
            layer = grouped_initial[z]
            if removed + len(layer) <= k_remove:
                for s in layer:
                    if s in temp:
                        temp.remove(s)
                removed += len(layer)
            else:
                n = k_remove - removed
                selected = rng.choice(layer, n, replace=False).tolist()
                for s in selected:
                    if s in temp:
                        temp.remove(s)
                break
        if len(temp) >= target_n_li:
            add_variant(sorted(temp[:target_n_li]))
        else:
            print(f"[WARN] Variant4: Skipped variant due to insufficient sites: len(temp)={len(temp)}")
        if max_variants and len(variants) >= max_variants:
            break

    # ---- 기존 전략 5: Full layer wipe (odd layers only) ----
    for _ in range(n_fullwipe_odd):
        temp = initial_indices.copy()
        removed = 0
        grouped_temp = group_by_layer(temp)
        z_layers_temp = list(grouped_temp.keys())
        for z in z_layers_temp[1::2]:
            layer = grouped_temp[z]
            if removed + len(layer) <= k_remove:
                for s in layer:
                    if s in temp:
                        temp.remove(s)
                removed += len(layer)
            else:
                n = k_remove - removed
                selected = rng.choice(layer, n, replace=False).tolist()
                for s in selected:
                    if s in temp:
                        temp.remove(s)
                break
        if len(temp) >= target_n_li:
            rng.shuffle(temp)
            add_variant(sorted(temp[:target_n_li]))
        else:
            print(f"[WARN] Variant5: Skipped variant due to insufficient sites: len(temp)={len(temp)}")
        if max_variants and len(variants) >= max_variants:
            break

    # =======================
    # 새 전략 A: Facing configuration
    # - per_layer_max = floor(target_n_li / 4) (최소 1)
    # - 맨 위(z 최대) 층에 per_layer_max 배치
    # - 이후 다른 층들에도 같은 (x,y) 좌표로 배치하면서 target_n_li 채움
    # =======================
    all_layers = all_indices_grouped_by_layer()
    layer_zs_sorted = list(all_layers.keys())            # 오름차순
    if len(layer_zs_sorted) >= 1:
        top_z = layer_zs_sorted[-1]                      # 맨 위 층
        top_layer_idxs = all_layers[top_z]
        top_map, top_kdt, top_pts = layer_xy_maps(top_layer_idxs)

        # top layer의 XY 키 리스트
        top_keys = list(top_map.keys())
        rng.shuffle(top_keys)

        per_layer_max = max(1, target_n_li // len(layer_zs_sorted))

        for _ in range(n_facing):
            cand = []
            # 1) 맨 위 층에 우선 배치
            chosen_keys = top_keys[:min(per_layer_max, len(top_keys))]
            for key in chosen_keys:
                cand.append(top_map[key])

            # 2) 아래(=z 작은) 층들에 동일 XY로 채우기 (필요 수만큼)
            if len(cand) < target_n_li:
                # top 아래로 내려가며 동일 XY key로 매칭
                # 오름차순 z에서 top 아래는 끝에서 두번째, 세번째 ...
                for z in reversed(layer_zs_sorted[:-1]):
                    layer_idxs = all_layers[z]
                    layer_map, layer_kdt, layer_pts = layer_xy_maps(layer_idxs)
                    for key in chosen_keys:
                        if len(cand) >= target_n_li:
                            break
                        idx = match_xy_in_layer(key, layer_map, layer_kdt, layer_pts, tol=5e-3)
                        if idx is not None:
                            cand.append(idx)
                    if len(cand) >= target_n_li:
                        break

            # 길이 정합
            if len(cand) >= target_n_li:
                cand = cand[:target_n_li]
                add_variant(sorted(cand))
            else:
                print(f"[WARN] Facing: could not reach target_n_li (got {len(cand)}/{target_n_li}).")

            if max_variants and len(variants) >= max_variants:
                break

    # =======================
    # 새 전략 B: Pair configuration
    # - per_layer_max = floor(target_n_li / 4) (최소 1)
    # - L1(맨 위): per_layer_max
    # - L3(세 번째 Li 층): L1과 다른 XY로 per_layer_max
    # - L2: L1과 동일 XY
    # - L4: L3과 동일 XY
    # =======================
    if len(layer_zs_sorted) >= 4:
        z1, z2, z3, z4 = layer_zs_sorted[-1], layer_zs_sorted[-2], layer_zs_sorted[-3], layer_zs_sorted[-4]
        L1_idxs, L2_idxs, L3_idxs, L4_idxs = all_layers[z1], all_layers[z2], all_layers[z3], all_layers[z4]
        L1_map, L1_kdt, L1_pts = layer_xy_maps(L1_idxs)
        L2_map, L2_kdt, L2_pts = layer_xy_maps(L2_idxs)
        L3_map, L3_kdt, L3_pts = layer_xy_maps(L3_idxs)
        L4_map, L4_kdt, L4_pts = layer_xy_maps(L4_idxs)

        keys_L1 = list(L1_map.keys())
        keys_L3 = list(L3_map.keys())
        rng.shuffle(keys_L1)
        rng.shuffle(keys_L3)

        per_layer_max = max(1, target_n_li // len(layer_zs_sorted))

        for _ in range(n_pair):
            cand = []

            # 1) L1에서 per_layer_max 선택
            K1 = keys_L1[:min(per_layer_max, len(keys_L1))]
            for key in K1:
                cand.append(L1_map[key])

            # 2) L3에서 K1과 겹치지 않는 XY로 per_layer_max 선택
            K3 = []
            for key in keys_L3:
                if len(K3) >= per_layer_max:
                    break
                if key not in set(K1):
                    K3.append(key)
            # 만약 충분치 않으면 근사라도 채워보자(겹치더라도 허용)
            i3 = 0
            while len(K3) < per_layer_max and i3 < len(keys_L3):
                if keys_L3[i3] not in K3:
                    K3.append(keys_L3[i3])
                i3 += 1

            for key in K3[:min(per_layer_max, len(K3))]:
                cand.append(L3_map[key])

            # 3) L2에 K1 동일 XY
            for key in K1:
                if len(cand) >= target_n_li:
                    break
                idx = match_xy_in_layer(key, L2_map, L2_kdt, L2_pts, tol=5e-3)
                if idx is not None:
                    cand.append(idx)

            # 4) L4에 K3 동일 XY
            for key in K3:
                if len(cand) >= target_n_li:
                    break
                idx = match_xy_in_layer(key, L4_map, L4_kdt, L4_pts, tol=5e-3)
                if idx is not None:
                    cand.append(idx)

            # 길이 정합
            if len(cand) >= target_n_li:
                cand = cand[:target_n_li]
                add_variant(sorted(cand))
            else:
                print(f"[WARN] Pair: could not reach target_n_li (got {len(cand)}/{target_n_li}).")

            if max_variants and len(variants) >= max_variants:
                break
    else:
        if n_pair > 0:
            print(f"[WARN] Pair: Not enough layers ({len(layer_zs_sorted)} < 4). Skipped.")

    # ---- 추가 전략: Elite perturbation (소규모 자리 치환 1~2개) ----
    def small_perturb(base):
        base = list(map(int, base))
        k = rng.integers(1, 3)  # 1~2자리
        pool = list(set(range(len(valid_sites))) - set(base))
        for _ in range(k):
            if not pool: break
            pos = rng.integers(len(base))
            base[pos] = int(rng.choice(pool))
        return sorted(base)

    seed_pool = [base_trim.tolist()]
    seed_pool += [v.tolist() for v in variants[:max(1, len(variants)//4)]]

    for _ in range(n_elite_perturb):
        src = seed_pool[rng.integers(len(seed_pool))]
        cand = small_perturb(src)
        add_variant(cand)
        if max_variants and len(variants) >= max_variants:
            break

#    print(f"[INFO] Generated {len(variants)} Li site variants from POSCAR-based initial_indices "
#          f"(min_symdiff={min_symdiff}).")
    return variants


def generate_li_variants_budgeted(
    initial_indices, valid_sites, target_n_li, *,
    budget,                    # 이번 iteration에서 필요한 variant 개수
    weights=None,              # 전략별 비율. 합이 1.0이면 좋음. None이면 균등
    min_symdiff=2,             # 유사도 필터
    rng=None
):
    """
    budget 만큼만 variant를 만들어 반환.
    내부적으로 여러 전략을 섞어 쓰고, budget을 채우면 즉시 종료.
    """
    import numpy as np
    if rng is None:
        rng = np.random.default_rng()

    # 기본 비중 (합≈1.0). 필요시 수정 가능
    if weights is None:
        weights = dict(
            random=0.20,
            layerwise=0.15,
            alt_odd=0.15,
            fullwipe=0.15,
            fullwipe_odd=0.15,
            facing=0.1,        # 새 전략 가중치
            pair=0.1,          # 새 전략 가중치
        )

    # 먼저 “씨앗” 세트(초기 트림 한 개)는 항상 포함
    base_variants = generate_li_variants(
        initial_indices, valid_sites, target_n_li,
        max_variants=None,
        n_random=1, n_layerwise=0, n_alt_odd=0, n_fullwipe=0, n_fullwipe_odd=0,
        n_elite_perturb=0, n_facing=0, n_pair=0,
        min_symdiff=min_symdiff, always_perturb_if_equal=True, rng=rng
    )
    seeds = [np.array(v, dtype=int) for v in base_variants]  # 보통 1개

    # 전략별 목표 개수 계산
    k = max(0, budget - len(seeds))
    if k == 0:
        return seeds[:budget]

    labels, probs = zip(*weights.items())
    alloc = [int(round(p * k)) for p in probs]
    diff = k - sum(alloc)
    order = np.argsort(probs)[::-1]
    for i in range(abs(diff)):
        alloc[order[i % len(alloc)]] += 1 if diff > 0 else -1
    target_per_strategy = dict(zip(labels, alloc))

    results = list(seeds)
    seen = {tuple(sorted(v)) for v in results}

    def add(v):
        key = tuple(sorted(v))
        if key in seen:
            return False
        if min_symdiff:
            for u in results:
                if len(set(u) ^ set(v)) < min_symdiff:
                    return False
        results.append(np.array(v, dtype=int))
        seen.add(key)
        return True

    def gen_random(cnt):
        return generate_li_variants(initial_indices, valid_sites, target_n_li,
                                    n_random=cnt, n_layerwise=0, n_alt_odd=0,
                                    n_fullwipe=0, n_fullwipe_odd=0, n_elite_perturb=0,
                                    n_facing=0, n_pair=0,
                                    max_variants=None, min_symdiff=min_symdiff, rng=rng)

    def gen_layerwise(cnt):
        return generate_li_variants(initial_indices, valid_sites, target_n_li,
                                    n_random=0, n_layerwise=cnt, n_alt_odd=0,
                                    n_fullwipe=0, n_fullwipe_odd=0, n_elite_perturb=0,
                                    n_facing=0, n_pair=0,
                                    max_variants=None, min_symdiff=min_symdiff, rng=rng)

    def gen_alt_odd(cnt):
        return generate_li_variants(initial_indices, valid_sites, target_n_li,
                                    n_random=0, n_layerwise=0, n_alt_odd=cnt,
                                    n_fullwipe=0, n_fullwipe_odd=0, n_elite_perturb=0,
                                    n_facing=0, n_pair=0,
                                    max_variants=None, min_symdiff=min_symdiff, rng=rng)

    def gen_fullwipe(cnt):
        return generate_li_variants(initial_indices, valid_sites, target_n_li,
                                    n_random=0, n_layerwise=0, n_alt_odd=0,
                                    n_fullwipe=cnt, n_fullwipe_odd=0, n_elite_perturb=0,
                                    n_facing=0, n_pair=0,
                                    max_variants=None, min_symdiff=min_symdiff, rng=rng)

    def gen_fullwipe_odd(cnt):
        return generate_li_variants(initial_indices, valid_sites, target_n_li,
                                    n_random=0, n_layerwise=0, n_alt_odd=0,
                                    n_fullwipe=0, n_fullwipe_odd=cnt, n_elite_perturb=0,
                                    n_facing=0, n_pair=0,
                                    max_variants=None, min_symdiff=min_symdiff, rng=rng)

    def gen_facing(cnt):
        return generate_li_variants(initial_indices, valid_sites, target_n_li,
                                    n_random=0, n_layerwise=0, n_alt_odd=0,
                                    n_fullwipe=0, n_fullwipe_odd=0, n_elite_perturb=0,
                                    n_facing=cnt, n_pair=0,
                                    max_variants=None, min_symdiff=min_symdiff, rng=rng)

    def gen_pair(cnt):
        return generate_li_variants(initial_indices, valid_sites, target_n_li,
                                    n_random=0, n_layerwise=0, n_alt_odd=0,
                                    n_fullwipe=0, n_fullwipe_odd=0, n_elite_perturb=0,
                                    n_facing=0, n_pair=cnt,
                                    max_variants=None, min_symdiff=min_symdiff, rng=rng)

    generators = dict(
        random=gen_random,
        layerwise=gen_layerwise,
        alt_odd=gen_alt_odd,
        fullwipe=gen_fullwipe,
        fullwipe_odd=gen_fullwipe_odd,
        facing=gen_facing,
        pair=gen_pair,
    )

    # 라운드로빈으로 각 전략에서 조금씩 가져오며 budget 채우기
    while len(results) < budget:
        progressed = False
        for name in labels:
            need = target_per_strategy.get(name, 0)
            if need <= 0:
                continue
            batch = min(3, need)
            vs = generators[name](batch)
            cnt_add = 0
            for v in vs:
                if add(v):
                    cnt_add += 1
                if len(results) >= budget:
                    break
            target_per_strategy[name] -= max(0, cnt_add)
            if cnt_add > 0:
                progressed = True
            if len(results) >= budget:
                break
        if not progressed:
            break

    return results[:budget]

# ----------------------------
# 평가 함수 (파일락 + 임시파일 고유화)
# ----------------------------
def energy_function_discrete_index(x, structure, n_li, r_min, z_layer_tol,
                                   pso_io_dir, valid_sites, max_per_layer, cell,
                                   success_counter=None, iteration=None, swarmindex=None):

    indices = np.round(x).astype(int)
    if len(set(indices)) < len(indices):
        return 1e5

    selected_sites = np.array([valid_sites[i] for i in indices])

    for i in range(n_li):
        for j in range(i + 1, n_li):
            dz = abs(selected_sites[i][2] - selected_sites[j][2])
            d  = np.linalg.norm(selected_sites[i] - selected_sites[j])
            if dz < z_layer_tol and d < r_min:
                return 1e5

    atoms_with_li = add_lithium_atoms(structure, selected_sites, valid_sites, snap=False)

    iter_str = f"{iteration:03d}" if iteration is not None else "000"
    in_name  = f"IN_{iter_str}.cfg"
    out_name = f"OUT_{iter_str}.cfg"
    in_path  = os.path.join(pso_io_dir, in_name)
    out_path = os.path.join(pso_io_dir, out_name)

    uid = f"{os.getpid()}_{uuid.uuid4().hex}"
    tmp_input  = os.path.join(pso_io_dir, f"__tmp_input_{iter_str}_{swarmindex}_{uid}.cfg")
    tmp_output = os.path.join(pso_io_dir, f"__tmp_output_{iter_str}_{swarmindex}_{uid}.cfg")

    _safe_append(in_path, "")  # ensure file exists
    with open(tmp_input, 'w') as f:
        write_mtp_cfg(atoms_with_li, f, swarmindex=swarmindex)

    with open(tmp_input, 'r') as rf:
        _safe_append(in_path, rf.read() + "\n")

    try:
        cmd = f"mlp calc-efs pot.mtp {tmp_input} {tmp_output}"
        subprocess.run(cmd, shell=True, check=True)

        success_file = os.path.join(pso_io_dir, f"IN_{iter_str}_success.cfg")
        with open(tmp_input, 'r') as rf:
            _safe_append(success_file, rf.read() + "\n")

        with open(tmp_output, 'r') as rf:
            relaxed_block = rf.read()
        _safe_append(out_path, relaxed_block)

        lines = relaxed_block.splitlines()
        for i, line in enumerate(lines):
            if "Energy" in line:
                try:
                    return float(lines[i + 1].strip())
                except (IndexError, ValueError):
                    return 1e6
    except subprocess.CalledProcessError:
        return 1e6
    finally:
        for p in (tmp_input, tmp_output):
            try:
                if os.path.exists(p):
                    os.remove(p)
            except:
                pass

# ----------------------------
# CFG 파서 / 저장 / 도우미들
# ----------------------------
def split_cfg_blocks(cfg_path):
    with open(cfg_path, 'r') as f:
        lines = f.readlines()
    blocks, current, inside = [], [], False
    for line in lines:
        if 'BEGIN_CFG' in line:
            inside = True
            current = [line]
        elif 'END_CFG' in line:
            current.append(line)
            blocks.append(''.join(current))
            inside = False
        elif inside:
            current.append(line)
    return blocks

def parse_cfg_block(block_str):
    lines = block_str.strip().splitlines()
    cell, positions, symbols = [], [], []
    read_cell = False
    read_atoms = False
    cell_lines_remaining = 0
    atom_lines_remaining = 0
    size_found = False
    for i, line in enumerate(lines):
        line = line.strip()
        if not line:
            continue
        if line.startswith("Size") and not size_found:
            try:
                atom_lines_remaining = int(lines[i + 1].strip())
                size_found = True
            except:
                raise ValueError("Failed to read Size from CFG block.")
        elif line.startswith("Supercell"):
            read_cell = True
            cell_lines_remaining = 3
            continue
        elif read_cell and cell_lines_remaining > 0:
            cell.append([float(x) for x in line.split()])
            cell_lines_remaining -= 1
            continue
        elif line.startswith("AtomData"):
            read_atoms = True
            continue
        elif read_atoms and atom_lines_remaining > 0:
            tokens = line.split()
            if len(tokens) < 5:
                continue
            _, atom_type, x, y, z = tokens[:5]
            positions.append([float(x), float(y), float(z)])
            symbols.append(int(atom_type))
            atom_lines_remaining -= 1
            if atom_lines_remaining == 0:
                read_atoms = False
    Z_map = {0: "Li", 1: "C"}
    elements = [Z_map.get(t, "X") for t in symbols]
    return Atoms(symbols=elements, positions=positions, cell=np.array(cell), pbc=True)

def find_best_structure_index(convergence_file):
    best_E = float("inf")
    best_iter, best_swarm = None, None
    with open(convergence_file, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split()
            if len(parts) > 3:
                continue
            iter_i, swarm_i, E = int(parts[0]), int(parts[1]), float(parts[2])
            if E < best_E:
                best_E = E
                best_iter = iter_i
                best_swarm = swarm_i
    return best_iter, best_swarm, best_E

def save_best_structure_from_convergence(cfg_dir, convergence_file, output_dir):
    best_iter, best_swarm, best_E = find_best_structure_index(convergence_file)
    print(f"[INFO] Best structure: Iter={best_iter}, Swarm={best_swarm}, E={best_E:.6f}", flush=True)

    in_file = os.path.join(cfg_dir, f"IN_{best_iter:03d}_relax_success.cfg")
    out_file = os.path.join(cfg_dir, f"OUT_{best_iter:03d}_relax.cfg")

    if not os.path.exists(in_file) or not os.path.exists(out_file):
        print(f"[WARNING] Missing IN/OUT file for iter {best_iter}", flush=True)
        return

    in_blocks = split_cfg_blocks(in_file)
    out_blocks = split_cfg_blocks(out_file)

    initial_atoms = parse_cfg_block(in_blocks[0])  # Relax에서 첫 블록이 top
    relaxed_atoms = parse_cfg_block(out_blocks[0])

    prefix = os.path.join(output_dir, f"best_iter{best_iter:03d}_swarm{best_swarm}_E{best_E:.6f}")
    write(prefix + "_initial.cif", initial_atoms, format="cif")
    write(prefix + "_relaxed.cif", relaxed_atoms, format="cif")
    print(f"[INFO] Saved best structure to:\n  {prefix}_initial.cif\n  {prefix}_relaxed.cif", flush=True)

def save_convergence_plot(history, all_energies, output_dir, filename="convergence_plot.png", data_filename=None):
    """
    PSO 에너지 수렴 그래프 저장
    - Best So Far 절대 상승 금지
    - New Best는 갱신 시점에만 별표
    """

    # === 각 iteration별 최소값 추출 ===
    per_iter_min = {}
    try:
        # all_energies = [(iter, swarm_idx, val), ...]
        for it, sw, val in all_energies:
            if val < 1e5:
                per_iter_min[it] = min(per_iter_min.get(it, np.inf), float(val))
        use_swarm_index = True
    except ValueError:
        # all_energies = [(iter, val), ...]
        for it, val in all_energies:
            if val < 1e5:
                per_iter_min[it] = min(per_iter_min.get(it, np.inf), float(val))
        use_swarm_index = False

    # === Iteration 정렬 후 최소값 배열 생성 ===
    iters_sorted = sorted(per_iter_min.keys())
    vals_min = np.array([per_iter_min[i] for i in iters_sorted], dtype=float)

    # === 누적 최소값(cummin) 계산 (절대 상승 금지) ===
    best_so_far = np.minimum.accumulate(vals_min)
    for k in range(1, len(best_so_far)):
        if best_so_far[k] > best_so_far[k-1]:
            best_so_far[k] = best_so_far[k-1]

    # === New Best 포인트만 추출 ===
    stars_x, stars_y = [], []
    cur = np.inf
    for i, v in zip(iters_sorted, best_so_far):
        if v < cur - 1e-12:  # 수치 오차 보정
            stars_x.append(i)
            stars_y.append(v)
            cur = v

    # === 산점도 데이터 ===
    try:
        gen_iters = [i for i, j, val in all_energies if val < 1e5]
        gen_vals  = [val for i, j, val in all_energies if val < 1e5]
    except ValueError:
        gen_iters = [i for i, val in all_energies if val < 1e5]
        gen_vals  = [val for i, val in all_energies if val < 1e5]

    # === 플롯 ===
    plt.figure()
    plt.scatter(gen_iters, gen_vals, label="All Candidates", alpha=0.5, s=10)
    plt.step(iters_sorted, best_so_far, where='post', color='red', linewidth=2, label="Best So Far")
    if stars_x:
        plt.plot(stars_x, stars_y, '*', color='gold', markersize=12, label="New Best")

    plt.xlabel("Iteration")
    plt.ylabel("Energy (eV)")
    plt.title("PSO Energy Convergence")
    plt.grid(True)
    plt.legend()

    os.makedirs(output_dir, exist_ok=True)
    path = os.path.join(output_dir, filename)
    plt.savefig(path, dpi=300)
    plt.close()
    print(f"[INFO] Convergence plot saved to: {path}", flush=True)

    # === 데이터 파일 저장 ===
    if data_filename:
        data_path = os.path.join(output_dir, data_filename)
        with open(data_path, "w") as f:
            if use_swarm_index:
                f.write("# Iteration   SwarmIndex   Energy\n")
                for iter_i, swarm_i, energy in all_energies:
                    f.write(f"{iter_i}   {swarm_i}   {energy:.6f}\n")
            else:
                f.write("# Iteration   Energy\n")
                for iter_i, energy in all_energies:
                    f.write(f"{iter_i}   {energy:.6f}\n")
        print(f"[INFO] Convergence data saved to: {data_path}", flush=True)

def save_convergence_plot_OLD(history, all_energies, output_dir, filename="convergence_plot.png", data_filename=None):
    iters, best_vals = zip(*history)
    try:
        gen_iters = [i for i, j, val in all_energies if val < 1e5]
        gen_vals  = [val for i, j, val in all_energies if val < 1e5]
        use_swarm_index = True
    except ValueError:
        gen_iters = [i for i, val in all_energies if val < 1e5]
        gen_vals  = [val for i, val in all_energies if val < 1e5]
        use_swarm_index = False

    iters_filtered, best_vals_filtered = [], []
    for i, v in zip(iters, best_vals):
        if v < 1e5:
            iters_filtered.append(i)
            best_vals_filtered.append(v)

    star_iters, star_vals = [], []
    best_so_far = float("inf")
    for i, v in zip(iters_filtered, best_vals_filtered):
        if v < best_so_far:
            best_so_far = v
            star_iters.append(i)
            star_vals.append(v)

    plt.figure()
    plt.scatter(gen_iters, gen_vals, label="All Candidates", alpha=0.5, s=10)
    plt.plot(iters_filtered, best_vals_filtered, '-', color='red', label="Best So Far", linewidth=2)
    plt.plot(star_iters, star_vals, '*', color='gold', markersize=12, label="New Best")
    plt.xlabel("Iteration")
    plt.ylabel("Energy (eV)")
    plt.title("PSO Energy Convergence")
    plt.grid(True)
    plt.legend()
    os.makedirs(output_dir, exist_ok=True)
    path = os.path.join(output_dir, filename)
    plt.savefig(path, dpi=300)
    plt.close()
    print(f"[INFO] Convergence plot saved to: {path}", flush=True)

    if data_filename:
        data_path = os.path.join(output_dir, data_filename)
        with open(data_path, "w") as f:
            if use_swarm_index:
                f.write("# Iteration   SwarmIndex   Energy\n")
                for iter_i, swarm_i, energy in all_energies:
                    f.write(f"{iter_i}   {swarm_i}   {energy:.6f}\n")
            else:
                f.write("# Iteration   Energy\n")
                for iter_i, energy in all_energies:
                    f.write(f"{iter_i}   {energy:.6f}\n")
        print(f"[INFO] Convergence data saved to: {data_path}", flush=True)

# ----------------------------
# Iteration 후 Relax 단계 (기존과 동일 로직)
# ----------------------------
def post_iteration_relaxation(iteration, output_dir, valid_sites):
    import heapq
    tmp_input = os.path.join(output_dir, "__tmp_relax_input.cfg")
    tmp_output = os.path.join(output_dir, "__tmp_relax_output.cfg")
    convergence_file = os.path.join(output_dir, "convergence_data.txt")
    in_file = os.path.join(output_dir, f"IN_{iteration:03d}_success.cfg")
    out_file = os.path.join(output_dir, f"OUT_{iteration:03d}.cfg")

    blocks = split_cfg_blocks(in_file)
    seen_energy = set()
    swarm_energy = []
    unique_structures = set()

    with open(convergence_file, "r") as f:
        for line in f:
            if line.startswith("#"): continue
            iter_i, swarm_i, energy = line.strip().split()[:3]
            iter_i, swarm_i = int(iter_i), int(swarm_i)
            energy = float(energy)
            if iter_i == iteration:
                if energy < 1e5 and energy not in seen_energy:
                    seen_energy.add(energy)
                    key = (round(energy, 6),)
                    if key not in unique_structures:
                        unique_structures.add(key)
                        swarm_energy.append((energy, swarm_i))

    top_k = max(1, len(swarm_energy) // 20)   # 5 %
    top_swarm = [i for _, i in heapq.nsmallest(top_k, swarm_energy)]

    in_relax = os.path.join(output_dir, f"IN_{iteration:03d}_relax.cfg")
    out_relax = os.path.join(output_dir, f"OUT_{iteration:03d}_relax.cfg")
    in_success = os.path.join(output_dir, f"IN_{iteration:03d}_relax_success.cfg")

    with open(in_relax, "w") as f:
        for idx in top_swarm:
            atoms = parse_cfg_block(blocks[idx])
            write_mtp_cfg(atoms, f, swarmindex=idx)
            f.write("\n")

    with open(in_success, "w") as fsucc, open(out_relax, "w") as fout:
        for idx in top_swarm:
            atoms = parse_cfg_block(blocks[idx])
            with open(tmp_input, "w") as ftmp:
                write_mtp_cfg(atoms, ftmp, swarmindex=idx)
            cmd = f"mlp relax mlip.ini --cfg-filename={tmp_input} --save-relaxed={tmp_output} --force-tolerance=0.04 > relax-output.txt"
            try:
                subprocess.run(cmd, shell=True, check=True)
                with open(tmp_output + "_0", 'r') as fin:
                    relaxed_block = fin.read()
                fout.write(relaxed_block + "\n")
                write_mtp_cfg(atoms, fsucc, swarmindex=idx)
                fsucc.write("\n")
                lines = relaxed_block.splitlines()
                for i, line in enumerate(lines):
                    if line.strip() == "Energy" and i + 1 < len(lines):
                        try:
                            energy = float(lines[i + 1].strip())
                            with open(convergence_file, "a") as fc:
                                fc.write(f"{iteration}   {idx}   {energy:.6f}  # relaxed\n")
                        except ValueError:
                            print(f"[ERROR] Failed to parse energy value for swarm {idx} at line {i+1}: '{lines[i + 1]}'", flush=True)
                        break
            except Exception as e:
                print(f"[ERROR] Relaxation failed for swarm {idx}: {e}", flush=True)
###########
# 병렬 계산 #
###########
def post_iteration_relaxation_parallel(iteration, output_dir, valid_sites, relax_workers=None):
    """
    iteration의 상위 후보들을 mlp relax로 병렬 처리.
    모든 relax 완료 후에 OUT_{iter}_relax.cfg, IN_{iter}_relax_success.cfg, convergence_data.txt 를 한 번에 갱신.
    반환: 성공 개수(int)
    """
    import heapq, os, shutil, subprocess, uuid
    from io import StringIO
    from concurrent.futures import ThreadPoolExecutor, as_completed

    # 경로
    convergence_file = os.path.join(output_dir, "convergence_data.txt")
    in_file_success  = os.path.join(output_dir, f"IN_{iteration:03d}_success.cfg")
    in_relax         = os.path.join(output_dir, f"IN_{iteration:03d}_relax.cfg")
    out_relax        = os.path.join(output_dir, f"OUT_{iteration:03d}_relax.cfg")
    in_relax_success = os.path.join(output_dir, f"IN_{iteration:03d}_relax_success.cfg")
    logs_dir         = os.path.join(output_dir, "relax_logs")

    # 로그 디렉토리 리셋
    if os.path.exists(logs_dir):
        shutil.rmtree(logs_dir)
    os.makedirs(logs_dir, exist_ok=True)

    # 입력 성공 블록 가드
    if not os.path.exists(in_file_success) or os.path.getsize(in_file_success) == 0:
        print(f"[INFO] No candidates to relax at iter {iteration} (missing or empty IN_success).")
        return 0

    # 후보 블록 로드 및 swarmindex 매핑
    blocks = split_cfg_blocks(in_file_success)

    def _extract_swarmindex(block_str: str):
        for line in block_str.splitlines():
            if "Feature" in line and "Structure_" in line:
                parts = line.replace("\t", " ").split()
                for tok in parts:
                    if tok.startswith("Structure_"):
                        try:
                            return int(tok.split("Structure_")[1])
                        except:
                            pass
        return None

    block_by_swarm = {}
    for b in blocks:
        si = _extract_swarmindex(b)
        if si is not None and si not in block_by_swarm:
            block_by_swarm[si] = b

    # 상위 후보 선별
    if not os.path.exists(convergence_file) or os.path.getsize(convergence_file) == 0:
        print(f"[INFO] No convergence data at iter {iteration}.")
        return 0

    seen_energy, unique_structures, swarm_energy = set(), set(), []
    with open(convergence_file, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split()
            if len(parts) < 3:
                continue
            iter_i, swarm_i, energy = int(parts[0]), int(parts[1]), float(parts[2])
            if iter_i == iteration and energy < 1e5 and energy not in seen_energy:
                seen_energy.add(energy)
                key = (round(energy, 6),)
                if key not in unique_structures:
                    unique_structures.add(key)
                    swarm_energy.append((energy, swarm_i))

    if not swarm_energy:
        print(f"[INFO] No candidates found for relaxation at iter {iteration}.")
        return 0

    top_k = max(1, len(swarm_energy) // 20)  # 상위 5%
    # block이 실제 존재하는 swarmindex만 채택
    top_swarm_all = [i for _, i in heapq.nsmallest(top_k, swarm_energy)]
    top_swarm = [i for i in top_swarm_all if i in block_by_swarm]

    if not top_swarm:
        print(f"[INFO] No top candidates for relax at iter {iteration} (no mapped blocks).")
        return 0

    # IN_{iter}_relax.cfg 생성 (입력 구조 기록)
    with open(in_relax, "w") as f:
        for idx in top_swarm:
            atoms = parse_cfg_block(block_by_swarm[idx])
            write_mtp_cfg(atoms, f, swarmindex=idx)
            f.write("\n")

    # --- 병렬 실행 준비 ---
    if relax_workers is None:
        relax_workers = max(1, (os.cpu_count() or 2) // 2)

    def _relax_one(idx):
        """단일 후보에 대해 mlp relax 실행하고, (idx, relaxed_block, energy, input_block_str) 반환"""
        uid = f"{os.getpid()}_{uuid.uuid4().hex}"
        tmp_input  = os.path.join(output_dir, f"__tmp_relax_input_{iteration:03d}_{idx}_{uid}.cfg")
        tmp_output = os.path.join(output_dir, f"__tmp_relax_output_{iteration:03d}_{idx}_{uid}.cfg")
        log_path   = os.path.join(logs_dir, f"relax_iter{iteration:03d}_idx{idx}_{uid}.log")

        blk = block_by_swarm.get(idx)
        if blk is None:
            with open(log_path, "w") as lf:
                lf.write(f"[WARN] No block mapped for swarmindex {idx}\n")
            return (idx, None, None, None)

        atoms = parse_cfg_block(blk)
        with open(tmp_input, "w") as ftmp:
            write_mtp_cfg(atoms, ftmp, swarmindex=idx)

        # mlip.ini 상대경로 보호
        ini_path = os.path.abspath("mlip.ini")
        ini_dir  = os.path.dirname(ini_path) if os.path.exists(ini_path) else None

        cmd = [
            "mlp", "relax", "mlip.ini",
            f"--cfg-filename={tmp_input}",
            f"--save-relaxed={tmp_output}",
            "--force-tolerance=0.02"
        ]
        try:
            with open(log_path, "w") as lf:
                proc = subprocess.run(
                    " ".join(cmd),
                    shell=True,
                    cwd=ini_dir if ini_dir else None,
                    stdout=lf,
                    stderr=lf,
                    check=False  # rc는 참고만
                )
                rc = proc.returncode

            # 일부 mlp는 _0 suffix 생성
            relaxed_path = tmp_output + "_0" if os.path.exists(tmp_output + "_0") else tmp_output

            # 성공 판정: 결과 파일 존재 + 크기 > 0
            if not os.path.exists(relaxed_path) or os.path.getsize(relaxed_path) == 0:
                with open(log_path, "a") as lf:
                    lf.write(f"\n[WARN] relaxed file missing/empty. rc={rc}, path={relaxed_path}\n")
                return (idx, None, None, None)

            with open(relaxed_path, 'r') as fin:
                relaxed_block = fin.read()

            # 에너지 파싱(없어도 OK)
            energy_val = None
            lines = relaxed_block.splitlines()
            for i, line in enumerate(lines):
                if line.strip() == "Energy" and i + 1 < len(lines):
                    try:
                        energy_val = float(lines[i + 1].strip())
                    except Exception:
                        energy_val = None
                    break

            # 입력 블록 문자열 저장(성공 기록용)
            buf = StringIO()
            write_mtp_cfg(atoms, buf, swarmindex=idx)
            input_block_str = buf.getvalue()

            with open(log_path, "a") as lf:
                lf.write(f"\n[OK] relaxed written. rc={rc}, size={os.path.getsize(relaxed_path)}\n")

            return (idx, relaxed_block, energy_val, input_block_str)

        except Exception as e:
            try:
                with open(log_path, "a") as lf:
                    lf.write(f"\n[EXC] {repr(e)}\n")
            except:
                pass
            return (idx, None, None, None)

        finally:
            for p in (tmp_input, tmp_output, tmp_output + "_0"):
                try:
                    if os.path.exists(p):
                        os.remove(p)
                except:
                    pass

    # --- 병렬 실행 ---
    results = []
    with ThreadPoolExecutor(max_workers=relax_workers) as ex:
        futs = [ex.submit(_relax_one, idx) for idx in top_swarm]
        for fut in as_completed(futs):
            results.append(fut.result())

    # --- 결과 쓰기(모아서 순서대로) ---
    succ = 0
    with open(out_relax, "w") as fout:
        for idx in top_swarm:
            recs = [r for r in results if r[0] == idx and r[1] is not None]
            if not recs:
                continue
            _, relaxed_block, _, _ = recs[0]
            fout.write(relaxed_block)
            if not relaxed_block.endswith("\n"):
                fout.write("\n")
            succ += 1

    with open(in_relax_success, "w") as fsucc:
        for idx in top_swarm:
            recs = [r for r in results if r[0] == idx and r[3] is not None]
            if not recs:
                continue
            _, _, _, input_block_str = recs[0]
            fsucc.write(input_block_str)
            if not input_block_str.endswith("\n"):
                fsucc.write("\n")

    # convergence_data.txt 에 relax 에너지 추가 기록
    if succ > 0:
        lines_to_append = []
        for idx in top_swarm:
            recs = [r for r in results if r[0] == idx and r[2] is not None]
            if not recs:
                continue
            _, _, energy_val, _ = recs[0]
            if energy_val is not None:
                lines_to_append.append(f"{iteration}   {idx}   {energy_val:.6f}  # relaxed\n")
        if lines_to_append:
            _safe_append(convergence_file, "".join(lines_to_append))

    print(f"[INFO] Parallel relaxation finished for iter {iteration} "
          f"({succ}/{len(top_swarm)} succeeded)", flush=True)
    return int(succ)

# ----------------------------
# DE 변이
# ----------------------------
def differential_index_mutation(x, swarm, F, lb, ub):
    ndim = len(x)
    r1, r2, r3 = np.random.choice(len(swarm), 3, replace=False)
    base = set(swarm[r1].astype(int))
    diff = set(swarm[r2].astype(int)) - set(swarm[r3].astype(int))
    trial = list(base.union(diff))
    if len(trial) > ndim:
        trial = np.random.choice(trial, ndim, replace=False)
    elif len(trial) < ndim:
        candidates = list(set(range(lb[0], ub[0] + 1)) - set(trial))
        trial += list(np.random.choice(candidates, ndim - len(trial), replace=False))
    return np.array(trial[:ndim])

# ----------------------------
# 초기 입자 생성
# ----------------------------
def initialize_discrete_particles(
    swarm_size, n_li, valid_sites, r_min, z_layer_tol, cell, max_per_layer,
    success_counter=None, initial_indices=None, base_structure=None, output_dir=None, variant_ratio=None
):
    particles = []
    valid_set = {tuple(np.round(site, 6)): i for i, site in enumerate(valid_sites)}

    budget = max(1, int(swarm_size * variant_ratio))  # 이번 iteration에 필요한 variant 수

    if initial_indices is not None:
        base_variants = generate_li_variants_budgeted(initial_indices, valid_sites, target_n_li=n_li, budget=budget)
        variant_output_path = os.path.join(output_dir, "variant_check_IN.cfg")
        with open(variant_output_path, "w") as f:
            for i, variant in enumerate(base_variants):
                selected_sites = np.array([valid_sites[idx] for idx in variant])
                atoms_with_li = add_lithium_atoms(base_structure, selected_sites, valid_sites, snap=False)
                write_mtp_cfg(atoms_with_li, f, swarmindex=i)
                f.write("\n")
        print(f"[INFO] POSCAR-based variants written to: {variant_output_path}", flush=True)
        for variant in base_variants:
            particles.append(np.array(variant, dtype=float))
        initial_count = len(base_variants)
    else:
        initial_count = 0

    max_trials = 10000
    trials = 0
    while len(particles) < swarm_size + initial_count and trials < max_trials:
        trials += 1
        try:
            selected_sites = select_li_sites_layerwise(
                valid_sites, n_li, r_min, max_per_layer,
                z_tol=z_layer_tol, cell=cell, success_counter=success_counter
            )
        except RuntimeError:
            print(f"[DEBUG] Trial {trials}: select_li_sites_layerwise failed", flush=True)
            continue

        indices = []
        for site in selected_sites:
            key = tuple(np.round(site, 6))
            if key in valid_set:
                indices.append(valid_set[key])
            else:
                found = False
                for k in valid_set:
                    if np.linalg.norm(np.array(k) - site) < 1e-3:
                        indices.append(valid_set[k]); found = True; break
                if not found:
                    print(f"[DEBUG] Trial {trials}: KeyError - no matching site found for {site}", flush=True)
                    indices = []
                    break
        if not indices:
            continue

        sorted_indices = sorted(indices)
        if any(np.array_equal(sorted_indices, sorted(p.astype(int))) for p in particles):
            print(f"[DEBUG] Trial {trials}: Duplicate structure detected", flush=True)
            continue
        particles.append(np.array(indices, dtype=float))

    if len(particles) == 0:
        print(f"[ERROR] No valid particles found using select_li_sites_layerwise after {max_trials} trials.", flush=True)
        raise RuntimeError("Initialization failed: No valid discrete particles generated.")

    if len(particles) < swarm_size:
        print(f"[WARNING] Only {len(particles)} particles generated via selector (target = {swarm_size}). Proceeding anyway.", flush=True)
    else:
        print(f"[INFO] Successfully initialized {len(particles)} discrete PSO particles using selector.", flush=True)
        print(f"       (POSCAR-based variants: {initial_count}, Randomly generated: {len(particles) - initial_count})", flush=True)

    return np.array(particles)

# ----------------------------
# 메인 PSO (ThreadPool로 병렬 평가, barrier 보장)
# ----------------------------

# 국소 변이#
def _local_mutation(particle, all_indices, g_best):
    """인덱스 기반 국소 변이: 스왑/치환/g_best로 당기기 중 1~2개 적용."""
    import numpy as np, random
    p = np.array(particle, dtype=int).copy()
    n = len(p)

    def op_swap(pp):
        if n < 2: return pp
        i, j = np.random.choice(n, 2, replace=False)
        pp[i], pp[j] = pp[j], pp[i]
        return pp

    def op_replace(pp):
        i = np.random.randint(n)
        candidates = list(set(all_indices) - set(pp))
        if candidates:
            pp[i] = np.random.choice(candidates)
        return pp

    def op_guided(pp):
        k = np.random.randint(n)
        pp[k] = int(g_best[k])
        return pp

    ops = [op_swap, op_replace, op_guided]
    for _ in range(np.random.randint(1, 3)):  # 1~2개 적용
        p = random.choice(ops)(p)
    return p

# 리히트 + 국소 변이 #
def custom_pso_with_hybrid_final(
    obj_func,              # 호환성 유지용(내부에서 직접 에너지 평가하므로 사용 안 함)
    valid_sites,
    r_min,
    z_layer_tol,
    cell,
    max_per_layer,
    n_li,
    swarmsize=20,
    max_iter=100,
    patience=25,
    save_progress_callback=None,
    return_history=False,
    save_structure_callback=True,
    log_file=None,
    output_dir=None,
    initial_indices=None,
    base_structure=None,
    F=None,
    use_de_prob=None,
    phase_ratio=None,
    explore_in_exploit_prob=0.15,
    debug=False,
    # variant 주입 관련
    variant_ratio=0.70,              # 한 iteration에서 스웜의 몇 %를 variant로 채울지
    variant_use_until=0.70,          # 전체 iter의 몇 %까지 variant 주입을 유지할지
    variant_refresh_every=1,         # 몇 iteration마다 새로 variant 풀 갱신할지
    variant_elite_share=0.50,        # 절반은 elite 기반 변형으로 생성
    min_symdiff=2,                   # variant 간 해밍거리(대칭차) 최소값 // 변형이 너무 비슷하면 min_symdiff를 3~5로 올려 유사 중복 억제.
    ##
    parallel_workers=None,   # ★ 동시 작업 수 옵션
    relax_workers=None,
    local_mutation_prob=0.30 # ★ 국소 변이 적용 확률
):
    import os, numpy as np, random
    from concurrent.futures import ThreadPoolExecutor, as_completed

    # 합리적 기본값 보정
    if F is None: F = 0.95
    if use_de_prob is None: use_de_prob = 0.9
    if phase_ratio is None: phase_ratio = 0.75

    if parallel_workers is None:
        parallel_workers = max(1, (os.cpu_count() or 2) // 2)
    if relax_workers is None:
        relax_workers = parallel_workers  # 별도 지정 없으면 동일 값 사용
        
    ndim = n_li
    lb = np.array([0] * ndim)
    ub = np.array([len(valid_sites) - 1] * ndim)
    all_indices = set(range(len(valid_sites)))

    log_fh = open(log_file, "w") if log_file else None
    def log(msg):
        if debug:
            print(msg)
        if log_fh:
            log_fh.write(msg + "\n")

    success_counter = {"count": 0}
    structure = base_structure  # 명시
    initial_particles = initialize_discrete_particles(
        swarm_size=swarmsize,
        n_li=n_li,
        valid_sites=valid_sites,
        r_min=r_min,
        z_layer_tol=z_layer_tol,
        cell=cell,
        max_per_layer=max_per_layer,
        success_counter=success_counter,
        initial_indices=initial_indices,
        base_structure=structure,
        output_dir=output_dir,
        variant_ratio=variant_ratio
    )

    x = np.array(initial_particles, dtype=object)
    actual_swarmsize = len(x)
    p_best = x.copy()
    p_best_vals = np.empty(actual_swarmsize, dtype=float)
    g_best = None
    g_best_val = float("inf")

    all_energies = []

    # ---- Iteration 0: 초기 파티클 병렬 평가 (barrier) ----
    jobs0 = [(x[i], 0, i) for i in range(actual_swarmsize)]
    def _eval(particle, iteration, swarmindex):
        return swarmindex, energy_function_discrete_index(
            np.array(particle, dtype=int), structure, n_li, r_min, z_layer_tol,
            pso_io_dir=output_dir, valid_sites=valid_sites, max_per_layer=max_per_layer,
            cell=cell, success_counter=success_counter, iteration=iteration, swarmindex=swarmindex
        ), particle

    with ThreadPoolExecutor(max_workers=parallel_workers) as ex:
        futures = [ex.submit(_eval, *job) for job in jobs0]
        for fut in as_completed(futures):
            idx, val, part = fut.result()
            p_best[idx] = part
            p_best_vals[idx] = val
            all_energies.append((0, idx, val))

    g_best = p_best[np.argmin(p_best_vals)]
    g_best_val = float(np.min(p_best_vals))
    best_val = g_best_val
    best_x = g_best.copy()
    no_improve_count = 0
    history = [(0, g_best_val)]

    convergence_data_path = os.path.join(output_dir, "convergence_data.txt")
    os.makedirs(output_dir, exist_ok=True)
    with open(convergence_data_path, "w") as f:
        f.write("# Iteration   SwarmIndex   Energy\n")
        for i in range(actual_swarmsize):
            f.write(f"0   {i}   {p_best_vals[i]:.6f}\n")

    # iter 0 relax: 성공시에만 저장 콜백
    succ0 = post_iteration_relaxation_parallel(0, output_dir=output_dir, valid_sites=valid_sites, relax_workers=relax_workers)
    succ0 = int(succ0 or 0)
    if succ0 > 0 and save_structure_callback:
        save_structure_callback(best_x, g_best_val, 0)

    seen_structures = set(tuple(sorted(p.astype(int))) for p in x)
    transition_iter = int(max_iter * phase_ratio) + 1

    # ---- 메인 루프 ----
    for iter in range(1, max_iter + 1):
        is_early = float(iter) / float(max_iter) < phase_ratio
        if iter == transition_iter:
            no_improve_count = 0
            log(f"[INFO] Exploitation phase begins at iter {iter}. Reset no_improve_count to 0.")

        # ---------- Variant injection (budgeted by ratio) ----------
        variant_candidates = []
        use_variant_now = (iter == 1) or (float(iter) / float(max_iter) <= variant_use_until)

        if use_variant_now:
            budget = max(1, int(swarmsize * variant_ratio))  # 이번 iteration에 필요한 variant 수
            variant_log_counts = {}  # ★ 전략별 생성 수 기록용

            # (옵션) remaining=0이어도 강제 섞기 켜기/비율
            force_mix_when_full = True
            force_mix_frac = 0.10  # budget의 10% 교체 시도

            # base_variants를 주기적으로 갱신
            if (iter == 1) or (iter % variant_refresh_every == 0) or (not base_variants_cache):
                try:
                    base_variants_cache = generate_li_variants_budgeted(
                        initial_indices, valid_sites, target_n_li=n_li, budget=budget
                    )
                except Exception:
                    base_variants_cache = []

            variant_candidates = []

            # 중복/유사 필터 포함
            def _add_variant(cand, source_label):
                key = tuple(sorted(map(int, cand)))
                if key in seen_structures:
                    return False
                if min_symdiff and min_symdiff > 0:
                    for v in variant_candidates:
                        if len(set(map(int, v)) ^ set(map(int, cand))) < min_symdiff:
                            return False
                variant_candidates.append(np.array(cand, dtype=float))
                seen_structures.add(key)
                variant_log_counts[source_label] = variant_log_counts.get(source_label, 0) + 1
                return True

            # 1) base variants로 먼저 채우기
            for v in base_variants_cache:
                if len(variant_candidates) >= budget:
                    break
                _add_variant(v, "base_variants_cache")

            # 2) elite 기반 소규모 변형으로 나머지 채우기
            remaining = budget - len(variant_candidates)      # <= 0일 수 있음
            if p_best_vals.size > 0:
                elite_k = max(1, swarmsize // 10)
                elite_order = np.argsort(p_best_vals)
                elites = [p_best[i].astype(int).copy() for i in elite_order[:elite_k]]

                def small_perturb(ind):
                    import numpy as np
                    ind = np.array(ind, dtype=int).copy()
                    k = np.random.randint(1, 3)  # 1~2자리 치환
                    pool = list(all_indices - set(ind))
                    for _ in range(k):
                        if not pool:
                            break
                        pos = np.random.randint(len(ind))
                        ind[pos] = np.random.choice(pool)
                    return np.array(sorted(ind), dtype=int)

                if remaining > 0:
                    target_elite = int(round(remaining * variant_elite_share))
                    trials = 0
                    # (a) elite 변형 먼저
                    while len(variant_candidates) < budget and trials < remaining * 20:
                        trials += 1
                        src = elites[np.random.randint(len(elites))]
                        cand = small_perturb(src)
                        _add_variant(cand, "elite_perturb")
                        if len(variant_candidates) >= budget or (
                            len(variant_candidates) >= budget - remaining + target_elite
                        ):
                            break

                    # (b) 아직 모자라면 base 변형/랜덤 치환으로 채우기
                    trials = 0
                    import numpy as np
                    while len(variant_candidates) < budget and trials < remaining * 30:
                        trials += 1
                        if base_variants_cache:
                            src = np.array(base_variants_cache[np.random.randint(len(base_variants_cache))], dtype=int)
                        else:
                            # fallback: g_best 살짝 섞기
                            src = np.array(g_best, dtype=int).copy()
                        cand = small_perturb(src)
                        _add_variant(cand, "random_perturb")

            # 3) ★ remaining==0 이어도 일부 강제 섞기(교체)
            if remaining == 0 and force_mix_when_full and p_best_vals.size > 0 and len(variant_candidates) > 0:
                import numpy as np
                force_k = max(1, int(round(budget * force_mix_frac)))
                force_k = min(force_k, len(variant_candidates))  # 안전 상한

                # 교체 대상 인덱스 뽑기
                replace_idx = np.random.choice(range(len(variant_candidates)), size=force_k, replace=False)

                replaced = 0
                trials = 0
                max_trials = force_k * 30

                while replaced < force_k and trials < max_trials:
                    trials += 1
                    ridx = replace_idx[replaced]  # 다음 교체 슬롯
                    src = elites[np.random.randint(len(elites))]
                    cand = small_perturb(src)

                    # 유사도/중복 체크
                    key = tuple(sorted(map(int, cand)))
                    if key in seen_structures:
                        continue
                    if min_symdiff and min_symdiff > 0:
                        too_similar = False
                        for v in variant_candidates:
                            if len(set(map(int, v)) ^ set(map(int, cand))) < min_symdiff:
                                too_similar = True
                                break
                        if too_similar:
                            continue

                    # 기존 것 제거 후 교체
                    old = variant_candidates[ridx]
                    old_key = tuple(sorted(map(int, old)))
                    if old_key in seen_structures:
                        seen_structures.remove(old_key)

                    variant_candidates[ridx] = np.array(cand, dtype=float)
                    seen_structures.add(key)
                    variant_log_counts["elite_forced_replace"] = variant_log_counts.get("elite_forced_replace", 0) + 1
                    replaced += 1

                if debug:
                    log(f"[VARIANT] Iter {iter}: forced mix applied (replaced {replaced}/{force_k})")

            # ★ 디버그 로그 출력
            if debug:
                log(f"[VARIANT] Iter {iter}: 총 {len(variant_candidates)}/{budget}개 생성")
                for src, cnt in variant_log_counts.items():
                    log(f"          - {src}: {cnt}개")
                 
        swarm_offset = len(variant_candidates)

        jobs = []
        # 1) variants 먼저 평가 큐에 넣기
        for vi, variant_particle in enumerate(variant_candidates):
            jobs.append((variant_particle, iter, vi))

        # 2) swarm 변형: 총량을 swarmsize로 제한
        quota = max(0, swarmsize - swarm_offset)  # 이번 iter에 swarm에서 만들 개수
        new_particles = []
        made = 0

        # swarm 변형
        new_particles = []
        for i in range(swarmsize):
            if made >= quota:           #  이번 iter에 swarm에서 만들 개수
                break
            for _ in range(20):
                particle = x[i].astype(int).copy()
                if is_early:
                    if np.random.rand() < use_de_prob:
                        particle = differential_index_mutation(x[i], x, F, lb, ub)
                    else:
                        j, k = np.random.choice(ndim, 2, replace=False)
                        particle[j], particle[k] = particle[k], particle[j]
                else:
                    if np.random.rand() < explore_in_exploit_prob:
                        particle = differential_index_mutation(x[i], x, F, lb, ub)
                    else:
                        shared = set(particle) & set(g_best)
                        if len(shared) > ndim:
                            shared = set(random.sample(shared, ndim))
                        remaining = list(all_indices - shared)
                        n_fill = ndim - len(shared)
                        if len(remaining) >= n_fill and n_fill > 0:
                            fill = np.random.choice(list(remaining), n_fill, replace=False)
                            particle = list(shared) + list(fill)
                        else:
                            for j in range(ndim):
                                if np.random.rand() < 0.3:
                                    particle[j] = g_best[j]
                        random.shuffle(particle)
                        particle = np.array(particle, dtype=int)

                    # ★ 국소 변이
                    if np.random.rand() < local_mutation_prob:
                        particle = _local_mutation(particle, all_indices, g_best)

                key = tuple(sorted(np.array(particle).astype(int)))
                if key not in seen_structures:
                    seen_structures.add(key)
                    new_particles.append((i, np.array(particle)))
                    break

        for i, particle in new_particles:
            full_index = swarm_offset + i
            jobs.append((particle, iter, full_index))

        # ---- 병렬 실행 ----
        iter_results = []  # (swarm_index_written, val, particle)
        with ThreadPoolExecutor(max_workers=parallel_workers) as ex:
            futures = [ex.submit(_eval, *job) for job in jobs]
            for fut in as_completed(futures):
                iter_results.append(fut.result())

        # ---- 결과 반영 ----
        with open(convergence_data_path, "a") as f:
            for swarm_idx_written, val, particle in sorted(iter_results, key=lambda t: t[0]):
                f.write(f"{iter}   {swarm_idx_written}   {val:.6f}\n")
                all_energies.append((iter, swarm_idx_written, val))

        # swarm의 p_best/x 갱신
        for i, particle in new_particles:
            full_index = swarm_offset + i
            found = [r for r in iter_results if r[0] == full_index]
            if not found:
                continue
            _, val, part = found[0]
            if val < p_best_vals[i]:
                p_best[i] = part
                p_best_vals[i] = val
            x[i] = part

        # g_best 갱신 (variants 포함)
        for swarm_idx_written, val, part in iter_results:
            if val < g_best_val:
                g_best_val = val
                g_best = part

        if g_best_val < best_val:
            log(f"[IMPROVED] Iter {iter:03d} | New Best Energy: {g_best_val:.6f} (Previous: {best_val:.6f})")
            best_val = g_best_val
            best_x = g_best.copy()
            no_improve_count = 0

            # 개선 시 relax → 성공시에만 저장 콜백
            succ = post_iteration_relaxation_parallel(iteration=iter, output_dir=output_dir, valid_sites=valid_sites, relax_workers=relax_workers)
            succ = int(succ or 0)
            if succ > 0 and save_structure_callback:
                save_structure_callback(best_x, g_best_val, iter)
        else:
            if not is_early:
                no_improve_count += 1
                log(f"[Iter {iter:03d}] Current Best: {g_best_val:.6f} | No improve count: {no_improve_count}/{patience}")
            else:
                log(f"[Iter {iter:03d}] Current Best: {g_best_val:.6f} (exploration phase – no patience check)")

        if return_history:
            history.append((iter, g_best_val))
        if save_progress_callback and iter % 5 == 0:
            save_progress_callback(history, all_energies, iter)

        # ---- 리히트(reheat): 개선 없이 patience 회수 지나면 재시작(엘리트 보존) ----
        if not is_early and no_improve_count >= patience:
            log(f"[REHEAT] No improvement for {no_improve_count} steps at iter {iter}. Reinitializing swarm with elites...")

            # 1) 엘리트 선별 (상위 10%)
            elite_k = max(1, swarmsize // 10)
            elite_order = np.argsort(p_best_vals)  # 낮은 에너지 우선
            elite_idx = elite_order[:elite_k]
            elites = [p_best[i].astype(int).copy() for i in elite_idx]

            # 2) 나머지 무작위 재초기화 (거리/층 제약 반영)
            valid_set = {tuple(np.round(site, 6)): i for i, site in enumerate(valid_sites)}
            new_particles = [np.array(e, dtype=float) for e in elites]

            tries = 0
            while len(new_particles) < swarmsize and tries < 10000:
                tries += 1
                try:
                    sel_sites = select_li_sites_layerwise(
                        valid_sites, n_li, r_min, max_per_layer,
                        z_tol=z_layer_tol, cell=cell, success_counter=None
                    )
                    # 좌표→인덱스
                    idxs = []
                    for s in sel_sites:
                        key = tuple(np.round(s, 6))
                        if key in valid_set:
                            idxs.append(valid_set[key])
                        else:
                            found = False
                            for k, vidx in valid_set.items():
                                if np.linalg.norm(np.array(k) - s) < 1e-3:
                                    idxs.append(vidx); found = True; break
                            if not found:
                                raise KeyError
                    cand = np.array(sorted(idxs), dtype=float)

                    # 중복 방지
                    if any(np.array_equal(cand, np.array(sorted(p.astype(int)))) for p in new_particles):
                        continue
                    new_particles.append(cand)
                except Exception:
                    continue

            if len(new_particles) < swarmsize:
                log(f"[REHEAT][WARN] Only reinitialized {len(new_particles)}/{swarmsize} particles.")

            # 3) 스웜 교체 및 상태 초기화
            x = np.array(new_particles, dtype=object)
            p_best = x.copy()
            p_best_vals = np.array([_eval(p, iter, i)[1] for i, p in enumerate(x)])
            g_best = p_best[np.argmin(p_best_vals)]
            g_best_val = np.min(p_best_vals)

            seen_structures = set(tuple(sorted(p.astype(int))) for p in x)
            best_val = min(best_val, g_best_val)
            no_improve_count = 0
            log(f"[REHEAT] New g_best after reheat: {g_best_val:.6f}. Swarm reset complete.")
            continue  # 다음 iteration

    if log_fh:
        log_fh.close()

    if return_history:
        return best_x, best_val, history
    else:
        return best_x, best_val

# ----------------------------
# 최상위 optimize + 아카이브
# ----------------------------
def optimize_li_with_min_distance(soc, POSCAR_file, stack_types, max_iter=100, swarmsize=20, patience=10, 
                                  mlp_model="pot.mtp", max_per_layer=1, use_input_c_layers=True, Ncell=2,
                                  struc_type=100, parallel_workers=None, relax_workers=None):
    def log(msg):
        print(msg, flush=True)

    start_time = time.time()
    success_counter = {"count": 0}
    base_structure = read(POSCAR_file, format='vasp')
    max_li = len([a for a in base_structure if a.symbol == "Li"])
    n_li = round(soc / 100 * max_li)
    print("max_li=",max_li, "n_li=",n_li, flush=True)
    if (n_li < 5):
        max_iter = int(max_iter * 0.5)
        swarmsize = int(swarmsize * 0.5)
        print("max_iter=",max_iter, "swarmsize=",swarmsize, flush=True)
    base_structure = base_structure[[atom.symbol != "Li" for atom in base_structure]]
    base_c_atoms = base_structure.copy()

    r_min, z_layer_tol = analyze_li_layer_distances(POSCAR_file, file_format='vasp')
    print("r_min=",r_min, flush=True)
    reference_valid_sites = None
    reference_bounds = None
    reference_rmin = None
    reference_ztol = None
    
    for stack_type in stack_types:
        if soc > 35 and stack_type != "AAAA":
            continue

        if use_input_c_layers:
            structure = generate_graphite_stack_from_input(base_c_atoms, stack_type, Ncell, struc_type)
        else:
            structure = generate_graphite_stack(stack_type, layers=4, spacing=3.35, Ncell=2, struc_type=100)
    
        output_dir = f"resultsOPT/SOC{soc}_Li{n_li}_{stack_type}"
        os.makedirs(output_dir, exist_ok=True)
    
        if stack_type == "AAAA":
            r_min, z_layer_tol = analyze_li_layer_distances(POSCAR_file, file_format='vasp')
            valid_sites = get_hex_center_sites(structure)
            if len(valid_sites) < n_li:
                print(f"[WARNING] Not enough valid Li sites ({len(valid_sites)}) for n_li={n_li}, skipping {stack_type}", flush=True)
                continue
            else:
                print(f"[DEBUG] valid_sites: {len(valid_sites)}", flush=True)
            np.random.shuffle(valid_sites)
            bounds = lithium_index_bounds(valid_sites, n_li)
            reference_valid_sites = valid_sites
            reference_bounds = bounds
            reference_rmin = r_min
            reference_ztol = z_layer_tol
        else:
            if reference_valid_sites is None:
                print(f"[Info] valid_sites generation using AAAA stacking", flush=True)
                structure_AAAA = generate_graphite_stack_from_input(base_c_atoms, "AAAA", Ncell, struc_type)
                r_min, z_layer_tol = analyze_li_layer_distances(POSCAR_file, file_format='vasp')
                valid_sites = get_hex_center_sites(structure_AAAA)
                if len(valid_sites) < n_li:
                    print(f"[WARNING] Not enough valid Li sites ({len(valid_sites)}) for n_li={n_li}, skipping {stack_type}", flush=True)
                    continue
                np.random.shuffle(valid_sites)
                bounds = lithium_index_bounds(valid_sites, n_li)
                structure = generate_graphite_stack_from_input(base_c_atoms, stack_type, Ncell, struc_type)
            else:
                valid_sites = reference_valid_sites
                bounds = reference_bounds
                r_min = reference_rmin
                z_layer_tol = reference_ztol
                print(f"[DEBUG] [Reuse] valid_sites: {len(valid_sites)}", flush=True)

        all_energies = []

        # 초기 Li indices from POSCAR
        atoms = read(POSCAR_file)
        li_positions = np.array([atom.position for atom in atoms if atom.symbol == "Li"])
        tree = cKDTree(valid_sites)
        initial_indices = []
        for pos in li_positions:
            dist, idx = tree.query(pos)
            if dist > 0.3:
                raise ValueError(f"Li at {pos} does not match any valid site (min dist = {dist:.3f})")
            initial_indices.append(idx)
        print(f"[INFO] Extracted {len(initial_indices)} Li site indices from {POSCAR_file} structure.", flush=True)

        x, val, history = custom_pso_with_hybrid_final(
            obj_func=None,
            swarmsize=swarmsize,
            max_iter=max_iter,
            patience=patience,
            log_file=os.path.join(output_dir, "log_pso.txt"),
            save_progress_callback=lambda history, all_energies, iter: save_convergence_plot(
                history, all_energies, output_dir, filename=f"convergence_iter{iter:03d}.png", data_filename=None
            ),
            save_structure_callback=lambda x_best, E_best, iter_stop: save_best_structure_from_convergence(
                cfg_dir=output_dir,
                convergence_file=os.path.join(output_dir, "convergence_data.txt"),
                output_dir=output_dir
            ),
            return_history=True,
            n_li=n_li,
            output_dir=output_dir,
            valid_sites=valid_sites,
            r_min=r_min,
            z_layer_tol=z_layer_tol,
            cell=structure.get_cell(),
            max_per_layer=max_per_layer,
            initial_indices=initial_indices,
            base_structure=structure,
            F=0.95,            # **큰 DE 보폭** → r2–r3 차이를 크게 반영하여 **멀리 떨어진 구조 생성**
            use_de_prob=0.9,   # DE 방식 선택 비율을 높여 **과감한 탐색 중심 전략** 사용
            phase_ratio=0.5,   # 0.6, 0.75
            explore_in_exploit_prob=0.15, # 0.05
            debug=True,
            parallel_workers=parallel_workers,  # ★ 동시 작업 수
            relax_workers=relax_workers
        )

        save_best_structure_from_convergence(
            cfg_dir=output_dir,
            convergence_file=os.path.join(output_dir, "convergence_data.txt"),
            output_dir=output_dir
        )

        import glob, tarfile, shutil
        def archive_cfg_files(output_dir):
            cfg_dir = os.path.join(output_dir, "CFG")
            os.makedirs(cfg_dir, exist_ok=True)
            for file in glob.glob(os.path.join(output_dir, "IN_*.cfg")) + glob.glob(os.path.join(output_dir, "OUT_???.cfg")):
                shutil.move(file, os.path.join(cfg_dir, os.path.basename(file)))
            tar_path = os.path.join(output_dir, "CFG.tar.gz")
            with tarfile.open(tar_path, "w:gz") as tar:
                tar.add(cfg_dir, arcname="CFG")
            shutil.rmtree(cfg_dir)
            print(f"[INFO] Archived CFG files to: {tar_path}")
        archive_cfg_files(output_dir)

        end_time = time.time()
        elapsed_time = end_time - start_time
        h, m, s = int(elapsed_time // 3600), int((elapsed_time % 3600) // 60), elapsed_time % 60
        log(f"[INFO] Total PSO run time: {h:02d}:{m:02d}:{s:05.2f} (hh:mm:ss)")

# ----------------------------
# CLI 실행부 (workers 옵션 지원)
# ----------------------------
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=None, help="동시 작업 수(N). 기본값 = CPU코어의 절반")
    parser.add_argument("--relax-workers", type=int, default=None, help="relax 병렬 작업 수(기본: workers)")
    args = parser.parse_args()

    cells = 2
    POSCAR_file = "POSCAR2x2_facing_P6MMM"
    stack_types = ["AAAA", "ABAA", "ABAB"]

    for soc in [ 50, 75, 25, 12]:                         # for 2x2 (4x4 = 16) 
        #       3*4  2*4  1*4  2
#    for soc in [67, 56, 50, 44, 33, 22, 11]:               # for 3x3 (9x4 = 36)
        #       6*4  5*4 18 4*4 3*4  2*4 1*4
#    for soc in [75,  50,  38,  31,  25,  19,  13,  6]:    # for 4x4 (16x4 = 64) 
        #      12*4  8*4  6*4  5*4  4*4  3*4  2*4  1*4
#    for soc in [75, 50,   33,   25,  20,  17,  14,  11,  8]:       # for 6x6 (36x4 = 144) 
        #      27*4  18*4  12*4  9*4  7*4  6*4  5*4  4*4  3*4
        optimize_li_with_min_distance(
            soc, POSCAR_file, stack_types,
            max_iter=200,
            swarmsize=200,
            patience=50,          
            mlp_model="pot.mtp",
            Ncell=cells,
            max_per_layer=cells*cells,
            use_input_c_layers=True,
            struc_type=100,             # SOC100 구조 -> 100 / SOC75구조 -> 75
            parallel_workers=args.workers,    # ★ 여기로 전달
            relax_workers=args.relax_workers
        )


"""
phase_ratio=0.5
swarmsize=200 (가능하면)
use_de_prob=0.9, F=0.95           0.8 // 0.85
explore_in_exploit_prob=0.15      0.05
patience=50 + 리히트 활성화

정체 해소(리히트)
patience ↓ : 60 → 25
조기수렴 시 재가열(restart): 상위 10% elite는 보존, 나머지 70%는 무작위 재초기화(+ 소수는 g_best 근방 가우시안 섭동).

 Variatn Types
    random=0.20,
    layerwise=0.15,     층마다 하나씩 
    alt_odd=0.15,       홀수층마다 하나씩
    fullwipe=0.15,     
    fullwipe_odd=0.15,
    facing=0.1,         모든 층 동일한 x,y 좌표 
    pair=0.1,           두층씩 동일한 x,y 좌표 // 다른 두층은 다른 x,y 좌표 

Variant 주입    
    variant_ratio=0.50       -> 0.70, # 한 iteration에서 스웜의 몇 %를 variant로 채울지
    variant_use_until=0.50,  -> 0.70  # 전체 iter의 몇 %까지 variant 주입을 유지할지
    variant_refresh_every=1,          # 몇 iteration마다 새로 variant 풀 갱신할지
    variant_elite_share=0.50,         # remaining의 절반은 elite 기반 변형으로 생성 => 0개 일 가능성 큼. 
                                   => remaining=0 일때 
                                      force_mix_frac = 0.10  # budget의 10% 교체 시도
    min_symdiff=2,                    # variant 간 해밍거리(대칭차) 최소값 // 변형이 너무 비슷하면 3~5로 올려 유사 중복 억제.
"""

"""
 python PSO_opt_parallel_var.py --workers 2 --relax-workers 2 >& LOG &

"""
