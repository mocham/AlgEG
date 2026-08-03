"""Self-contained elliptic Weyl orbit computation for all exceptional types.

Generates the data that ``twisted_cycles.py`` needs directly from root data.
Covers the elliptic Levi/Weyl configuration enumeration for
G2, F4, E6, E7, E8.
"""

from itertools import combinations, permutations, product

from sage.all import RootSystem, WeylGroup, identity_matrix


def _get_standard_cartan_matrix(letter, rank):
    I = list(range(1, rank + 1))
    C = {i: {j: (2 if i == j else 0) for j in I} for i in I}

    def add_edge(i, j, val_ij, val_ji):
        C[i][j] = val_ij
        C[j][i] = val_ji

    if letter == "A":
        for i in range(1, rank):
            add_edge(i, i + 1, -1, -1)
    elif letter == "B":
        for i in range(1, rank - 1):
            add_edge(i, i + 1, -1, -1)
        if rank >= 2:
            add_edge(rank - 1, rank, -2, -1)
    elif letter == "C":
        for i in range(1, rank - 1):
            add_edge(i, i + 1, -1, -1)
        if rank >= 2:
            add_edge(rank - 1, rank, -1, -2)
    elif letter == "D":
        for i in range(1, rank - 1):
            add_edge(i, i + 1, -1, -1)
        if rank >= 3:
            add_edge(rank - 2, rank, -1, -1)
    elif letter == "E":
        add_edge(1, 3, -1, -1)
        add_edge(2, 4, -1, -1)
        add_edge(3, 4, -1, -1)
        if rank >= 5:
            add_edge(4, 5, -1, -1)
        if rank >= 6:
            add_edge(5, 6, -1, -1)
        if rank >= 7:
            add_edge(6, 7, -1, -1)
        if rank >= 8:
            add_edge(7, 8, -1, -1)
    elif letter == "F":
        add_edge(1, 2, -1, -1)
        add_edge(2, 3, -2, -1)
        add_edge(3, 4, -1, -1)
    elif letter == "G":
        add_edge(1, 2, -1, -3)
    return C, I


ELLIPTIC_CACHE = {}


def _get_standard_elliptic_words(type_tuple):
    if type_tuple not in ELLIPTIC_CACHE:
        W = WeylGroup(type_tuple)
        rank = type_tuple[1]
        words = []
        for c in W.conjugacy_classes():
            w = c.representative()
            mat = w.matrix()
            if (mat - identity_matrix(mat.nrows())).rank() == rank:
                words.append(w.reduced_word())
        ELLIPTIC_CACHE[type_tuple] = words
    return ELLIPTIC_CACHE[type_tuple]


def _get_submatrix(C_amb, indices):
    return {u: {v: C_amb[u][v] for v in indices} for u in indices}


def _identify_and_embed(comp_indices, C_amb):
    k = len(comp_indices)
    M_C = _get_submatrix(C_amb, comp_indices)
    if k == 1:
        candidates = [("A", 1)]
    elif k == 2:
        candidates = [("A", 2), ("B", 2), ("C", 2), ("G", 2)]
    elif k == 3:
        candidates = [("A", 3), ("B", 3), ("C", 3)]
    elif k == 4:
        candidates = [("A", 4), ("B", 4), ("C", 4), ("D", 4), ("F", 4)]
    else:
        candidates = [("A", k), ("B", k), ("C", k), ("D", k)]
        if k in (6, 7, 8):
            candidates.append(("E", k))
    for letter, rank in candidates:
        C_std, I_std = _get_standard_cartan_matrix(letter, rank)
        for p in permutations(comp_indices):
            match = True
            for i_std, i_amb in zip(I_std, p):
                for j_std, j_amb in zip(I_std, p):
                    if C_std[i_std][j_std] != M_C[i_amb][j_amb]:
                        match = False
                        break
                if not match:
                    break
            if match:
                return (letter, rank), {i_std: i_amb for i_std, i_amb in zip(I_std, p)}
    raise ValueError(f"Could not identify Cartan submatrix for {comp_indices}")


def get_elliptic_levi_representatives(cartan_type):
    amb_letter, amb_rank = cartan_type
    C_amb, I_amb = _get_standard_cartan_matrix(amb_letter, amb_rank)
    result = {}
    for r in range(len(I_amb)):
        for J in combinations(I_amb, r):
            if not J:
                continue
            J_list = list(J)
            adj = {u: [] for u in J_list}
            for u in J_list:
                for v in J_list:
                    if u != v and C_amb[u][v] != 0:
                        adj[u].append(v)
            visited = set()
            components = []
            for u in J_list:
                if u not in visited:
                    comp = []
                    q = [u]
                    visited.add(u)
                    while q:
                        curr = q.pop(0)
                        comp.append(curr)
                        for neighbor in adj[curr]:
                            if neighbor not in visited:
                                visited.add(neighbor)
                                q.append(neighbor)
                    components.append(comp)
            comp_words = []
            for C in components:
                T_tuple, mapping = _identify_and_embed(C, C_amb)
                std_words = _get_standard_elliptic_words(T_tuple)
                mapped = [tuple(mapping[x] for x in w) for w in std_words]
                comp_words.append(mapped)
            combined_reps = []
            for word_tuple in product(*comp_words):
                combined = []
                for w in word_tuple:
                    combined.extend(w)
                combined_reps.append(tuple(int(x) for x in combined))
            result["".join(str(x) for x in J)] = combined_reps
    return result


def compute_poset_w_orbits(cartan_type, reduced_word):
    R = RootSystem(cartan_type)
    L = R.root_lattice()
    W = L.weyl_group()
    w = W.from_reduced_word(reduced_word)
    str_to_lattice = {str(root): root for root in L.roots()}
    n = len(L.roots()[0].to_vector())
    poset_roots = []
    for alpha in R.root_poset():
        lattice_alpha = str_to_lattice[str(alpha)]
        poset_roots.append(lattice_alpha)
    seen = set()
    orbits = []
    for alpha in poset_roots:
        if alpha not in seen:
            orbit = []
            current = alpha
            good = True
            while current not in seen:
                vec_str = "".join(str(c) for c in current.to_vector())
                if len(vec_str) != n:
                    good = False
                orbit.append(vec_str)
                seen.add(current)
                current = w.action(current)
            if good:
                orbits.append(orbit)
    return orbits


EXCEPTIONAL_TYPES = [("G", 2), ("F", 4), ("E", 6), ("E", 7), ("E", 8)]


def compute_all_elliptic_weyl_orbits():
    ells = {}
    for cts in EXCEPTIONAL_TYPES:
        key = f"{cts[0]}{cts[1]}"
        ells[key] = get_elliptic_levi_representatives(cts)
    orbits = {}
    cnt = 0
    for ct, items in ells.items():
        orbits[ct] = {}
        for levi, rws in items.items():
            orbits[ct][levi] = {}
            for rw in rws:
                rwstr = "".join(str(x) for x in rw)
                cnt += 1
                orbits[ct][levi][rwstr] = compute_poset_w_orbits(
                    (ct[0], int(ct[1])), rw)
    return orbits


def compute_nonobstructing_orbit_pairs(levi, orbits, Phi):
    root_diff = lambda r1, r2: "".join(
        str(int(d1) - int(d2)) for d1, d2 in zip(r1, r2))
    root_supp = lambda root: set(i + 1 for i, d in enumerate(root) if d != '0')
    levi_supp = set(int(i) for i in levi)
    obstructing_matrix = {}
    for orbit_id1, orbit1 in enumerate(orbits):
        for orbit_id2, orbit2 in enumerate(orbits):
            obstructing = True
            for r1 in orbit1:
                for r2 in orbit2:
                    r_diff = root_diff(r1, r2)
                    if r_diff not in Phi:
                        continue
                    r_diff_supp = root_supp(r_diff)
                    if len(r_diff_supp) > 0 and r_diff_supp.issubset(levi_supp):
                        obstructing = False
                        break
                if not obstructing:
                    break
            obstructing_matrix[(orbit_id1, orbit_id2)] = obstructing
    return [pair for pair, check in obstructing_matrix.items() if not check]


def compute_orbit_sums(levi, orbits):
    total_ht = lambda root: sum(int(digit) for digit in root)
    levi_ht = lambda root: sum(int(root[int(d) - 1]) for d in levi)
    root_sum = lambda r1, r2: "".join(
        str(int(d1) + int(d2)) for d1, d2 in zip(r1, r2))
    orbit_dict = {}
    simple_orbits = []
    sum_orbits = {}
    for orbit_id, roots in enumerate(orbits):
        for root in roots:
            orbit_dict[root] = orbit_id
        if total_ht(roots[0]) - levi_ht(roots[0]) == 1:
            simple_orbits.append(orbit_id)
    for orbit_id1, orbit1 in enumerate(orbits):
        for orbit_id2, orbit2 in enumerate(orbits):
            sum_orbit_ids = set()
            for r1 in orbit1:
                for r2 in orbit2:
                    r_sum = root_sum(r1, r2)
                    if r_sum in orbit_dict:
                        sum_orbit_ids.add(orbit_dict[r_sum])
            if sum_orbit_ids:
                sum_orbits[f"{orbit_id1}+{orbit_id2}"] = sorted(sum_orbit_ids)
    return simple_orbits, sum_orbits


def get_less_than(orbits):
    orbit_dict = {}
    roots = []
    for orb_id, orbit in enumerate(orbits):
        for root in orbit:
            orbit_dict[root] = orb_id
        roots.extend(orbit)
    less_thans = set()
    for r1 in roots:
        for r2 in roots:
            if all(d1 <= d2 for d1, d2 in zip(r1, r2)):
                less_thans.add((orbit_dict[r1], orbit_dict[r2]))
    return less_thans


def _directed_heights(tree):
    """Intrinsic heights for a DAG whose vertices need not be root tuples."""
    parents = {vertex: set() for vertex in tree}
    for source, targets in tree.items():
        for target in targets:
            parents.setdefault(target, set()).add(source)
    heights = {}

    def visit(vertex):
        if vertex not in heights:
            heights[vertex] = 1 + max((visit(parent)
                                       for parent in parents.get(vertex, ())),
                                      default=0)
        return heights[vertex]

    for vertex in sorted(parents):
        visit(vertex)
    return heights


def compute_non_borel_strata(orbits):
    from roots import positive_roots_set, root_less_than, downward_closed_subsets

    records = {}
    for ct, items in orbits.items():
        cartan_type = (ct[0], int(ct[1:]))
        Phi = positive_roots_set(cartan_type)
        records[ct] = {}
        for levi, rw_items in items.items():
            records[ct][levi] = {}
            for rw, orbit_list in rw_items.items():
                simples, sum_dict = compute_orbit_sums(levi, orbit_list)
                records[ct][levi][rw] = {
                    'orbits': orbit_list,
                    'simple_orbits': simples,
                    'orbit_sum': sum_dict,
                    'nonobstructing_orbit_pairs': compute_nonobstructing_orbit_pairs(
                        levi, orbit_list, {''.join(str(c) for c in r) for r in Phi}),
                }

    brecords = {}
    for ctx_key, ctx_val in records.items():
        cartan_type = (ctx_key[0], int(ctx_key[1:]))
        Phi = positive_roots_set(cartan_type)
        brecords[ctx_key] = {}
        for levi, rw_items in ctx_val.items():
            brecords[ctx_key][levi] = {}
            for rw, rec in rw_items.items():
                orbit_dict = {}
                orbit_str = {}
                tree = {}
                orbital_nonobstructing_pairs = set(
                    tuple(sorted(pair)) for pair in rec['nonobstructing_orbit_pairs'])
                for orb_id, orbit in enumerate(rec['orbits']):
                    st = "-".join(orbit)
                    orbit_dict[st] = orb_id
                    orbit_str[orb_id] = st
                for orb_id in range(len(rec['orbits'])):
                    tree[orbit_str[orb_id]] = []
                    for simple_id in rec['simple_orbits']:
                        sum_key = f"{simple_id}+{orb_id}"
                        if sum_key in rec['orbit_sum']:
                            for sum_id in rec['orbit_sum'][sum_key]:
                                tree[orbit_str[orb_id]].append(orbit_str[sum_id])
                nonobstructing_roots = {
                    v for k, v in orbit_str.items()
                    if (k, k) in orbital_nonobstructing_pairs}
                less_thans = get_less_than(rec['orbits'])

                layers = {}
                ht = _directed_heights(tree) if tree else {}
                for src, targets in tree.items():
                    h = ht.get(src, 0)
                    if h not in layers:
                        layers[h] = {}
                    layers[h][src] = tuple(
                        t for t in targets
                        if ht.get(t, -1) == h + 1 and t not in nonobstructing_roots)
                layers = {h: l for h, l in layers.items() if l}

                strata_info = None
                if ctx_key == 'F4':
                    from roots import downward_closed_subsets
                    Ks = downward_closed_subsets(
                        list(range(len(rec['orbits']))),
                        lambda a, b: (a, b) in less_thans)
                    Psi2s_list = []
                    orbit_indices = list(range(len(rec['orbits'])))
                    for basis in range(1, 1 << len(orbit_indices)):
                        basis_set = [i for i in orbit_indices if basis & (1 << i)]
                        if not basis_set:
                            continue
                        nonobs = False
                        for i in basis_set:
                            for j in basis_set:
                                if (i, j) in orbital_nonobstructing_pairs:
                                    nonobs = True
                            if nonobs:
                                break
                        if nonobs:
                            continue
                        Psi2s_list.append(basis_set)
                    strata_info = {'Ks': [[orbit_str[oid] for oid in K]
                                          for K in Ks if K],
                                   'Psi2s': [[orbit_str[oid] for oid in p]
                                              for p in Psi2s_list]}

                brecords[ctx_key][levi][rw] = {
                    'tree': tree,
                    'nonobstructing_roots': sorted(nonobstructing_roots),
                    'orbit_dict': orbit_dict,
                    'orbit_str': orbit_str,
                    'nonobstructing_orbit_pairs': sorted(
                        orbital_nonobstructing_pairs),
                    'less_thans': sorted(less_thans),
                    'bipartites': {str(h): v for h, v in layers.items()},
                    'strata': strata_info,
                }
    return brecords
