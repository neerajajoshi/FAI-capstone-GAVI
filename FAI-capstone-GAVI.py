import random
import time


class CSP:
    def __init__(self, n, L, T, R, solution):
        self.n        = n
        self.L        = L
        self.T        = T
        self.R        = R
        self.ground_truth = solution

    @property
    def m(self):
        return len(self.T)

    @property
    def density(self):
        return self.m / self.n

    @property
    def R_mean(self):
        return sum(len(r) for r in self.R) / self.m


def generate_csp(n, m, L_size, R_mean_target, seed=None):
    """Generate a random binary CSP guaranteed to have at least one solution."""
    rng = random.Random(seed)
    L = list(range(L_size))

    solution = [rng.choice(L) for _ in range(n)]

    all_pairs = [(p, q) for p in range(n) for q in range(p+1, n)]
    T = rng.sample(all_pairs, m)

    all_val_pairs = [(a, b) for a in L for b in L]

    R = []
    for (p, q) in T:
        sol_pair = (solution[p], solution[q])

        k = max(1, round(rng.gauss(R_mean_target, 0.5)))
        k = min(k, len(all_val_pairs))

        extras = [vp for vp in all_val_pairs if vp != sol_pair]
        rng.shuffle(extras)
        Rj = set([sol_pair] + extras[:k - 1])
        R.append(Rj)

    return CSP(n, L, T, R, solution)


def num_conflicts(a, csp):
    return sum(1 for j,(p,q) in enumerate(csp.T) if (a[p],a[q]) not in csp.R[j])

def fitness(a, csp):
    return 1.0 - num_conflicts(a, csp) / csp.m

def is_solution(a, csp):
    return num_conflicts(a, csp) == 0

def argmax(lst):
    return max(range(len(lst)), key=lambda i: lst[i])

def listmean(lst):
    return sum(lst)/len(lst) if lst else None

def roulette(pop, fits, rng):
    weights = [f*f for f in fits]
    total = sum(weights)
    if total == 0:
        return rng.choice(pop)[:]
    r = rng.uniform(0, total)
    cum = 0.0
    for ind, w in zip(pop, weights):
        cum += w
        if cum >= r:
            return ind[:]
    return pop[-1][:]

def uniform_cx(p1, p2, rng):
    return [p1[i] if rng.random() < 0.5 else p2[i] for i in range(len(p1))]


def ihc(csp, max_hc_steps=500, max_restarts=200, rng=None):
    if rng is None: rng = random.Random()

    for _ in range(max_restarts):
        a = [rng.choice(csp.L) for _ in range(csp.n)]

        for _ in range(max_hc_steps):
            if is_solution(a, csp):
                return a, True

            cvars = set()
            for j,(p,q) in enumerate(csp.T):
                if (a[p],a[q]) not in csp.R[j]:
                    cvars.add(p)
                    cvars.add(q)
            if not cvars:
                return a, True

            var = rng.choice(list(cvars))

            bv, bc = a[var], num_conflicts(a, csp)
            vals = csp.L[:]
            rng.shuffle(vals)
            for val in vals:
                a[var] = val
                c = num_conflicts(a, csp)
                if c < bc:
                    bc, bv = c, val
            a[var] = bv

    return None, False


def ga(csp, pop_size=500, max_gen=200, mutation_rate=0.02, rng=None):
    if rng is None: rng = random.Random()

    pop = [[rng.choice(csp.L) for _ in range(csp.n)] for _ in range(pop_size)]

    for gen in range(max_gen):
        fits = [fitness(ind, csp) for ind in pop]

        bi = argmax(fits)
        if fits[bi] == 1.0:
            return pop[bi], True

        elite = pop[bi][:]
        new_pop = [elite]

        while len(new_pop) < pop_size:
            child = uniform_cx(roulette(pop, fits, rng),
                               roulette(pop, fits, rng), rng)

            for i in range(csp.n):
                if rng.random() < mutation_rate:
                    child[i] = rng.choice(csp.L)

            new_pop.append(child)

        pop = new_pop

    fits = [fitness(ind, csp) for ind in pop]
    bi = argmax(fits)
    return (pop[bi], True) if fits[bi] == 1.0 else (None, False)


def gavi(csp, pop_size=500, max_gen=200, infection_rate=0.4, rng=None):
    if rng is None: rng = random.Random()

    virus_pool = [((p,q),(vp,vq))
                  for j,(p,q) in enumerate(csp.T)
                  for (vp,vq) in csp.R[j]]

    if not virus_pool:
        return ga(csp, pop_size, max_gen, rng=rng)

    random.shuffle(virus_pool)

    infectivity = [1.0] * len(virus_pool)

    pop = [[rng.choice(csp.L) for _ in range(csp.n)] for _ in range(pop_size)]

    n_infect = max(1, int(pop_size * infection_rate))

    for gen in range(max_gen):
        fits = [fitness(ind, csp) for ind in pop]

        bi = argmax(fits)
        if fits[bi] == 1.0:
            return pop[bi], True

        elite = pop[bi][:]
        new_pop = [elite]

        while len(new_pop) < pop_size:
            child = uniform_cx(roulette(pop, fits, rng),
                               roulette(pop, fits, rng), rng)
            new_pop.append(child)
        pop = new_pop

        vt = sum(infectivity)

        for k in rng.sample(range(pop_size), min(n_infect, pop_size)):

            r = rng.uniform(0, vt)
            cum = 0.0; vi = 0
            for idx, inf in enumerate(infectivity):
                cum += inf
                if cum >= r:
                    vi = idx
                    break

            (pp,qq),(vp,vq) = virus_pool[vi]
            old_f = fitness(pop[k], csp)

            pop[k][pp] = vp
            pop[k][qq] = vq

            new_f = fitness(pop[k], csp)

            delta = 0.1 if new_f > old_f else -0.1
            infectivity[vi] = max(0.01, infectivity[vi] + delta)

            vt = sum(infectivity)

    fits = [fitness(ind, csp) for ind in pop]
    bi = argmax(fits)
    return (pop[bi], True) if fits[bi] == 1.0 else (None, False)


def run_experiment(n_csps, n, m, L_size, R_mean, algorithms, seed_base=42):
    """
    Run each algorithm on n_csps randomly generated CSP instances.
    Returns aggregated success rate and mean solution time per algorithm.
    """
    results = {name: {"solved": 0, "times": []} for name in algorithms}

    for i in range(n_csps):
        csp = generate_csp(n, m, L_size, R_mean, seed=seed_base + i)

        for name, algo in algorithms.items():
            rng = random.Random(seed_base + i + 99999)

            t0 = time.perf_counter()
            sol, ok = algo(csp, rng=rng)
            elapsed = time.perf_counter() - t0

            if ok and sol is not None and is_solution(sol, csp):
                results[name]["solved"] += 1
                results[name]["times"].append(elapsed)

    for r in results.values():
        r["success_pct"] = 100.0 * r["solved"] / n_csps
        r["mean_time"]   = listmean(r["times"])

    return results


if __name__ == "__main__":

    CSP_COUNTS = [10, 30, 50]  # Run experiments for all three sizes

    algos = {
        "IHC":  lambda csp, rng: ihc(csp,  max_hc_steps=500, max_restarts=200, rng=rng),
        "GA":   lambda csp, rng: ga(csp,   pop_size=500, max_gen=200, rng=rng),
        "GAVI": lambda csp, rng: gavi(csp, pop_size=500, max_gen=200,
                                      infection_rate=0.4, rng=rng),
    }

    density_configs = [
        (49,  0.98),
        (100, 2.0),
        (150, 3.0),
        (200, 4.0),
        (250, 5.0),
    ]
    R_means = [2, 4, 6]

    # Store all results for comparative analysis: 
    all_results = {}

    for n_csps in CSP_COUNTS:
        all_results[n_csps] = {}

        for R_mean in R_means:
            all_results[n_csps][R_mean] = {}

            print(f"\n{'='*70}")
            print(f"  N_CSPS = {n_csps}  |  R_mean = {R_mean}  |  (n=50, |L|=4)")
            print(f"{'='*70}")
            print(f"  {'d':>5}  {'Algo':>6}  {'Success%':>10}  {'MeanTime(s)':>13}")
            print(f"  {'-'*45}")

            for m_val, d in density_configs:
                all_results[n_csps][R_mean][d] = {}

                res = run_experiment(n_csps, n=50, m=m_val, L_size=4,
                                     R_mean=R_mean, algorithms=algos, seed_base=0)

                for name, r in res.items():
                    all_results[n_csps][R_mean][d][name] = r
                    t = f"{r['mean_time']:.4f}" if r['mean_time'] is not None else "    N/A"
                    print(f"  {d:>5}  {name:>6}  {r['success_pct']:>10.1f}  {t:>13}")
                print()

    # Comparative Analysis 
    print(f"\n{'='*70}")
    print("  COMPARATIVE ANALYSIS: Effect of N_CSPS on Success Rate (%)")
    print(f"{'='*70}")

    for R_mean in R_means:
        print(f"\n  R_mean = {R_mean}")
        print(f"  {'d':>5}  {'Algo':>6}  {'N=10':>8}  {'N=30':>8}  {'N=50':>8}")
        print(f"  {'-'*45}")

        for m_val, d in density_configs:
            for name in algos:
                row = f"  {d:>5}  {name:>6}"
                for n_csps in CSP_COUNTS:
                    pct = all_results[n_csps][R_mean][d][name]["success_pct"]
                    row += f"  {pct:>7.1f}%"
                print(row)
            print()

    print(f"\n{'='*70}")
    print("  COMPARATIVE ANALYSIS: Effect of N_CSPS on Mean Solve Time (s)")
    print(f"{'='*70}")

    for R_mean in R_means:
        print(f"\n  R_mean = {R_mean}")
        print(f"  {'d':>5}  {'Algo':>6}  {'N=10':>10}  {'N=30':>10}  {'N=50':>10}")
        print(f"  {'-'*52}")

        for m_val, d in density_configs:
            for name in algos:
                row = f"  {d:>5}  {name:>6}"
                for n_csps in CSP_COUNTS:
                    mt = all_results[n_csps][R_mean][d][name]["mean_time"]
                    row += f"  {mt:.4f}s" if mt is not None else f"  {'N/A':>9}"
                print(row)
            print()

    print("Done!")