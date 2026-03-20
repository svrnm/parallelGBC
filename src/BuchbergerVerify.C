/**
 * Buchberger's criterion verification implementation.
 *
 * This file is part of parallelGBC.
 */
#include "../include/BuchbergerVerify.H"
#include <atomic>
#include <cstddef>
#include <ctime>
#include <iostream>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace parallelGBC {

static Polynomial spolynomial(const Polynomial& f, const Polynomial& g,
                              const CoeffField* field, const TOrdering* O) {
	Term lcm = f.LT().lcm(g.LT());
	Term term_f = lcm.div(f.LT());
	Term term_g = lcm.div(g.LT());

	coeffType inv_lc_f = field->inv(f.LC());
	coeffType inv_lc_g = field->inv(g.LC());

	Polynomial s1 = f.mul(term_f);
	s1.mulBy(inv_lc_f, field);

	Polynomial s2 = g.mul(term_g);
	s2.mulBy(inv_lc_g, field);

	return Polynomial::sub(s1, s2, field, O);
}

static Polynomial reduce(const Polynomial& p, const std::vector<Polynomial>& G,
                         const CoeffField* field, const TOrdering* O) {
	Polynomial r = p;
	r.bringIn(field, false);
	r.order(O);

	while(!r.isZero()) {
		bool reduced = false;
		for(size_t i = 0; i < r.size() && !reduced; i++) {
			Monomial mi = r[i];
			for(size_t k = 0; k < G.size(); k++) {
				if(G[k].isZero()) {
					continue;
				}
				if(mi.second.isDivisibleBy(G[k].LT())) {
					Term m = mi.second.div(G[k].LT());
					coeffType c = field->div(mi.first, G[k].LC());
					Polynomial g_scaled = G[k].mul(m);
					g_scaled.mulBy(c, field);
					r = Polynomial::sub(r, g_scaled, field, O);
					r.bringIn(field, false);
					r.order(O);
					reduced = true;
					break;
				}
			}
		}
		if(!reduced) {
			break;
		}
	}
	return r;
}

/** Map linear index 0 .. n*(n-1)/2 - 1 to pairs (i,j) with i < j in row-major order. */
static void pairFromLinearIndex(size_t n, size_t idx, size_t& i, size_t& j) {
	for(i = 0; i < n; i++) {
		size_t rowLen = n - 1 - i;
		if(idx < rowLen) {
			j = i + 1 + idx;
			return;
		}
		idx -= rowLen;
	}
}

bool verifyBuchbergerCriterion(const std::vector<Polynomial>& G,
                               const CoeffField* field,
                               const TOrdering* O,
                               bool showProgress,
                               size_t progressStepPercent) {
	std::vector<Polynomial> G_ordered;
	G_ordered.reserve(G.size());
	for(size_t k = 0; k < G.size(); k++) {
		Polynomial g = G[k];
		g.bringIn(field, false);
		g.order(O);
		G_ordered.push_back(g);
	}
	const size_t n = G_ordered.size();
	const size_t totalPairs = (n > 1) ? ((n * (n - 1)) / 2) : 0;

	if(showProgress && totalPairs == 0) {
		std::cerr << "Buchberger verify: 0/0 pairs (trivial)\n";
		return true;
	}

#ifdef _OPENMP
	const int maxThreads = omp_get_max_threads();
	const bool useParallel = (totalPairs > 8 && maxThreads > 1);
	if(showProgress && useParallel) {
		std::cerr << "Buchberger verify: " << totalPairs << " S-pairs, OpenMP max_threads="
		          << maxThreads << "\n";
	}
	if(useParallel) {
		std::atomic<bool> ok(true);
		#pragma omp parallel for schedule(dynamic, 4)
		for(ptrdiff_t p = 0; p < static_cast<ptrdiff_t>(totalPairs); p++) {
			if(!ok.load(std::memory_order_relaxed)) {
				continue;
			}
			size_t i = 0;
			size_t j = 0;
			pairFromLinearIndex(n, static_cast<size_t>(p), i, j);
			if(G_ordered[i].isZero() || G_ordered[j].isZero()) {
				continue;
			}
			Polynomial s = spolynomial(G_ordered[i], G_ordered[j], field, O);
			if(!s.isZero()) {
				Polynomial r = reduce(s, G_ordered, field, O);
				if(!r.isZero()) {
					ok.store(false, std::memory_order_relaxed);
					if(showProgress) {
						#pragma omp critical(buchberger_verify_fail)
						std::cerr << "Buchberger verify: FAILED at pair (" << i << "," << j << ")\n";
					}
				}
			}
		}
		if(showProgress && ok.load()) {
			std::cerr << "Buchberger verify: finished " << totalPairs << " pairs\n";
		}
		return ok.load();
	}
#endif

	// Serial path (small instance, no OpenMP, or single thread)
	size_t checkedPairs = 0;
	size_t nextProgressPercent = (progressStepPercent == 0) ? 100 : progressStepPercent;
	time_t lastHeartbeat = std::time(nullptr);
	if(showProgress) {
		std::cerr << "Buchberger verify: 0/" << totalPairs << " pairs (0%)\n";
	}
	for(size_t i = 0; i < G_ordered.size(); i++) {
		for(size_t j = i + 1; j < G_ordered.size(); j++) {
			if(showProgress) {
				time_t now = std::time(nullptr);
				if(now - lastHeartbeat >= 5) {
					std::cerr << "Buchberger verify: working on pair (" << i << "," << j
					          << "), completed " << checkedPairs << "/" << totalPairs << " pairs\n";
					lastHeartbeat = now;
				}
			}
			checkedPairs++;
			if(G_ordered[i].isZero() || G_ordered[j].isZero()) {
				continue;
			}
			Polynomial s = spolynomial(G_ordered[i], G_ordered[j], field, O);
			if(!s.isZero()) {
				Polynomial r = reduce(s, G_ordered, field, O);
				if(!r.isZero()) {
					if(showProgress) {
						std::cerr << "Buchberger verify: FAILED at pair (" << i << "," << j
						          << "), checked " << checkedPairs << "/" << totalPairs << " pairs\n";
					}
					return false;
				}
			}
			if(showProgress && totalPairs > 0) {
				size_t percent = (checkedPairs * 100) / totalPairs;
				if(percent >= nextProgressPercent || checkedPairs == totalPairs) {
					std::cerr << "Buchberger verify: " << checkedPairs << "/" << totalPairs
					          << " pairs (" << percent << "%)\n";
					while(nextProgressPercent <= percent) {
						nextProgressPercent += (progressStepPercent == 0) ? 100 : progressStepPercent;
					}
				}
			}
		}
	}
	return true;
}

}
