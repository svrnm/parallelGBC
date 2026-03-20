/**
 * Buchberger's criterion verification for Gröbner bases.
 *
 * A set G is a Gröbner basis iff for all pairs (g_i, g_j), the S-polynomial
 * S(g_i, g_j) reduces to 0 modulo G.
 *
 * This file is part of parallelGBC.
 */
#ifndef BUCHBERGER_VERIFY_H
#define BUCHBERGER_VERIFY_H
#include "Polynomial.H"
#include "TOrdering.H"

namespace parallelGBC {

/**
 * Verify that the given set of polynomials forms a Gröbner basis using
 * Buchberger's criterion. Returns true if all S-polynomials reduce to zero.
 * S-pair checks are parallelized with OpenMP when available (read-only G).
 */
bool verifyBuchbergerCriterion(const std::vector<Polynomial>& G,
                               const CoeffField* field,
                               const TOrdering* O,
                               bool showProgress = false,
                               size_t progressStepPercent = 5);

}
#endif
