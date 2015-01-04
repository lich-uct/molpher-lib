/*
 Copyright (c) 2012 Peter Szepe

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "FingerprintStrategy.h"

// include strategies ..
#include "AtomPairsFngpr.hpp"
#include "MorganFngpr.hpp"
#include "TopolLayeredFngpr1.hpp"
#include "TopolLayeredFngpr2.hpp"
#include "TopolSingleFngpr.hpp"
#include "TopolTorsFngpr.hpp"

Fingerprint *GetFingerprint(RDKit::ROMol *mol, FingerprintSelector fp)
{
    Fingerprint *result;

    FingerprintStrategy *strategy;
    switch (fp) {
    case FP_TOPOLOGICAL:
        strategy = new TopolSingleFngpr();
        break;
    case FP_TOPOLOGICAL_LAYERED_1:
        strategy = new TopolLayeredFngpr1();
        break;
    case FP_TOPOLOGICAL_LAYERED_2:
        strategy = new TopolLayeredFngpr2();
        break;
    case FP_ATOM_PAIRS:
        strategy = new AtomPairsFngpr();
        break;
    case FP_TOPOLOGICAL_TORSION:
        strategy = new TopolTorsFngpr();
        break;
    case FP_MORGAN:
        strategy = new MorganFngpr();
        break;
    default:
        strategy = new MorganFngpr();
        break;
    }

    result = strategy->GetFingerprint(mol);
    delete strategy;
    return result;
}