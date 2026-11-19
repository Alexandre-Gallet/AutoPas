/**
 * @file VerletListsReferencesTest.cpp
 * @author Alexandre Gallet
 * @date 18.11.2026
 */

#include "VerletListsReferencesTest.h"

#include <string>
#include <type_traits>

#include "autopas/cells/FullParticleCell.h"
#include "autopas/cells/ReferenceParticleCell.h"
#include "autopas/containers/linkedCells/LinkedCellsReferences.h"
#include "autopas/containers/verletListsCellBased/verletLists/VerletLists.h"
#include "testingHelpers/commonTypedefs.h"

namespace {
using DefaultVerletLists = autopas::VerletLists<ParticleFP64>;
using ReferenceVerletLists = autopas::VerletLists<ParticleFP64, autopas::LinkedCellsReferences<ParticleFP64>>;

static_assert(std::is_same_v<DefaultVerletLists::ParticleCellType, autopas::FullParticleCell<ParticleFP64>>);
static_assert(std::is_same_v<ReferenceVerletLists::ParticleCellType, autopas::ReferenceParticleCell<ParticleFP64>>);
}  // namespace

// @todo Include VerletListsReferences.h once the container exists and construct it in each test.
// Use boxMin, boxMax, cutoff, skin, and GetParam() for its geometry and cell-size factor.
// commonTypedefs.h supplies ParticleFP64, Molecule, mock functors, LJFunctorType, and ATMFunctorType.
// Derive traversal cell types from the new container's ParticleCellType rather than assuming FPCell.
// Replace each GTEST_SKIP() with setup, operation, and assertions when implementing that test.

TEST_P(VerletListsReferencesTest, testBuildAndIteratePairwise) {
  // @todo Rebuild and traverse a two-particle list; check both directed calls with Newton3 disabled.
  GTEST_SKIP() << "Pairwise neighbor-list test not implemented yet.";
}

TEST_P(VerletListsReferencesTest, testBuildAndIterateTriwise) {
  // @todo Rebuild and traverse a three-particle list; check triplet multiplicity and force-update ownership.
  GTEST_SKIP() << "Three-body neighbor-list test not implemented yet.";
}

TEST_P(VerletListsReferencesTest, testCutoffAndSkin) {
  // @todo Distinguish neighbor-list inclusion within the skin from force evaluation within the cutoff.
  GTEST_SKIP() << "Cutoff and skin test not implemented yet.";
}

TEST_P(VerletListsReferencesTest, testRebuildTwice) {
  // @todo Rebuild twice and verify that obsolete or duplicate neighbor entries do not accumulate.
  GTEST_SKIP() << "Repeated-rebuild test not implemented yet.";
}

TEST_P(VerletListsReferencesTest, testReferencesAfterGrowth) {
  // @todo Force particle-storage growth, perform the required rebuild, and verify references and interactions.
  GTEST_SKIP() << "Storage-growth reference test not implemented yet.";
}

TEST_P(VerletListsReferencesTest, testReferencesAfterDeletion) {
  // @todo Delete particles, perform the required update/rebuild, and verify surviving IDs and interactions.
  GTEST_SKIP() << "Deletion reference test not implemented yet.";
}

TEST_P(VerletListsReferencesTest, testHaloOwnership) {
  // @todo Combine owned and halo particles; check interaction counting and which forces may be updated.
  GTEST_SKIP() << "Halo ownership test not implemented yet.";
}

TEST_P(VerletListsReferencesTest, testLJForces) {
  // @todo Compare forces by particle ID against a trusted pairwise reference using a nontrivial configuration.
  GTEST_SKIP() << "LJ force comparison not implemented yet.";
}

TEST_P(VerletListsReferencesTest, testATMForces) {
  // @todo Compare forces by particle ID against a trusted three-body reference using a non-collinear triplet.
  GTEST_SKIP() << "ATM force comparison not implemented yet.";
}

INSTANTIATE_TEST_SUITE_P(Generated, VerletListsReferencesTest, ::testing::Values(1.0, 2.0),
                         [](const ::testing::TestParamInfo<double> &info) {
                           return "CellSizeFactor_" + std::to_string(static_cast<int>(info.param));
                         });
