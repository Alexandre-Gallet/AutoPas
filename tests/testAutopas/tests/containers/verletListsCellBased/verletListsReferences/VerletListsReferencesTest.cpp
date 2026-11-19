/**
 * @file VerletListsReferencesTest.cpp
 * @author Alexandre Gallet
 * @date 18.11.2026
 */

#include "VerletListsReferencesTest.h"

#include <gmock/gmock.h>

#include <string>
#include <type_traits>

#include "autopas/cells/FullParticleCell.h"
#include "autopas/cells/ReferenceParticleCell.h"
#include "autopas/containers/linkedCells/LinkedCellsReferences.h"
#include "autopas/containers/verletListsCellBased/verletLists/VerletLists.h"
#include "autopas/containers/verletListsCellBased/verletListsReferences/VerletListsReferences.h"
#include "autopas/tuning/selectors/ContainerSelector.h"
#include "autopas/tuning/selectors/TraversalSelector.h"
#include "autopas/utils/StaticContainerSelector.h"
#include "testingHelpers/commonTypedefs.h"

namespace {
using DefaultVerletLists = autopas::VerletLists<ParticleFP64>;
using ReferenceVerletLists = autopas::VerletLists<ParticleFP64, autopas::LinkedCellsReferences<ParticleFP64>>;

static_assert(std::is_same_v<DefaultVerletLists::ParticleCellType, autopas::FullParticleCell<ParticleFP64>>);
static_assert(std::is_same_v<ReferenceVerletLists::ParticleCellType, autopas::ReferenceParticleCell<ParticleFP64>>);
}  // namespace

template <class Functor_T>
void VerletListsReferencesTest::buildAndTraverse(Functor_T &functor, autopas::InteractionTypeOption interactionType,
                                                 std::initializer_list<std::array<double, 3>> positions) const {
  const autopas::Configuration config{autopas::ContainerOption::verletListsReferences,
                                      GetParam(),
                                      autopas::TraversalOption::vl_list_iteration,
                                      autopas::LoadEstimatorOption::none,
                                      dataLayout,
                                      autopas::Newton3Option::disabled,
                                      interactionType,
                                      autopas::VectorizationPatternOption::NA};
  ASSERT_TRUE(config.hasCompatibleValues());
  const autopas::ContainerSelectorInfo info{boxMin, boxMax, cutoff, GetParam(), skin, 64, 8, 100, config.loadEstimator};
  auto container = autopas::ContainerSelector<ParticleFP64>::generateContainer(config.container, info);
  ASSERT_TRUE(autopas::withStaticContainerType(*container, [](auto &selected) {
    return std::is_same_v<std::remove_cvref_t<decltype(selected)>, autopas::VerletListsReferences<ParticleFP64>>;
  }));
  size_t id = 0;
  for (const auto &position : positions) {
    container->addParticle(ParticleFP64{position, {}, id++});
  }
  auto traversal = autopas::TraversalSelector::generateTraversalFromConfig<ParticleFP64>(
      config, functor, container->getTraversalSelectorInfo());
  ASSERT_NE(traversal, nullptr);
  container->rebuildNeighborLists(traversal.get());
  container->computeInteractions(traversal.get());
}

TEST_P(VerletListsReferencesTest, testBuildAndIteratePairwise) {
  MPairwiseFunctor functor;
  const auto first = ::testing::Property(&ParticleFP64::getID, 0ul);
  const auto second = ::testing::Property(&ParticleFP64::getID, 1ul);
  EXPECT_CALL(functor, AoSFunctor(first, second, newton3));
  EXPECT_CALL(functor, AoSFunctor(second, first, newton3));

  // At both cell-size factors the particles lie on opposite sides of the cell boundary at x = 3.5.
  buildAndTraverse(functor, autopas::InteractionTypeOption::pairwise, {{3.25, 3., 3.}, {3.75, 3., 3.}});
}

TEST_P(VerletListsReferencesTest, testBuildAndIterateTriwise) {
  MTriwiseFunctor functor;
  for (size_t centerID = 0; centerID < 3; ++centerID) {
    EXPECT_CALL(functor,
                AoSFunctor(::testing::Property(&ParticleFP64::getID, centerID), ::testing::_, ::testing::_, newton3))
        .WillOnce([](ParticleFP64 &center, ParticleFP64 &neighbor1, ParticleFP64 &neighbor2, bool) {
          EXPECT_THAT((std::array{center.getID(), neighbor1.getID(), neighbor2.getID()}),
                      ::testing::UnorderedElementsAre(0ul, 1ul, 2ul));
        });
  }

  buildAndTraverse(functor, autopas::InteractionTypeOption::triwise, {{3.25, 3., 3.}, {3.75, 3., 3.}, {3.5, 3.5, 3.}});
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
