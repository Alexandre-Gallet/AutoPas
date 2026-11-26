/**
 * @file VerletListsReferencesTest.cpp
 * @author Alexandre Gallet
 * @date 18.11.2026
 */

#include "VerletListsReferencesTest.h"

#include <gmock/gmock.h>

#include <cstdint>
#include <map>
#include <set>
#include <string>
#include <type_traits>
#include <vector>

#include "autopas/LogicHandler.h"
#include "autopas/cells/FullParticleCell.h"
#include "autopas/cells/ReferenceParticleCell.h"
#include "autopas/containers/linkedCells/LinkedCellsReferences.h"
#include "autopas/containers/verletListsCellBased/verletLists/VerletLists.h"
#include "autopas/containers/verletListsCellBased/verletLists/traversals/VLListIterationTraversal.h"
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

using Backend = autopas::VerletListsReferences<ParticleFP64>;
using PairTraversal = autopas::VLListIterationTraversal<Backend::ParticleCellType, MPairwiseFunctor>;
using NeighborRows = std::map<unsigned long, std::multiset<unsigned long>>;

// Decode the CRS indices using current addresses from the owning vector. Keep duplicates visible.
void expectNeighborRows(Backend &container, const NeighborRows &expected) {
  const auto &list = container.getNeighborList();
  const auto &index = container.getParticleIndex();
  ASSERT_EQ(list.size(), expected.size());
  ASSERT_EQ(index.size(), list.size());
  std::vector<const ParticleFP64 *> byIndex(list.size(), nullptr);
  container.forEach(
      [&](const auto &p) {
        ASSERT_TRUE(index.contains(&p));
        const size_t i = index.at(&p);
        ASSERT_LT(i, byIndex.size());
        ASSERT_EQ(byIndex[i], nullptr);
        byIndex[i] = &p;
      },
      autopas::IteratorBehavior::ownedOrHaloOrDummy);
  ASSERT_THAT(byIndex, ::testing::Each(::testing::NotNull()));
  ASSERT_EQ(list.offsets.front(), 0ul);
  ASSERT_EQ(list.offsets.back(), list.indices.size());
  NeighborRows actual;
  for (size_t i = 0; i < list.size(); ++i) {
    ASSERT_LE(list.offsets[i], list.offsets[i + 1]);
    ASSERT_LE(list.offsets[i + 1], list.indices.size());
    auto &neighbors = actual[byIndex[i]->getID()];
    for (auto j : list.getNeighbors(i)) {
      ASSERT_LT(j, byIndex.size());
      neighbors.insert(byIndex[j]->getID());
    }
  }
  EXPECT_EQ(actual, expected);
}
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
  // Check that AutoPas allows this combination before testing how it behaves.
  ASSERT_TRUE(config.hasCompatibleValues());
  const autopas::ContainerSelectorInfo info{boxMin, boxMax, cutoff, GetParam(), skin, 64, 8, 100, config.loadEstimator};
  // ParticleFP64 is the particle type. The container option selects the storage backend.
  // Check that the factory created the reference-backed Verlet container.
  auto container = autopas::ContainerSelector<ParticleFP64>::generateContainer(config.container, info);
  ASSERT_TRUE(autopas::withStaticContainerType(*container, [](auto &selected) {
    return std::is_same_v<std::remove_cvref_t<decltype(selected)>, autopas::VerletListsReferences<ParticleFP64>>;
  }));
  // Assign IDs in position order so the mock functor can identify each particle it receives.
  size_t id = 0;
  for (const auto &position : positions) {
    container->addParticle(ParticleFP64{position, {}, id++});
  }
  auto traversal = autopas::TraversalSelector::generateTraversalFromConfig<ParticleFP64>(
      config, functor, container->getTraversalSelectorInfo());
  ASSERT_NE(traversal, nullptr);
  // Build the neighbor lists from the particle positions, then traverse them.
  // The traversal makes the pair or triplet calls expected by the calling test.
  container->rebuildNeighborLists(traversal.get());
  container->computeInteractions(traversal.get());
}

TEST_P(VerletListsReferencesTest, testBuildAndIteratePairwise) {
  MPairwiseFunctor functor;
  const auto first = ::testing::Property(&ParticleFP64::getID, 0ul);
  const auto second = ::testing::Property(&ParticleFP64::getID, 1ul);
  // With Newton3 disabled, each particle must be visited once as the center.
  // Checking both ID orders catches a duplicated direction as well as a missing one.
  EXPECT_CALL(functor, AoSFunctor(first, second, newton3));
  EXPECT_CALL(functor, AoSFunctor(second, first, newton3));

  // At both cell-size factors the particles lie on opposite sides of the cell boundary at x = 3.5.
  // buildAndTraverse calls computeInteractions leading to the above functor calls
  buildAndTraverse(functor, autopas::InteractionTypeOption::pairwise, {{3.25, 3., 3.}, {3.75, 3., 3.}});
}

TEST_P(VerletListsReferencesTest, testBuildAndIterateTriwise) {
  MTriwiseFunctor functor;
  // One geometric triplet gives three calls with Newton3 disabled, one per center.
  // The two neighbors may appear in either order, but all three IDs must be distinct.
  for (size_t centerID = 0; centerID < 3; ++centerID) {
    EXPECT_CALL(functor,
                AoSFunctor(::testing::Property(&ParticleFP64::getID, centerID), ::testing::_, ::testing::_, newton3))
        .WillOnce([](ParticleFP64 &center, ParticleFP64 &neighbor1, ParticleFP64 &neighbor2, bool) {
          EXPECT_THAT((std::array{center.getID(), neighbor1.getID(), neighbor2.getID()}),
                      ::testing::UnorderedElementsAre(0ul, 1ul, 2ul));
        });
  }
  // At both cell-size factors the particles lie on opposite sides of the cell boundary at x = 3.5 and the center particle lies directly on it.
  // buildAndTraverse calls computeInteractions leading to the above functor calls
  buildAndTraverse(functor, autopas::InteractionTypeOption::triwise, {{3.25, 3., 3.}, {3.75, 3., 3.}, {3.5, 3.5, 3.}});
}

TEST_P(VerletListsReferencesTest, testCutoffAndSkin) {
  // @todo Distinguish neighbor-list inclusion within the skin from force evaluation within the cutoff.
  GTEST_SKIP() << "Cutoff and skin test not implemented yet.";
}

// Move particle 30 into range of particle 20 and check that rebuilding adds the new neighbor
// relationship. Rebuild again at the same positions to catch missing or duplicated entries.
TEST_P(VerletListsReferencesTest, testRebuildTwice) {
  Backend container(boxMin, boxMax, cutoff, skin, Backend::VerletSoA, GetParam());
  MPairwiseFunctor functor;
  PairTraversal traversal(functor, dataLayout, newton3);
  ASSERT_THROW(container.computeInteractions(&traversal), autopas::utils::ExceptionHandler::AutoPasException);
  container.addParticle({{3.2, 3., 3.}, {}, 10});
  container.addParticle({{3.8, 3., 3.}, {}, 20});
  container.addParticle({{5.4, 3., 3.}, {}, 30});
  container.rebuildNeighborLists(&traversal);
  expectNeighborRows(container, {{10, {20}}, {20, {10}}, {30, {}}});

  for (auto &p : container) {
    if (p.getID() == 30) p.setR({4.7, 3., 3.});
  }
  ASSERT_TRUE(container.updateContainer(false).empty());
  const NeighborRows moved{{10, {20}}, {20, {10, 30}}, {30, {20}}};
  container.rebuildNeighborLists(&traversal);
  expectNeighborRows(container, moved);
  container.rebuildNeighborLists(&traversal);
  expectNeighborRows(container, moved);
}

// Add particles until we observe that particle 10 has changed address. The old Verlet list
// must reject traversal. After cleanup and rebuilding, each of the three owners has two neighbors.
TEST_P(VerletListsReferencesTest, testReferencesAfterGrowth) {
  Backend container(boxMin, boxMax, cutoff, skin, Backend::VerletSoA, GetParam());
  MPairwiseFunctor functor;
  PairTraversal traversal(functor, dataLayout, newton3);
  container.addParticle({{3.2, 3., 3.}, {}, 10});
  container.addParticle({{3.8, 3., 3.}, {}, 20});
  container.rebuildNeighborLists(&traversal);
  expectNeighborRows(container, {{10, {20}}, {20, {10}}});

  const auto firstOwnerAddress = [&] {
    const auto owner = container.begin(autopas::IteratorBehavior::owned);
    EXPECT_EQ(owner->getID(), 10ul);
    return reinterpret_cast<std::uintptr_t>(&*owner);
  };
  const auto original = firstOwnerAddress();
  container.addParticle({{3.5, 3.3, 3.}, {}, 30});
  unsigned long haloID = 100;
  while (firstOwnerAddress() == original) {
    ASSERT_LT(haloID, 1024ul);
    container.addHaloParticle({{.95, 2., 2.}, {}, haloID++, autopas::OwnershipState::halo});
  }
  ASSERT_THROW(container.computeInteractions(&traversal), autopas::utils::ExceptionHandler::AutoPasException);
  ASSERT_TRUE(container.updateContainer(false).empty());
  container.rebuildNeighborLists(&traversal);
  expectNeighborRows(container, {{10, {20, 30}}, {20, {10, 30}}, {30, {10, 20}}});
  // With Newton3 disabled, each of the three owners visits both neighbors: six directed calls.
  EXPECT_CALL(functor, AoSFunctor(::testing::_, ::testing::_, false)).Times(6);
  container.computeInteractions(&traversal);
}

// Erasing the interleaved halo shifts vector slots and invalidates the list. Marking owner 20
// as a dummy preserves its slot, so owners 10 and 30 can still interact through the existing list.
TEST_P(VerletListsReferencesTest, testReferencesAfterDeletion) {
  Backend container(boxMin, boxMax, cutoff, skin, Backend::VerletSoA, GetParam());
  MPairwiseFunctor functor;
  PairTraversal traversal(functor, dataLayout, newton3);
  container.addParticle({{3.1, 3., 3.}, {}, 10});
  container.addHaloParticle({{.95, 3., 3.}, {}, 11, autopas::OwnershipState::halo});
  container.addParticle({{3.3, 3., 3.}, {}, 20});
  container.addParticle({{3.7, 3., 3.}, {}, 30});
  container.rebuildNeighborLists(&traversal);

  container.deleteHaloParticles();
  ASSERT_THROW(container.computeInteractions(&traversal), autopas::utils::ExceptionHandler::AutoPasException);
  container.rebuildNeighborLists(&traversal);
  expectNeighborRows(container, {{10, {20, 30}}, {20, {10, 30}}, {30, {10, 20}}});

  auto deleted = container.begin();
  while (deleted.isValid() and deleted->getID() != 20) ++deleted;
  ASSERT_TRUE(deleted.isValid());
  const auto *slot = &*deleted;
  const auto row = container.getParticleIndex().at(slot);
  EXPECT_FALSE(container.deleteParticle(*deleted));
  ASSERT_TRUE(slot->isDummy());
  ASSERT_EQ(container.getParticleIndex().at(slot), row);
  // Existing rows can still present dummy neighbors to the force functor.
  EXPECT_CALL(functor, AoSFunctor(::testing::_, ::testing::_, false)).Times(::testing::AnyNumber());
  const auto first = ::testing::Property(&ParticleFP64::getID, 10ul);
  const auto last = ::testing::Property(&ParticleFP64::getID, 30ul);
  EXPECT_CALL(functor, AoSFunctor(first, last, false));
  EXPECT_CALL(functor, AoSFunctor(last, first, false));
  container.computeInteractions(&traversal);
  // Physically remove the dummy. The next rebuild must contain only the surviving owners, 10 and 30.
  ASSERT_TRUE(container.updateContainer(false).empty());
  ASSERT_THROW(container.computeInteractions(&traversal), autopas::utils::ExceptionHandler::AutoPasException);
  container.rebuildNeighborLists(&traversal);
  expectNeighborRows(container, {{10, {30}}, {30, {10}}});
  container.deleteAllParticles();
  container.rebuildNeighborLists(&traversal);
  expectNeighborRows(container, {});
}

// Exercise halo exchange while reusing a Verlet list. Returning halo 2 must recover its old slot
// with updated coordinates. New halo 3 goes into the buffer and must still contribute to LJ and ATM forces.
TEST_P(VerletListsReferencesTest, testHaloOwnership) {
  constexpr auto pairwise = autopas::InteractionTypeOption::pairwise;
  constexpr auto triwise = autopas::InteractionTypeOption::triwise;
  constexpr unsigned int rebuildFrequency = 10;
  const autopas::LogicHandlerInfo geometry{.boxMin = boxMin, .boxMax = boxMax, .cutoff = cutoff, .verletSkin = skin};
  const autopas::AutoTunerInfo tunerInfo{};
  auto tuners = std::make_shared<autopas::TuningManager>(tunerInfo);
  for (auto interaction : {pairwise, triwise}) {
    const std::set<autopas::Configuration> configs{
        {autopas::ContainerOption::verletListsReferences, GetParam(), autopas::TraversalOption::vl_list_iteration,
         autopas::LoadEstimatorOption::none, dataLayout, autopas::Newton3Option::disabled, interaction,
         autopas::VectorizationPatternOption::NA}};
    autopas::AutoTuner::TuningStrategiesListType strategies;
    tuners->addAutoTuner(std::make_unique<autopas::AutoTuner>(strategies, configs, tunerInfo, rebuildFrequency, ""),
                         interaction);
  }
  autopas::LogicHandler<Molecule> handler(tuners, geometry, rebuildFrequency, "", tunerInfo.aosSortingThreshold,
                                          tunerInfo.soaSortingThreshold);
  std::array<Molecule, 4> particles{Molecule{{1.35, 3., 3.}, {}, 0}, Molecule{{1.25, 3.38, 3.02}, {}, 1},
                                    Molecule{{.85, 3.1, 3.15}, {}, 2}, Molecule{{.9, 3.53, 3.21}, {}, 3}};
  particles[2].setOwnershipState(autopas::OwnershipState::halo);
  particles[3].setOwnershipState(autopas::OwnershipState::halo);
  handler.addParticle(particles[0]);
  handler.addParticle(particles[1]);
  ASSERT_TRUE(handler.updateContainer().empty());
  handler.addHaloParticle(particles[2]);
  LJFunctorType<> lj(cutoff);
  lj.setParticleProperties(24., .16);
  ATMFunctorType<> atm(cutoff);
  atm.setParticleProperties(.004);
  handler.computeInteractionsPipeline(&lj, pairwise);
  handler.computeInteractionsPipeline(&atm, triwise);

  auto &container = handler.getContainer();
  ASSERT_EQ(container.getContainerType(), autopas::ContainerOption::verletListsReferences);
  auto halo = container.begin(autopas::IteratorBehavior::halo);
  ASSERT_TRUE(halo.isValid());
  const auto *returningSlot = &*halo;
  ASSERT_TRUE(handler.updateContainer().empty());
  ASSERT_TRUE(handler.neighborListsAreValid());
  ASSERT_TRUE(returningSlot->isDummy());
  particles[2].addR({.02, 0., 0.});
  handler.addHaloParticle(particles[2]);
  handler.addHaloParticle(particles[3]);
  halo = container.begin(autopas::IteratorBehavior::halo);
  ASSERT_TRUE(halo.isValid());
  EXPECT_EQ(&*halo, returningSlot);
  EXPECT_EQ(halo->getR(), particles[2].getR());
  ++halo;
  EXPECT_FALSE(halo.isValid());
  std::vector<unsigned long> bufferedIDs;
  for (const auto &cell : std::get<1>(handler.getParticleBuffers())) {
    for (const auto &p : cell) bufferedIDs.push_back(p.getID());
  }
  ASSERT_THAT(bufferedIDs, ::testing::ElementsAre(3ul));

  const auto resetForces = [&]() {
    for (auto iter = handler.begin(autopas::IteratorBehavior::ownedOrHalo); iter.isValid(); ++iter) iter->setF({});
    for (auto &p : particles) p.setF({});
  };
  const auto checkOwnedForces = [&]() {
    std::multiset<unsigned long> owners;
    for (auto iter = handler.begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
      ASSERT_LT(iter->getID(), 2ul);
      owners.insert(iter->getID());
      for (size_t d = 0; d < 3; ++d) {
        EXPECT_NEAR(iter->getF()[d], particles[iter->getID()].getF()[d], 1e-10) << "owner " << iter->getID();
      }
    }
    EXPECT_THAT(owners, ::testing::ElementsAre(0ul, 1ul));
  };

  // Compute the expected forces by calling the same functors directly on the known pairs and triplets.
  // Agreement for each owner checks that the list and the buffer together account for the interactions.
  resetForces();
  handler.computeInteractionsPipeline(&lj, pairwise);
  for (size_t i = 0; i < 2; ++i) {
    lj.AoSFunctor(particles[i], particles[1 - i], false);
    lj.AoSFunctor(particles[i], particles[2], false);
    lj.AoSFunctor(particles[i], particles[3], false);
  }
  checkOwnedForces();

  resetForces();
  handler.computeInteractionsPipeline(&atm, triwise);
  for (size_t i = 0; i < 2; ++i) {
    atm.AoSFunctor(particles[i], particles[1 - i], particles[2], false);
    atm.AoSFunctor(particles[i], particles[1 - i], particles[3], false);
    atm.AoSFunctor(particles[i], particles[2], particles[3], false);
  }
  checkOwnedForces();
  for (auto iter = handler.begin(autopas::IteratorBehavior::halo); iter.isValid(); ++iter) {
    EXPECT_EQ(iter->getF(), (std::array<double, 3>{}));
  }
  EXPECT_TRUE(handler.neighborListsAreValid());
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
