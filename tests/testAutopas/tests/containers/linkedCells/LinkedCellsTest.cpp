/**
 * @file LinkedCellsTest.cpp
 * @author seckler
 * @date 27.04.18
 */

#include "LinkedCellsTest.h"

#include <cstdint>
#include <map>

#include "autopas/particles/OwnershipState.h"
#include "autopas/utils/ArrayUtils.h"

using LinkedCellsReferencesStorageTest = AutoPasTestBase;

TEST_F(LinkedCellsReferencesStorageTest, deletingInterleavedHaloKeepsOwnedCellReferences) {
  // construct LinkedCellsReferences object and add 3 owned and one halo particles in the following structure:
  // owned|halo|owned|owned
  //  10  | 11 |  12 | 13
  autopas::LinkedCellsReferences<ParticleFP64> storage({1., 1., 1.}, {6., 6., 6.}, 1., .2);
  storage.addParticle({{1.1, 2., 2.}, {}, 10});
  storage.addHaloParticle({{.95, 2., 2.}, {}, 11, autopas::OwnershipState::halo});
  storage.addParticle({{3., 2., 2.}, {}, 12});
  storage.addParticle({{5., 2., 2.}, {}, 13});

  // Map each owned particle ID to the cell index expected from its position.
  const std::map<unsigned long, size_t> expectedCells{{10, storage.getCellBlock().get1DIndexOfPosition({1.1, 2., 2.})},
                                                      {12, storage.getCellBlock().get1DIndexOfPosition({3., 2., 2.})},
                                                      {13, storage.getCellBlock().get1DIndexOfPosition({5., 2., 2.})}};
  // we remove halo particles which changes the structure:
  // owned|owned|owned
  //  10  | 12  | 13
  // Erasing halo 11 shifts particles 12 and 13 one slot left in the owning vector.
  // Their original cells must point to their new addresses.
  storage.deleteHaloParticles();

  // Read the owning vector first, so even a dangling cell pointer can be detected without dereferencing it.
  std::map<const ParticleFP64 *, unsigned long> live;
  storage.forEach([&](const auto &p) { live.emplace(&p, p.getID()); }, autopas::IteratorBehavior::ownedOrHaloOrDummy);
  // check correct count i.e. that one particle was deleted
  ASSERT_EQ(live.size(), expectedCells.size());
  std::map<unsigned long, size_t> actualCells;
  for (size_t i = 0; i < storage.getCells().size(); ++i) {
    for (const auto *p : storage.getCells()[i]._particles) {
      //checks that the pointer p is in the live map of particles without dereferencing it incase of a dangling pointer
      ASSERT_TRUE(live.contains(p));
      // Insertion succeeds only the first time this particle ID is found across all cells.
      EXPECT_TRUE(actualCells.emplace(live.at(p), i).second);
    }
  }
  EXPECT_EQ(actualCells, expectedCells);
}

TEST_F(LinkedCellsReferencesStorageTest, growthKeepsLeavingOwnerInItsCellUntilCollection) {
  // Record particle 10's cell and address, then move its position just outside the domain. Its original cell still holds its pointer.
  autopas::LinkedCellsReferences<ParticleFP64> storage({1., 1., 1.}, {6., 6., 6.}, 1., .2);
  storage.addParticle({{1.02, 2., 2.}, {}, 10});
  auto &originalCell = storage.getCellBlock().getContainingCell({1.02, 2., 2.});
  const auto originalAddress = reinterpret_cast<std::uintptr_t>(&originalCell[0]);
  originalCell[0].setR({.99, 2., 2.});

  // Add particles in another cell until we observe particle 10 at a new address through its original cell.
  unsigned long id = 11;
  do {
    // Fail if we reach the ID limit without observing an address change.
    ASSERT_LT(id, 1000ul);
    // add a new  particle in a different cell from particle p_id = 10
    ASSERT_NO_THROW(storage.addParticle({{4., 2., 2.}, {}, id++}));
    // sanity check that the cell with particle p_id = 10 still only has one particle
    ASSERT_EQ(originalCell.size(), 1ul);
  } while (reinterpret_cast<std::uintptr_t>(&originalCell[0]) == originalAddress);
  EXPECT_EQ(originalCell[0].getID(), 10ul);
  // Collect particle 10 and remove its entry from storage. A second update should find no departures.
  const auto leaving = storage.updateContainer(false);
  ASSERT_EQ(leaving.size(), 1ul);
  EXPECT_EQ(leaving.front().getID(), 10ul);
  EXPECT_TRUE(leaving.front().isOwned());
  EXPECT_TRUE(originalCell.isEmpty());
  EXPECT_EQ(storage.getNumberOfParticles(autopas::IteratorBehavior::owned), id - 11);
  EXPECT_TRUE(storage.updateContainer(false).empty());
}

TYPED_TEST_SUITE_P(LinkedCellsTest);

TYPED_TEST_P(LinkedCellsTest, testUpdateContainer) {
  using namespace autopas::utils::ArrayMath::literals;

  const std::array<double, 3> zero{0., 0., 0.};
  const std::array<double, 3> boxMin = zero;
  const std::array<double, 3> boxMax{4.5, 4.5, 4.5};
  // set values so we have 3x3x3 cells + halo = 5x5x5
  const double cutoff{1.0};
  const double skin{0.5};
  typename TestFixture::LinkedCellsType linkedCells(boxMin, boxMax, cutoff, skin);

  // create owned particles
  const std::vector<ParticleFP64> ownedParticles{
      // clang-format off
      {{0.5, 0.5, 0.5}, zero, 0},
      {{1.5, 1.5, 1.5}, zero, 1},
      {{1.6, 1.5, 1.5}, zero, 2},
      {{4.4, 1.5, 1.5}, zero, 3},
      {{4.0, 4.0, 4.0}, zero, 4},
      // clang-format on
  };

  // These are going to be halo particles
  const std::vector<ParticleFP64> haloParticles{
      {{-0.5, +1.5, +1.5}, zero, 5, autopas::OwnershipState::halo},
      {{+5.0, +1.5, +1.5}, zero, 6, autopas::OwnershipState::halo},
      {{+1.5, -0.5, +1.5}, zero, 7, autopas::OwnershipState::halo},
      {{+1.5, +1.5, -0.5}, zero, 8, autopas::OwnershipState::halo},
  };

  // calculate the cell IDs for each particle
  std::unordered_map<size_t, size_t> particleIdToCellIdMap;
  for (const auto &particleCollection : {ownedParticles, haloParticles}) {
    for (const auto &particle : particleCollection) {
      const auto cellID = linkedCells.getCellBlock().get1DIndexOfPosition(particle.getR());
      particleIdToCellIdMap[particle.getID()] = cellID;
    }
  }
  ASSERT_EQ(particleIdToCellIdMap.size(), ownedParticles.size() + haloParticles.size())
      << "There should be exactly one cellID per particleID.\n"
         "Either the test is set up wrong or get1DIndexOfPosition is broken.";

  // we insert owned and halo particles alternating. This way we can check if references are updated correctly when
  // using LinkedCellsReferences
  linkedCells.addParticle(ownedParticles[0]);
  linkedCells.addHaloParticle(haloParticles[0]);
  linkedCells.addParticle(ownedParticles[1]);
  linkedCells.addHaloParticle(haloParticles[1]);
  linkedCells.addParticle(ownedParticles[2]);
  linkedCells.addHaloParticle(haloParticles[2]);
  linkedCells.addParticle(ownedParticles[3]);
  linkedCells.addHaloParticle(haloParticles[3]);
  linkedCells.addParticle(ownedParticles[4]);

  this->checkParticleIDsInCells(
      linkedCells,
      {{particleIdToCellIdMap[8], {{8, autopas::OwnershipState::halo}}},
       {particleIdToCellIdMap[0], {{0, autopas::OwnershipState::owned}}},
       {particleIdToCellIdMap[7], {{7, autopas::OwnershipState::halo}}},
       {particleIdToCellIdMap[5], {{5, autopas::OwnershipState::halo}}},
       {particleIdToCellIdMap[3], {{3, autopas::OwnershipState::owned}}},
       {particleIdToCellIdMap[1], {{1, autopas::OwnershipState::owned}, {2, autopas::OwnershipState::owned}}},
       {particleIdToCellIdMap[6], {{6, autopas::OwnershipState::halo}}},
       {particleIdToCellIdMap[4], {{4, autopas::OwnershipState::owned}}}},
      true, __LINE__);

  // // new locations for owned particles
  linkedCells.getCells()[particleIdToCellIdMap[0]].begin()->addR({+2.0, +0.0, +0.0});  // move to {1.5, 0.5, 0.5}
  linkedCells.getCells()[particleIdToCellIdMap[1]].begin()->addR({-1.0, -0.0, -0.0});  // move to {-0.5, 1.5, 1.5}
  linkedCells.getCells()[particleIdToCellIdMap[3]].begin()->addR({+0.2, +0.0, -1.0});  // move to {5.0, 1.5, 0.5}
  linkedCells.getCells()[particleIdToCellIdMap[4]].begin()->addR({-0.9, -2.0, -2.0});  // move to {1.6, 0.5, 0.5}

  std::vector<ParticleFP64> invalidParticles;
  EXPECT_NO_THROW(invalidParticles = linkedCells.updateContainer(this->_keepListsValid));
  EXPECT_EQ(invalidParticles.size(), 1);
  EXPECT_EQ(invalidParticles[0].getID(), 3);

  if (this->_keepListsValid) {
    // if the lists are kept valid, particles are NOT moved between cells!
    // halo particles should now be dummies
    // particle 3 should be a leaving particle and therefore a dummy
    this->checkParticleIDsInCells(
        linkedCells,
        {{particleIdToCellIdMap[8], {{8, autopas::OwnershipState::dummy}}},
         {particleIdToCellIdMap[0], {{0, autopas::OwnershipState::owned}}},
         {particleIdToCellIdMap[7], {{7, autopas::OwnershipState::dummy}}},
         {particleIdToCellIdMap[5], {{5, autopas::OwnershipState::dummy}}},
         {particleIdToCellIdMap[3], {{3, autopas::OwnershipState::dummy}}},
         {particleIdToCellIdMap[1], {{1, autopas::OwnershipState::owned}, {2, autopas::OwnershipState::owned}}},
         {particleIdToCellIdMap[6], {{6, autopas::OwnershipState::dummy}}},
         {particleIdToCellIdMap[4], {{4, autopas::OwnershipState::owned}}}},
        true, __LINE__);
  } else {
    // if the lists are not kept valid, particles should be moved between cells, so update the cells!
    // halo particles should be removed by updateContainer() at this point
    this->checkParticleIDsInCells(
        linkedCells,
        {
            {particleIdToCellIdMap[0] + 1, {{0, autopas::OwnershipState::owned}}},  // moved one cell to the right
            {particleIdToCellIdMap[1] - 1, {{1, autopas::OwnershipState::owned}}},  // moved one cell to the left
            {particleIdToCellIdMap[2], {{2, autopas::OwnershipState::owned}}},      // didn't change cell
            {particleIdToCellIdMap[4] - (0 + 1 * 5 + 1 * 5 * 5),
             {{4, autopas::OwnershipState::owned}}},  // moved one cell to the front and one down
        },
        false /*here, we do not know the order!*/, __LINE__);
  }
}

TYPED_TEST_P(LinkedCellsTest, testUpdateContainerCloseToBoundary) {
  const std::array<double, 3> boxMin{0., 0., 0.};
  const std::array<double, 3> boxMax{10., 10., 10.};
  const double cutoff{1.5};
  const double skin{1.};  // particles are moved by up to 0.5 and lists might be kept valid
  typename TestFixture::LinkedCellsType linkedCells(boxMin, boxMax, cutoff, skin);

  int id = 1;
  for (const double x : {0., 5., 9.999}) {
    for (const double y : {0., 5., 9.999}) {
      for (const double z : {0., 5., 9.999}) {
        const ParticleFP64 p({x, y, z}, {0., 0., 0.}, id++);
        EXPECT_NO_THROW(linkedCells.addParticle(p));  // inside, therefore ok!
      }
    }
  }
  std::set<unsigned long> movedIDs;
  // we move particles that are close to the boundary to outside the container and remember their IDs
  for (auto iter = linkedCells.begin(); iter.isValid(); ++iter) {
    for (unsigned short dim = 0; dim < 3; ++dim) {
      if (iter->getR()[dim] < 0.5) {
        auto r = iter->getR();
        // smallest double smaller than 0
        r[dim] = std::nexttoward(0., -1.);
        iter->setR(r);
        movedIDs.insert(iter->getID());
      }
      if (iter->getR()[dim] > 9.5) {
        auto r = iter->getR();
        r[dim] = 10.;
        iter->setR(r);
        movedIDs.insert(iter->getID());
      }
    }
  }

  // now update the container!
  const auto invalidParticles = linkedCells.updateContainer(this->_keepListsValid);
  // the particles should no longer be in the inner cells!
  for (auto iter = linkedCells.begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
    EXPECT_EQ(movedIDs.count(iter->getID()), 0)
        << "Particle " << iter->getID() << " at " << autopas::utils::ArrayUtils::to_string(iter->getR())
        << " is still in an inner cell although it was moved!";
  }

  // the particles should now be inside the invalidParticles vector!
  EXPECT_EQ(movedIDs.size(), invalidParticles.size());
  for (const auto &particle : invalidParticles) {
    EXPECT_EQ(movedIDs.count(particle.getID()), 1)
        << "Particle " << particle.getID() << " at " << autopas::utils::ArrayUtils::to_string(particle.getR())
        << " was not returned by updateContainer()!";
  }
}

REGISTER_TYPED_TEST_SUITE_P(LinkedCellsTest, testUpdateContainer, testUpdateContainerCloseToBoundary);

// Workaround for storing two types.
// Currently, clang produces bugs if one tries to store this using a tuple or a pair.
// This problem is described in P0641R2 and occurs if an explicitly defaulted constructor cannot be instantiated.
// This is fixed in c++20.
template <typename first, typename second>
struct two_values {
  using first_t = first;
  using second_t = second;
};

// defines the types of linkedCells and _keepListsValid
struct LC_KeepListsValid : two_values<autopas::LinkedCells<ParticleFP64>, std::true_type> {};
struct LC_DontKeepListsValid : two_values<autopas::LinkedCells<ParticleFP64>, std::false_type> {};
struct LCRef_KeepListsValid : two_values<autopas::LinkedCellsReferences<ParticleFP64>, std::true_type> {};
struct LCRef_DontKeepListsValid : two_values<autopas::LinkedCellsReferences<ParticleFP64>, std::false_type> {};

using MyTypes =
    ::testing::Types<LC_KeepListsValid, LC_DontKeepListsValid, LCRef_KeepListsValid, LCRef_DontKeepListsValid>;

/// @todo c++20: replace with:
// using MyTypes = ::testing::Types<std::tuple<autopas::LinkedCells<ParticleFP64>, std::true_type>,
//                                  std::tuple<autopas::LinkedCells<ParticleFP64>, std::false_type>,
//                                  std::tuple<autopas::LinkedCellsReferences<ParticleFP64>, std::true_type>,
//                                  std::tuple<autopas::LinkedCellsReferences<ParticleFP64>, std::false_type> >;

INSTANTIATE_TYPED_TEST_SUITE_P(GeneratedTyped, LinkedCellsTest, MyTypes);
