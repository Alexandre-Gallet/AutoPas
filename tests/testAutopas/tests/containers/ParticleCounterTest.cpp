/**
 * @file ParticleCounterTest.cpp
 * @author D. Martin
 * @date 10.09.23
 */

#include "ParticleCounterTest.h"

#include "autopas/tuning/selectors/ContainerSelector.h"

/**
 * Lambda to generate a readable string out of the parameter of this test.
 */
static auto toString = [](const auto &info) {
  auto containerOption = info.param;
  std::stringstream resStream;
  resStream << containerOption.to_string();
  std::string res = resStream.str();
  return res;
};

/**
 * Tests for all containers if the method getNumberOfParticles(IteratorBehavior behavior) returns the correct number of
 * particles internally stored in cells with respect to different IteratorBehaviors.
 *
 */
TEST_P(ParticleCounterTest, testGetNumberOfParticles) {
  const auto containerOption = GetParam();

  // Construct container
  constexpr double skin = _cutoff * 0.1;

  const auto containerInfo = autopas::ContainerSelectorInfo{
      _boxMin, _boxMax, _cutoff, _cellSizeFactor, skin, 32, 8, 8, autopas::LoadEstimatorOption::none};
  auto container = autopas::ContainerSelector<Molecule>::generateContainer(containerOption, containerInfo);

  // add owned particles
  for (int i = 0; i < _numOwnedParticles; i++) {
    Molecule p({1., 1., 1.}, {0., 0., 0.}, i);
    p.setOwnershipState(autopas::OwnershipState::owned);
    container->addParticle(p);
  }

  // add halo particles
  for (int i = 0; i < _numHaloParticles; i++) {
    Molecule p({_boxMax[0] + 1., _boxMax[1] + 1., _boxMax[2] + 1.}, {0., 0., 0.}, _numOwnedParticles + i);
    p.setOwnershipState(autopas::OwnershipState::halo);
    container->addHaloParticle(p);
  }

  for (int i = 0; i < _numDummyParticles; i++) {
    Molecule p({1., 1., 1.}, {0., 0., 0.}, _numOwnedParticles + _numHaloParticles + i);
    p.setOwnershipState(autopas::OwnershipState::dummy);
    container->addParticle(p);
  }

  // LinkedCellsReferences omits particles that are already dummy when inserting their cell references.
  const bool usesReferenceCells = containerOption == autopas::ContainerOption::linkedCellsReferences or
                                  containerOption == autopas::ContainerOption::verletListsReferences;
  const auto expectedDummyParticles = usesReferenceCells ? 0u : _numDummyParticles;

  // by default getNumberOfParticles() should return the number of owned particles
  EXPECT_EQ(container->getNumberOfParticles(), _numOwnedParticles);
  EXPECT_EQ(container->getNumberOfParticles(autopas::options::IteratorBehavior::owned), _numOwnedParticles);
  EXPECT_EQ(container->getNumberOfParticles(autopas::options::IteratorBehavior::halo), _numHaloParticles);
  EXPECT_EQ(container->getNumberOfParticles(autopas::options::IteratorBehavior::ownedOrHalo),
            _numOwnedParticles + _numHaloParticles);
  EXPECT_EQ(container->getNumberOfParticles(autopas::options::IteratorBehavior::dummy), expectedDummyParticles);
  EXPECT_EQ(container->getNumberOfParticles(autopas::options::IteratorBehavior::ownedOrHaloOrDummy),
            _numOwnedParticles + _numHaloParticles + expectedDummyParticles);
}

INSTANTIATE_TEST_SUITE_P(Generated, ParticleCounterTest, ::testing::ValuesIn(autopas::ContainerOption::getAllOptions()),
                         toString);
