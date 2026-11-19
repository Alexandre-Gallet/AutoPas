/**
 * @file ContainerSelectorTest.cpp
 * @author F. Gratl
 * @date 22.06.18
 */

#include "ContainerSelectorTest.h"

#include <set>

#include "autopas/tuning/selectors/ContainerSelector.h"
#include "testingHelpers/GenerateValidConfigurations.h"
#include "testingHelpers/commonTypedefs.h"

TEST_F(ContainerSelectorTest, testSelectAndGetCurrentContainer) {
  autopas::ContainerSelectorInfo containerInfo(bBoxMin, bBoxMax, cutoff, cellSizeFactor, verletSkin, 64, 8, 8,
                                               autopas::LoadEstimatorOption::none);

  // test all individual options
  for (auto containerOp : autopas::ContainerOption::getAllOptions()) {
    const auto container = autopas::ContainerSelector<ParticleFP64>::generateContainer(containerOp, containerInfo);

    EXPECT_EQ(containerOp, container->getContainerType());
  }
}

TEST_F(ContainerSelectorTest, testReferenceVerletSearchSpace) {
  constexpr autopas::ContainerOption container = autopas::ContainerOption::verletListsReferences;
  EXPECT_FALSE(autopas::ContainerOption::getMostOptions().contains(container));
  EXPECT_TRUE(autopas::ContainerOption::getAllPairwiseOptions().contains(container));
  EXPECT_TRUE(autopas::ContainerOption::getAllTriwiseOptions().contains(container));

  for (const auto interaction : {autopas::InteractionTypeOption::pairwise, autopas::InteractionTypeOption::triwise}) {
    std::set<autopas::Configuration> expected;
    for (double cellSizeFactor : {1., 1.5}) {
      expected.emplace(container, cellSizeFactor, autopas::TraversalOption::vl_list_iteration,
                       autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::aos,
                       autopas::Newton3Option::disabled, interaction, autopas::VectorizationPatternOption::NA);
    }
    EXPECT_EQ(generateAllValidConfigurations(interaction, {container}), expected);
  }
}
