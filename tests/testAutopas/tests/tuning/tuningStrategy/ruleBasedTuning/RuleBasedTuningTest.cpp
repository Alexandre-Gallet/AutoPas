/**
 * @file RuleBasedTuningTest.cpp
 * @author humig
 * @date 30.06.2021
 */

#include "RuleBasedTuningTest.h"

TEST_F(RuleBasedTuningTest, testParser) {
#ifdef AUTOPAS_ENABLE_RULES_BASED_AND_FUZZY_TUNING
  // This rule names both containers directly, so the parser needs no predefined variables.
  std::vector<std::pair<std::string, autopas::RuleSyntax::Define>> definitions;
  autopas::RuleSyntax::RuleBasedProgramParser parser{definitions};
  auto [tree, context] = parser.parse(R"([container="VerletListsReferences"] >= [container="VerletLists"];)");
  ASSERT_FALSE(tree.generateCode(context).instructions.empty());
  // Check that parsing preserved the direction of the rule: the reference-backed container
  // is on the greater side of >= and ordinary VerletLists is on the smaller side.
  const auto &orders = context.getConfigurationOrders();
  ASSERT_EQ(orders.size(), 1);
  EXPECT_EQ(orders.front().greater._containers,
            std::set<autopas::ContainerOption>{autopas::ContainerOption::verletListsReferences});
  EXPECT_EQ(orders.front().smaller._containers,
            std::set<autopas::ContainerOption>{autopas::ContainerOption::verletLists});
#else
  GTEST_SKIP() << "Rule-based tuning is disabled in this build.";
#endif
}
