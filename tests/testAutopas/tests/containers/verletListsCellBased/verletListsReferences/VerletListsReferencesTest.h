/**
 * @file VerletListsReferencesTest.h
 * @autho Alexandre Gallet
 * @date 18.11.2026
 */

#pragma once

#include <array>
#include <initializer_list>

#include "AutoPasTestBase.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/options/InteractionTypeOption.h"

/**
 * Test fixture for the reference-backed Verlet container.
 *
 * GetParam() supplies the cell-size factor. The initial backend targets AoS with Newton3 disabled.
 */
class VerletListsReferencesTest : public AutoPasTestBase, public ::testing::WithParamInterface<double> {
 protected:
  /** Lower corner of the test domain. */
  static constexpr std::array<double, 3> boxMin{1., 1., 1.};
  /** Upper corner of the test domain. */
  static constexpr std::array<double, 3> boxMax{6., 6., 6.};
  /** Interaction cutoff. */
  static constexpr double cutoff{1.};
  /** Verlet skin. */
  static constexpr double skin{0.2};
  /** Initial force-evaluation layout. */
  static constexpr autopas::DataLayoutOption dataLayout{autopas::DataLayoutOption::aos};
  /** Initial Newton3 setting. */
  static constexpr bool newton3{false};

  template <class Functor_T>
  void buildAndTraverse(Functor_T &functor, autopas::InteractionTypeOption interactionType,
                        std::initializer_list<std::array<double, 3>> positions) const;
};
