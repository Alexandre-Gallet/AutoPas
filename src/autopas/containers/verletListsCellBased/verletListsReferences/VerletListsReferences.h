/**
 * @file VerletListsReferences.h
 */

#pragma once

#include "autopas/containers/linkedCells/LinkedCellsReferences.h"
#include "autopas/containers/verletListsCellBased/verletLists/VerletLists.h"

namespace autopas {

/**
 * Verlet lists with a single particle vector and linked cells holding references into that vector.
 * @tparam Particle_T Particle type stored in the container.
 */
template <class Particle_T>
class VerletListsReferences : public VerletLists<Particle_T, LinkedCellsReferences<Particle_T>> {
 public:
  /** Inherit the geometry and neighbor-list build options from VerletLists. */
  using VerletLists<Particle_T, LinkedCellsReferences<Particle_T>>::VerletLists;

  /** @copydoc ParticleContainerInterface::getContainerType() */
  [[nodiscard]] ContainerOption getContainerType() const override { return ContainerOption::verletListsReferences; }
};

}  // namespace autopas
