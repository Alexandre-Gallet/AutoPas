/**
 * @file VerletListsReferences.h
 */

#pragma once

#include "autopas/containers/linkedCells/LinkedCellsReferences.h"
#include "autopas/containers/verletListsCellBased/verletLists/VerletLists.h"
#include "autopas/utils/ExceptionHandler.h"

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

  /** Refresh cell references before the Verlet builder assigns indices to particle addresses. */
  void rebuildNeighborLists(TraversalInterface *traversal) override {
    this->_neighborListIsValid.store(false, std::memory_order_relaxed);
    this->_linkedCells.updateDirtyParticleReferences();
    VerletLists<Particle_T, LinkedCellsReferences<Particle_T>>::rebuildNeighborLists(traversal);
  }

  /** Traverse the last successful rebuild while its particle addresses and indices remain valid. */
  void computeInteractions(TraversalInterface *traversal) override {
    if (not this->_neighborListIsValid.load(std::memory_order_relaxed)) {
      utils::ExceptionHandler::exception(
          "VerletListsReferences::computeInteractions() requires a valid neighbor list. "
          "Call rebuildNeighborLists() after changing particle storage.");
      return;
    }
    VerletLists<Particle_T, LinkedCellsReferences<Particle_T>>::computeInteractions(traversal);
  }

  /** @copydoc ParticleContainerInterface::getContainerType() */
  [[nodiscard]] ContainerOption getContainerType() const override { return ContainerOption::verletListsReferences; }
};

}  // namespace autopas
