/**
 * @file ParticleVector.h
 *
 * @author lunaticcoding
 * @date 04.05.2020
 */

#pragma once

#include <algorithm>
#include <mutex>
#include <span>
#include <vector>

#include "autopas/utils/WrapOpenMP.h"

/**
 * ParticleVector class.
 * This class uses a std::vector to store particles and allows access to its references.
 * It also keeps track, which references are no longer valid, which happens
 * when operations like delete and resize are performed.
 * @tparam Particle_T of the Particle that is being stored
 */
template <class Particle_T>
class ParticleVector {
 public:
  ParticleVector() = default;

  /**
   * Returns the dirty flag, indicating whether Particles exist in the vector that are not stored in a cell yet.
   * @return True if dirty, false otherwise
   */
  bool isDirty() { return _dirty; }

  /**
   * Marks the ParticleVector as clean. Should be called after updating the references of dirty particles.
   */
  void markAsClean() {
    _dirty = false;
    _dirtyIndex = _particleListImp.size();
  }

  /**
   * Remove all particles from the container and mark it as dirty.
   */
  void clearAllParticles() {
    _particleListImp.clear();
    _dirty = true;
    _dirtyIndex = 0;
  }

  /**
   * Remove all halo particles from the container and mark it as dirty.
   */
  void clearHaloParticles() {
    eraseIf([](const auto &particle) { return particle.isHalo(); });
  }

  /**
   * Remove all dummy particles from the container and mark it as dirty.
   */
  void deleteDummyParticles() {
    eraseIf([](const auto &particle) { return particle.isDummy(); });
  }

  /**
   * Compact the vector in order and request a full rebuild of the cell references.
   * @param remove Predicate selecting particles to erase.
   */
  template <class Predicate>
  void eraseIf(Predicate remove) {
    std::erase_if(_particleListImp, remove);
    _dirty = true;
    _dirtyIndex = 0;
  }

  /**
   * Append a particle while holding the same lock for the capacity check and the append.
   * After a successful append, mark the new slot dirty. If the vector reallocated,
   * every existing particle has a new address, so the dirty range starts at slot zero.
   * @param particle Particle to append.
   */
  void push_back(const Particle_T &particle) {
    std::lock_guard guard(_particleListLock);
    const bool relocates = willReallocateOnAppend();
    _particleListImp.push_back(particle);
    _dirty = true;
    if (relocates) {
      _dirtyIndex = 0;
    }
  }

  /** @return Whether appending the next particle will relocate existing particles. */
  [[nodiscard]] bool willReallocateOnAppend() const { return _particleListImp.size() == _particleListImp.capacity(); }

  /** @return A view of the current particle slots, valid until the vector is structurally modified. */
  [[nodiscard]] std::span<Particle_T> particles() { return _particleListImp; }

  /**
   * Roll back the last append if inserting the new particle's cell pointer failed.
   * @pre The last particle has no cell pointer. All earlier cell pointers are valid.
   * Removing the last slot preserves the addresses of all earlier particles.
   */
  void undoLastAppend() {
    _particleListImp.pop_back();
    markAsClean();
  }

  /**
   * Get the number of Particles in the data structure.
   * @return Total number of Particles
   */
  int totalSize() { return _particleListImp.size(); }

  /**
   * Get the number of dirty Particles in the data structure.
   * @return Number of dirty Particles
   */
  int dirtySize() { return totalSize() - _dirtyIndex; }

  /**
   * Indicates, whether References already stored in cells need to be updated.
   * @return Boolean indicating the above
   */
  bool needsRebuild() { return _dirtyIndex == 0; }

  /**
   * Begin of the iterator over dirty Particles
   * @return Start of the iterator
   */
  auto beginDirty() { return _particleListImp.begin() + _dirtyIndex; }
  /**
   * End of the iterator over dirty Particles
   * @return End of the iterator
   */
  auto endDirty() { return _particleListImp.end(); }

  /**
   * Iterate over all particles and execute lambda function on them.
   * @param forEachLambda lambda function that defines how to modify particles
   */
  template <typename Lambda>
  void forEach(Lambda forEachLambda) {
    for (Particle_T &p : _particleListImp) {
      forEachLambda(p);
    }
  }

  /**
   * Iterate over all particles and execute lambda function on them.
   * @param reduceLambda lambda function that defines how to modify particles
   * @param result reference to variable where result will be stored after operation
   */
  template <typename Lambda, typename A>
  void reduce(Lambda reduceLambda, A &result) {
    for (Particle_T &p : _particleListImp) {
      reduceLambda(p, result);
    }
  }

 private:
  /**
   * Flag indicating whether there are out-of-date references in the vector.
   */
  bool _dirty = false;
  /**
   * Index of the first out-of-date reference.
   */
  size_t _dirtyIndex{0};
  autopas::AutoPasLock _particleListLock;
  std::vector<Particle_T> _particleListImp;
};
