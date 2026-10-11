/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#ifndef INFOMAP_OPTIMIZER_BASE_H_
#define INFOMAP_OPTIMIZER_BASE_H_

#include "InfomapBase.h"
#include "InfoNode.h"
#include "FlowData.h"
#include <vector>

namespace infomap {

class StateNetwork;

class InfomapOptimizerBase {
  friend class InfomapBase;
  using FlowDataType = FlowData;

public:
  InfomapOptimizerBase() = default;

  virtual ~InfomapOptimizerBase() = default;

  virtual void init(InfomapBase* infomap) = 0;

  // ===================================================
  // IO
  // ===================================================

  virtual std::ostream& toString(std::ostream& out) const = 0;

  // ===================================================
  // Getters
  // ===================================================

  virtual double getCodelength() const = 0;

  virtual double getIndexCodelength() const = 0;

  virtual double getModuleCodelength() const = 0;

  virtual double getMetaCodelength(bool /*unweighted*/ = false) const { return 0.0; }

#if INFOMAP_FEATURE_LOSSY_MAP_EQUATION
  // Lossy map equation metrics; zero for every other objective.
  virtual double getLossyRate() const { return 0.0; } // L_lossy (bits/step)
  virtual double getLossyDistortion() const { return 0.0; } // D = sum over noise modules of H_i
  virtual double getLossyOneLevelLossless() const { return 0.0; } // L_1 = H(p_alpha), lambda-independent
#endif

  // Forward per-network properties (e.g. total degree, node count) to the objective.
  // Only meaningful for objectives that use them (BiasedMapEquation); a no-op otherwise.
  // Stored per-instance so parallel trial workers don't share mutable state.
  virtual void setNetworkProperties(const StateNetwork& /*network*/) {}

  // Propagate already-computed network properties from a parent optimizer to a sub/super
  // Infomap instance (which has no full StateNetwork of its own). No-op unless the objective
  // uses them (BiasedMapEquation).
  virtual void inheritNetworkPropertiesFrom(const InfomapOptimizerBase& /*parent*/) {}

  // Propagate the objective's own parameters -- what it minimises, as opposed to the network
  // it minimises over -- from a parent optimizer. Needed by the super-level instance, which is
  // built without memory and therefore from a default Config. No-op unless the objective has
  // such parameters (BiasedMapEquation). See getSuperInfomap.
  virtual void inheritObjectiveParametersFrom(const InfomapOptimizerBase& /*parent*/) {}

  // Whether this objective implements the parameters above: the entropy bias correction and the
  // preferred-module-count cost. Only BiasedMapEquation does, so --entropy-corrected and
  // --preferred-number-of-modules are inert on every other objective. Declared by the objective
  // rather than derived from the input flags so the run path cannot disagree with
  // InfomapBase::initOptimizer about which objective is in use (#904).
  virtual bool implementsObjectiveParameters() const { return false; }

protected:
  virtual unsigned int numActiveModules() const = 0;

  // ===================================================
  // Run: Init: *
  // ===================================================

  // Init terms that is constant for the whole network
  virtual void initTree() = 0;

  virtual void initNetwork() = 0;

  virtual void initSuperNetwork() = 0;

  virtual double calcCodelength(const InfoNode& parent) const = 0;

  // ===================================================
  // Run: Partition: *
  // ===================================================

  virtual void initPartition() = 0;

  virtual void moveActiveNodesToPredefinedModules(std::vector<unsigned int>& modules) = 0;

  virtual unsigned int optimizeActiveNetwork() = 0;

  virtual unsigned int tryMoveEachNodeIntoBestModule() = 0;

  // virtual unsigned int tryMoveEachNodeIntoBestModuleLocal() = 0;

  virtual unsigned int tryMoveEachNodeIntoBestModuleInParallel() = 0;

  virtual void consolidateModules(bool replaceExistingModules = true) = 0;

  virtual bool restoreConsolidatedOptimizationPointIfNoImprovement(bool forceRestore = false) = 0;

  // ===================================================
  // Debug: *
  // ===================================================

  virtual void printDebug() = 0;
};

} /* namespace infomap */

#endif // INFOMAP_OPTIMIZER_BASE_H_
