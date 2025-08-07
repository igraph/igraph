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

  void setInterruptionHandler(interruptionHandlerFn *interruptionHandler) {
    this->interruptionHandler = interruptionHandler;
  }

protected:
  virtual unsigned int numActiveModules() const = 0;

  void checkInterruption() {
    if (interruptionHandler) {
      if (interruptionHandler()) {
        throw infomap::InterruptException();
      }
    }
  }

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

private:

  interruptionHandlerFn *interruptionHandler = NULL;

};

} /* namespace infomap */

#endif // INFOMAP_OPTIMIZER_BASE_H_
