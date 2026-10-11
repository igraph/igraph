/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#include "OutputPlan.h"
#include "InfomapError.h"
#include "Network.h"
#include "OutputFormats.h"
#include "SafeFile.h"
#include "../core/InfomapBase.h"
#include "../utils/Console.h"
#include "../utils/convert.h"
#include "../utils/format.h"

#include <algorithm>
#include <stdexcept>
#include <utility>

namespace infomap {

namespace {

  std::string foldSeparators(const std::string& path)
  {
    std::string folded;
    folded.reserve(path.size());
    for (char c : path) {
      const char separator = c == '\\' ? '/' : c;
#ifdef _WIN32
      // Keep a leading "//". On Windows that prefix is the whole difference
      // between a UNC path ("\\server\share\x") and a root-relative one
      // ("\server\share\x"), which name different files; folding them together
      // would refuse a run that overwrites nothing. POSIX has no such
      // distinction -- Linux and macOS resolve "//x" to "/x" -- so there the
      // collapse stays unconditional.
      const bool atUncPrefix = folded.size() == 1 && folded[0] == '/';
#else
      const bool atUncPrefix = false;
#endif
      if (separator == '/' && !folded.empty() && folded.back() == '/' && !atUncPrefix)
        continue;
      folded.push_back(separator);
    }
    return folded;
  }

  bool isAbsolutePath(const std::string& folded)
  {
    if (folded.empty())
      return false;
    if (folded[0] == '/') // POSIX, and a Windows root-relative path
      return true;
#ifdef _WIN32
    // "C:/...". Only on Windows: on POSIX "a:/b" is a relative path into a
    // directory named "a:", and calling it absolute would skip the anchoring.
    const char drive = folded[0];
    const bool isDriveLetter = (drive >= 'A' && drive <= 'Z') || (drive >= 'a' && drive <= 'z');
    return folded.size() >= 3 && isDriveLetter && folded[1] == ':' && folded[2] == '/';
#else
    return false;
#endif
  }

  // Relative paths on both sides of the comparison resolve against the working
  // directory, and the preflight runs before anything can chdir, so reading it
  // once is enough.
  const std::string& workingDirectory()
  {
    static const std::string cwd = currentWorkingDirectory();
    return cwd;
  }

  // Path comparison for "would this output overwrite an input?". std::filesystem
  // is C++17 and this translation unit is C++14; realpath/stat are POSIX-only and
  // the project also builds on Windows. So: fold separators, collapse repeated
  // ones, resolve a relative path against the working directory, and drop "./"
  // segments, then compare. Anchoring to the working directory is what lets the
  // two sides be given in different forms -- a bare `ml.net` input against an
  // absolute output directory names the same file, and every combination of that
  // (absolute input with `.`, absolute --summary-json path) used to slip through.
  // It does not resolve symlinks, "..", Windows case-insensitivity, or the
  // extended-length "\\?\C:\" prefix, so it can miss an exotic aliasing; it never
  // reports a collision that is not one --
  // note that anchoring only ever makes two paths that already named the same
  // file compare equal, since the anchor is the same string for both sides. If
  // the working directory cannot be read at all, the comparison degrades to the
  // paths as given, which is a miss and never a false refusal.
  std::string normalizePathForComparison(const std::string& path)
  {
    std::string folded = foldSeparators(path);
    if (!folded.empty() && !isAbsolutePath(folded)) {
      const std::string& cwd = workingDirectory();
      if (!cwd.empty())
        folded = foldSeparators(cwd + "/" + folded);
    }

    std::string result;
    result.reserve(folded.size());
    for (std::size_t i = 0; i < folded.size();) {
      const bool atSegmentStart = i == 0 || folded[i - 1] == '/';
      if (atSegmentStart && folded.compare(i, 2, "./") == 0) {
        i += 2;
        continue;
      }
      result.push_back(folded[i]);
      ++i;
    }
    return result;
  }

  // The option that owns a report path, so the refusal can name it. The
  // basename guidance does not apply to these: the caller pointed the option at
  // a file directly, and --out-name does not move it.
  std::string reportOptionFlag(const std::string& reportKey)
  {
    if (reportKey == "summary_json")
      return "--summary-json";
    if (reportKey == "timing_json")
      return "--timing-json";
    if (reportKey == "run_manifest")
      return "--manifest-json";
    if (reportKey == "trial_results")
      return "--trial-results";
    return {};
  }

  // Every file the run reads. Writing a result over any of them destroys input.
  std::vector<std::string> inputPaths(const Config& config)
  {
    std::vector<std::string> inputs;
    inputs.push_back(config.networkFile);
    inputs.push_back(config.clusterDataFile);
    inputs.push_back(config.metaDataFile);
    inputs.insert(inputs.end(), config.additionalInput.begin(), config.additionalInput.end());
    return inputs;
  }

  OutputArtifact artifact(const Config& config,
                          const std::string& basename,
                          OutputPhase phase,
                          OutputKind kind,
                          const std::string& resultKey,
                          const std::string& label,
                          bool states = false,
                          bool printFlow = false)
  {
    OutputArtifact output;
    output.resultKey = resultKey;
    output.label = label;
    output.filename = outputFilenameForResultKey(basename, resultKey);
    output.phase = phase;
    output.kind = kind;
    output.states = states;
    output.cluLevel = config.cluLevel;
    output.printFlow = printFlow;
    return output;
  }

  void addPhysicalAndStateArtifacts(std::vector<OutputArtifact>& artifacts,
                                    const Config& config,
                                    const std::string& basename,
                                    OutputKind kind,
                                    const std::string& physicalResultKey,
                                    const std::string& stateResultKey,
                                    const std::string& label,
                                    const std::string& physicalLabel,
                                    const std::string& stateLabel)
  {
    if (!config.printStates()) {
      artifacts.push_back(artifact(config, basename, OutputPhase::AfterPartition, kind, physicalResultKey, label));
      return;
    }

    artifacts.push_back(artifact(config, basename, OutputPhase::AfterPartition, kind, physicalResultKey, physicalLabel));
    artifacts.push_back(artifact(config, basename, OutputPhase::AfterPartition, kind, stateResultKey, stateLabel, true));
  }

  std::pair<std::string, std::string> splitExtension(const std::string& filename)
  {
    const auto slashPos = filename.find_last_of("/\\");
    const auto dotPos = filename.find_last_of('.');
    if (dotPos == std::string::npos || (slashPos != std::string::npos && dotPos < slashPos)) {
      return { filename, "" };
    }
    return { filename.substr(0, dotPos), filename.substr(dotPos + 1) };
  }

  std::vector<std::string> summarizeOutputFiles(const std::vector<std::pair<std::string, std::string>>& outputFiles)
  {
    std::vector<std::pair<std::string, std::vector<std::string>>> groups;
    for (const auto& outputFile : outputFiles) {
      const auto parts = splitExtension(outputFile.second);
      auto it = std::find_if(groups.begin(), groups.end(), [&parts](const std::pair<std::string, std::vector<std::string>>& group) {
        return group.first == parts.first;
      });
      if (it == groups.end()) {
        groups.push_back({ parts.first, { parts.second } });
      } else {
        it->second.push_back(parts.second);
      }
    }

    std::vector<std::string> summaries;
    summaries.reserve(groups.size());
    for (const auto& group : groups) {
      if (group.second.size() == 1) {
        summaries.push_back(group.second.front().empty() ? group.first : fmt::format(FMT_STRING("{}.{}"), group.first, group.second.front()));
        continue;
      }
      summaries.push_back(fmt::format(FMT_STRING("{}.{{{}}}"), group.first, io::stringify(group.second, ",")));
    }
    return summaries;
  }

  void writeModularOutput(InfomapBase& infomap, const OutputArtifact& output)
  {
    switch (output.kind) {
    case OutputKind::Tree:
      infomap.writeTree(output.filename, output.states);
      break;
    case OutputKind::FlowTree:
      infomap.writeFlowTree(output.filename, output.states);
      break;
    case OutputKind::Newick:
      infomap.writeNewickTree(output.filename, output.states);
      break;
    case OutputKind::Json:
      infomap.writeJsonTree(output.filename, output.states, output.writeLinks);
      break;
    case OutputKind::Csv:
      infomap.writeCsvTree(output.filename, output.states);
      break;
    case OutputKind::Clu:
      infomap.writeClu(output.filename, output.states, output.cluLevel);
      break;
    default:
      throw std::logic_error("Output artifact is not a modular result");
    }
  }

  void writeNetworkOutput(Network& network, const OutputArtifact& output)
  {
    switch (output.kind) {
    case OutputKind::StateNetwork:
      network.writeStateNetwork(output.filename);
      break;
    case OutputKind::PajekNetwork:
    case OutputKind::FlowNetwork:
      network.writePajekNetwork(output.filename, output.printFlow);
      break;
    default:
      throw std::logic_error("Output artifact is not a network result");
    }
  }

} // namespace

std::string outputPlanBasename(const Config& config, int trial)
{
  std::string basename = config.outDirectory + config.outName;

  if (config.printAllTrials && trial != -1 && config.numTrials > 1) {
    basename += fmt::format(FMT_STRING("_trial_{}"), trial);
  }

  return basename;
}

std::vector<OutputArtifact> planOutputArtifacts(const Config& config, const std::string& basename, OutputPhase phase)
{
  if (config.noFileOutput) {
    return {};
  }

  std::vector<OutputArtifact> artifacts;

  if (phase == OutputPhase::BeforeFlow) {
    if (config.printStateNetwork) {
      artifacts.push_back(artifact(config, basename, phase, OutputKind::StateNetwork, "states", "state network", true));
    }
    if (config.printPajekNetwork) {
      artifacts.push_back(config.printStates()
                              ? artifact(config, basename, phase, OutputKind::PajekNetwork, "states_as_physical", "state network as Pajek")
                              : artifact(config, basename, phase, OutputKind::PajekNetwork, "net", "Pajek network"));
    }
    return artifacts;
  }

  if (phase == OutputPhase::AfterFlow) {
    if (config.printFlowNetwork) {
      artifacts.push_back(config.printStates()
                              ? artifact(config, basename, phase, OutputKind::FlowNetwork, "flow_as_physical", "flow state network as Pajek", true, true)
                              : artifact(config, basename, phase, OutputKind::FlowNetwork, "flow", "flow network", false, true));
    }
    return artifacts;
  }

  if (phase != OutputPhase::AfterPartition) {
    return artifacts;
  }

  if (config.printTree) {
    addPhysicalAndStateArtifacts(artifacts, config, basename, OutputKind::Tree, "tree", "tree_states", "tree", "physical tree", "state tree");
  }
  if (config.printFlowTree) {
    addPhysicalAndStateArtifacts(artifacts, config, basename, OutputKind::FlowTree, "ftree", "ftree_states", "flow tree", "physical flow tree", "state flow tree");
  }
  if (config.printNewick) {
    addPhysicalAndStateArtifacts(artifacts, config, basename, OutputKind::Newick, "newick", "newick_states", "Newick tree", "physical Newick tree", "state Newick tree");
  }
  if (config.printJson) {
    addPhysicalAndStateArtifacts(artifacts, config, basename, OutputKind::Json, "json", "json_states", "JSON tree", "physical JSON tree", "state JSON tree");
  }
  if (config.printCsv) {
    addPhysicalAndStateArtifacts(artifacts, config, basename, OutputKind::Csv, "csv", "csv_states", "CSV tree", "physical CSV tree", "state CSV tree");
  }
  if (config.printClu) {
    addPhysicalAndStateArtifacts(artifacts, config, basename, OutputKind::Clu, "clu", "clu_states", "node modules", "physical node modules", "state node modules");
  }

  return artifacts;
}

std::vector<OutputArtifact> planOutputArtifacts(const Config& config, OutputPhase phase, int trial)
{
  return planOutputArtifacts(config, outputPlanBasename(config, trial), phase);
}

std::vector<std::pair<std::string, std::string>> planReportArtifacts(const Config& config)
{
  std::vector<std::pair<std::string, std::string>> reports;
  const auto add = [&](const std::string& key, const std::string& path) {
    if (!path.empty() && path != "-")
      reports.emplace_back(key, path);
  };
  add("summary_json", config.summaryJsonPath);
  add("timing_json", config.timingJsonPath);
  add("run_manifest", config.runManifestPath);
  // The shard results file is an output the run writes, so it belongs here for
  // both consumers: the preflight collision check, and the manifest's list of
  // produced artifacts. Its writer already honours the overwrite policy, so
  // listing it only moves that refusal ahead of the work.
  add("trial_results", config.trialResultsPath);
  return reports;
}

std::vector<std::string> planAllOutputPaths(const Config& config, HigherOrderInput higherOrder)
{
  std::vector<std::string> paths;

  // Config::stateOutput is still false when the pre-flight runs -- it is set by
  // configureNetworkMode() once the network is read -- so the classification is
  // supplied instead. Every write phase now happens after that call, so one value
  // covers them all. Getting it wrong in either direction is a live defect: too
  // few paths lets a run destroy its own input (#1018), too many refuses a run
  // that would have been fine.
  //
  // The copy preserves a stateOutput the caller set itself -- the field is public
  // and the Python and R bindings expose a setter -- and the classification only
  // ever turns it on, never off. So the plan agrees with the writers on that case
  // too, which is what the forced first-order plan this replaced got wrong.
  Config planConfig = config;
  if (higherOrder == HigherOrderInput::Yes)
    planConfig.setStateOutput();

  const auto collectPhase = [&](OutputPhase phase, int trial) {
    for (const auto& artifact : planOutputArtifacts(planConfig, phase, trial))
      paths.push_back(artifact.filename);
  };

  // Network sidecars are written once, independent of trial.
  collectPhase(OutputPhase::BeforeFlow, -1);
  collectPhase(OutputPhase::AfterFlow, -1);

  // The final modular result is written once with the canonical basename, and
  // additionally per trial when --print-all-trials uses separate files.
  //
  // The trial numbers are the writer's global ones -- `trialOffset + i + 1`, the
  // same expression both writeResult call sites use -- not a local 1..numTrials.
  // A shard planned `_trial_1.._trial_4` while writing `_trial_11.._trial_14`,
  // which left both consumers of this plan looking at files that never exist:
  // --no-overwrite could only discover a real collision by hitting it mid-run,
  // after earlier artifacts were already on disk, and the input-overwrite check
  // above could not see that a run was about to destroy an input named like one
  // of its own per-trial outputs.
  collectPhase(OutputPhase::AfterPartition, -1);
  if (config.printAllTrials && config.numTrials > 1) {
    for (unsigned int i = 0; i < config.numTrials; ++i)
      collectPhase(OutputPhase::AfterPartition, static_cast<int>(config.trialOffset + i + 1));
  }

  for (const auto& report : planReportArtifacts(config))
    paths.push_back(report.second);

  return paths;
}

void preflightOutputTargets(const Config& config, HigherOrderInput higherOrder)
{
  const auto plannedPaths = planAllOutputPaths(config, higherOrder);

  // Checked before the overwrite policy, and deliberately not subject to it:
  // writing a result over the run's own input destroys the input, and
  // --no-overwrite cannot be the mitigation because it rejects every ordinary
  // re-run as well. Reachable from a plain command whenever the output goes to
  // the input's own directory and a requested format shares its extension --
  // `Infomap net.json . -o json` replaced the network with the result tree, and
  // `Infomap ml.net . -o network` replaced a multilayer input with its
  // single-layer projection, both exiting 0.
  // planAllOutputPaths mixes two kinds of target: artifacts named from
  // outDirectory + outName, and report paths the caller gave directly. The
  // mitigation differs, so the message has to know which one collided --
  // suggesting --out-name for a --summary-json path would be useless advice.
  //
  const auto reports = planReportArtifacts(config);
  const auto reportOwnerOf = [&reports](const std::string& path) {
    for (const auto& report : reports) {
      if (report.second == path)
        return reportOptionFlag(report.first);
    }
    return std::string();
  };

  const auto inputs = inputPaths(config);
  for (const auto& path : plannedPaths) {
    const auto normalizedOutput = normalizePathForComparison(path);
    for (const auto& input : inputs) {
      if (input.empty())
        continue;
      if (normalizePathForComparison(input) != normalizedOutput)
        continue;

      const auto ownerFlag = reportOwnerOf(path);
      const auto guidance = ownerFlag.empty()
          ? std::string("Write to a different directory, or pass --out-name to give the result a different basename.")
          : fmt::format(FMT_STRING("Point {} at a different file."), ownerFlag);

      throw InfomapError(ExitCode::OutputError,
                         fmt::format(FMT_STRING("Refusing to write output '{}' over the input file '{}'. {}"),
                                     path,
                                     input,
                                     guidance));
    }
  }

  if (config.overwriteOutput())
    return;

  for (const auto& path : plannedPaths) {
    if (pathExists(path))
      throw InfomapError(ExitCode::OutputError, fmt::format(FMT_STRING("Output file already exists: '{}'"), path));
  }
}

void writeOutputArtifact(InfomapBase& infomap, Network& network, const OutputArtifact& output)
{
  if (output.phase == OutputPhase::AfterPartition) {
    writeModularOutput(infomap, output);
    return;
  }

  writeNetworkOutput(network, output);
}

void writeOutputArtifacts(InfomapBase& infomap, Network& network, OutputPhase phase, int trial)
{
  const auto artifacts = planOutputArtifacts(infomap, phase, trial);
  std::vector<std::pair<std::string, std::string>> prettyOutputFiles;

  for (const auto& output : artifacts) {
    writeOutputArtifact(infomap, network, output);
    if (phase != OutputPhase::AfterPartition) {
      Console().status("Output", fmt::format(FMT_STRING("{} -> {}"), output.label, output.filename));
      continue;
    }
    prettyOutputFiles.emplace_back(output.label, output.filename);
  }

  if (prettyOutputFiles.empty()) {
    return;
  }

  if (prettyOutputFiles.size() == 1) {
    Console().status("Output", fmt::format(FMT_STRING("{} -> {}"), prettyOutputFiles.front().first, prettyOutputFiles.front().second));
    return;
  }

  const auto summaries = summarizeOutputFiles(prettyOutputFiles);
  for (unsigned int i = 0; i < summaries.size(); ++i) {
    const std::string prefix = i == 0 ? fmt::format(FMT_STRING("{} files -> "), prettyOutputFiles.size()) : "         ";
    Console().status("Output", fmt::format(FMT_STRING("{}{}"), prefix, summaries[i]));
  }
}

} // namespace infomap
