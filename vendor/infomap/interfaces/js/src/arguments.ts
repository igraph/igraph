type OutputFormats =
  | "clu"
  | "tree"
  | "ftree"
  | "newick"
  | "json"
  | "csv"
  | "network"
  | "states"
  | "flow";

export interface Arguments
  extends Partial<{
    // input
    clusterData: string;
    noInfomap: boolean;
    skipAdjustBipartiteFlow: boolean;
    bipartiteTeleportation: boolean;
    weightThreshold: number;
    noSelfLinks: boolean;
    nodeLimit: number;
    matchableMultilayerIds: number;
    assignToNeighbouringModule: boolean;
    metaData: string;
    metaDataRate: number;
    metaDataUnweighted: boolean;
    // output
    tree: boolean;
    ftree: boolean;
    clu: boolean;
    verbose: 1 | 2 | 3;
    silent: boolean;
    outName: string;
    noFileOutput: boolean;
    cluLevel: number;
    output: OutputFormats | OutputFormats[];
    hideBipartiteNodes: boolean;
    printAllTrials: boolean;
    // algorithm
    twoLevel: boolean;
    flowModel: "undirected" | "directed" | "undirdir" | "outdirdir" | "rawdir";
    directed: boolean;
    recordedTeleportation: boolean;
    useNodeWeightsAsFlow: boolean;
    toNodes: boolean;
    teleportationProbability: number;
    regularized: boolean;
    regularizationStrength: number;
    entropyCorrected: boolean;
    entropyCorrectionStrength: number;
    markovTime: number;
    variableMarkovTime: boolean;
    variableMarkovDamping: number;
    preferredNumberOfModules: number;
    multilayerRelaxRate: number;
    multilayerRelaxLimit: number;
    multilayerRelaxLimitUp: number;
    multilayerRelaxLimitDown: number;
    multilayerRelaxByJsd: boolean;
    // accuracy
    seed: number;
    numTrials: number;
    coreLoopLimit: number;
    coreLevelLimit: number;
    tuneIterationLimit: number;
    coreLoopCodelengthThreshold: number;
    tuneIterationRelativeThreshold: number;
    fastHierarchicalSolution: 1 | 2 | 3;
    preferModularSolution: boolean;
    innerParallelization: boolean;
    // about
    help: boolean | "advanced";
    version: boolean;
  }> { }

export default function argumentsToString(args: Arguments) {
  let result = "";

  if (args.clusterData) result += " --cluster-data " + args.clusterData;

  if (args.noInfomap) result += " --no-infomap";

  if (args.skipAdjustBipartiteFlow) result += " --skip-adjust-bipartite-flow";

  if (args.bipartiteTeleportation) result += " --bipartite-teleportation";

  if (args.weightThreshold != null)
    result += " --weight-threshold " + args.weightThreshold;

  if (args.noSelfLinks) result += " --no-self-links";

  if (args.nodeLimit != null) result += " --node-limit " + args.nodeLimit;

  if (args.matchableMultilayerIds != null)
    result += " --matchable-multilayer-ids " + args.matchableMultilayerIds;

  if (args.assignToNeighbouringModule)
    result += " --assign-to-neighbouring-module";

  if (args.metaData) result += " --meta-data " + args.metaData;

  if (args.metaDataRate != null)
    result += " --meta-data-rate " + args.metaDataRate;

  if (args.metaDataUnweighted) result += " --meta-data-unweighted";

  if (args.tree) result += " --tree";

  if (args.ftree) result += " --ftree";

  if (args.clu) result += " --clu";

  if (args.verbose) result += " -" + "v".repeat(args.verbose);

  if (args.silent) result += " --silent";

  if (args.outName) result += " --out-name " + args.outName;

  if (args.noFileOutput) result += " --no-file-output";

  if (args.cluLevel != null) result += " --clu-level " + args.cluLevel;

  if (args.output != null) {
    if (typeof args.output === "string") result += " -o " + args.output;
    else result += " -o " + args.output.join(",");
  }

  if (args.hideBipartiteNodes) result += " --hide-bipartite-nodes";

  if (args.printAllTrials) result += " --print-all-trials";

  if (args.twoLevel) result += " --two-level";

  if (args.flowModel) result += " --flow-model " + args.flowModel;

  if (args.directed) result += " --directed";

  if (args.recordedTeleportation) result += " --recorded-teleportation";

  if (args.useNodeWeightsAsFlow) result += " --use-node-weights-as-flow";

  if (args.toNodes) result += " --to-nodes";

  if (args.teleportationProbability != null)
    result += " --teleportation-probability " + args.teleportationProbability;

  if (args.regularized) result += " --regularized";

  if (args.regularizationStrength != null)
    result += " --regularization-strength " + args.regularizationStrength;

  if (args.entropyCorrected) result += " --entropy-corrected";

  if (args.entropyCorrectionStrength != null)
    result +=
      " --entropy-correction-strength " + args.entropyCorrectionStrength;

  if (args.markovTime) result += " --markov-time " + args.markovTime;

  if (args.variableMarkovTime) result += " --variable-markov-time";

  if (args.variableMarkovDamping != null)
    result +=
      " --variable-markov-damping " + args.variableMarkovDamping;

  if (args.preferredNumberOfModules != null)
    result += " --preferred-number-of-modules " + args.preferredNumberOfModules;

  if (args.multilayerRelaxRate != null)
    result += " --multilayer-relax-rate " + args.multilayerRelaxRate;

  if (args.multilayerRelaxLimit != null)
    result += " --multilayer-relax-limit " + args.multilayerRelaxLimit;

  if (args.multilayerRelaxLimitUp != null)
    result += " --multilayer-relax-limit-up " + args.multilayerRelaxLimitUp;

  if (args.multilayerRelaxLimitDown != null)
    result += " --multilayer-relax-limit-down " + args.multilayerRelaxLimitDown;

  if (args.multilayerRelaxByJsd != null)
    result += " --multilayer-relax-by-jsd " + args.multilayerRelaxByJsd;

  if (args.seed != null) result += " --seed " + args.seed;

  if (args.numTrials != null) result += " --num-trials " + args.numTrials;

  if (args.coreLoopLimit != null)
    result += " --core-loop-limit " + args.coreLoopLimit;

  if (args.coreLevelLimit != null)
    result += " --core-level-limit " + args.coreLevelLimit;

  if (args.tuneIterationLimit != null)
    result += " --tune-iteration-limit " + args.tuneIterationLimit;

  if (args.coreLoopCodelengthThreshold != null)
    result +=
      " --core-loop-codelength-threshold " + args.coreLoopCodelengthThreshold;

  if (args.tuneIterationRelativeThreshold != null)
    result +=
      " --tune-iteration-relative-threshold " +
      args.tuneIterationRelativeThreshold;

  if (args.fastHierarchicalSolution)
    result += " -" + "F".repeat(args.fastHierarchicalSolution);

  if (args.preferModularSolution) result += " --prefer-modular-solution";

  if (args.help) result += args.help === "advanced" ? " -hh" : " -h";

  if (args.version) result += " --version";

  return result;
}
