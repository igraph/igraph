import type {
  CluNode,
  CluStateNode,
  NodeBase,
  TreeNode as JsonTreeNode,
  TreeStateNode as JsonTreeStateNode,
} from "@mapequation/infomap/filetypes";
import type { Header as JsonHeader, Module } from "@mapequation/infomap";

type Optional<T, K extends keyof T> = Pick<Partial<T>, K> & Omit<T, K>;

type Header = Optional<JsonHeader, "directed"> & { cluLevel?: number };

// Support stree files
type Path = { path: string | number[] };
type TreeNode = Omit<JsonTreeNode, "modules" | "mec" | "path"> & Path;
type TreeStateNode = Omit<JsonTreeStateNode, "modules" | "mec" | "path"> & Path;

export type Result<NodeType extends NodeBase> = Header & {
  nodes: NodeType[];
  modules?: Module[];
};

export type Extension = "clu" | "tree";

export function extension(filename: string) {
  return filename.split(".").pop() ?? null;
}

export function lines(file: string) {
  // split file on \r\n or \n
  return file.split(/\r?\n/);
}

export function readFile(file: File): Promise<string> {
  return new Promise((resolve, reject) => {
    const reader = new FileReader();
    reader.onload = () => resolve(reader.result as string);
    reader.onerror = reject;
    reader.readAsText(file);
  });
}

export function parse(
  file: string | string[],
  type?: Extension,
  parseLinks = false,
  strictHeader = true
) {
  if (typeof file === "string") {
    file = lines(file);
  }

  const nodeHeader = parseNodeHeader(file);

  if ((type && type === "clu") || nodeHeader.includes("module")) {
    if (nodeHeader[0] === "node_id") {
      return parseClu(file, nodeHeader, strictHeader);
    } else if (nodeHeader[0] === "state_id") {
      return parseClu<CluStateNode>(file, nodeHeader, strictHeader);
    }
  } else if ((type && type === "tree") || nodeHeader.includes("path")) {
    if (nodeHeader[nodeHeader.length - 2] === "name") {
      return parseTree(file, nodeHeader, parseLinks, strictHeader);
    } else if (
      nodeHeader[nodeHeader.length - 2] === "state_id" ||
      nodeHeader[nodeHeader.length - 1] === "layer_id"
    ) {
      return parseTree<TreeStateNode>(
        file,
        nodeHeader,
        parseLinks,
        strictHeader
      );
    }
  }
  throw new Error("File must be either a tree or a clu file");
}

const map = {
  node_id: "id",
  module: "moduleId",
  flow: "flow",
  name: "name",
  state_id: "stateId",
  layer_id: "layerId",
  path: "path",
} as const;

export function parseClu<NodeType extends CluNode>(
  file: string | string[],
  nodeHeader?: string[],
  strictHeader = true
): Result<NodeType> {
  // First order
  // # node_id module flow
  // 10 1 0.0384615
  // 11 1 0.0384615
  // 12 1 0.0384615

  // States
  // # state_id module flow node_id
  // 1 1 0.166667 1
  // 2 1 0.166667 2
  // 3 1 0.166667 3

  // Multilayer
  // # state_id module flow node_id layer_id
  // 3 1 0.166667 1 2
  // 4 1 0.166667 2 2
  // 5 1 0.166667 3 2

  if (typeof file === "string") {
    file = lines(file);
  }

  const header = parseHeader(file, strictHeader);

  if (!nodeHeader) {
    nodeHeader = parseNodeHeader(file);
  }

  const nodes: NodeType[] = [];

  for (const [i, line] of nodeSection(file)) {
    const fields = line.split(" ").map(Number);

    if (fields.length < nodeHeader.length) {
      continue;
    }

    const node: Partial<NodeType> = {};

    for (let i = 0; i < nodeHeader.length; ++i) {
      const field = nodeHeader[i] as keyof typeof map;
      const key = map[field] as keyof NodeType;
      // @ts-ignore
      node[key] = fields[i];
    }

    if (
      node.id === undefined ||
      node.moduleId === undefined ||
      node.flow === undefined
    ) {
      throw new Error(`Invalid node at line ${i + 1}: ${line}`);
    }

    nodes.push(node as NodeType);
  }

  return {
    ...header,
    nodes,
  };
}

export function parseTree<NodeType extends TreeNode>(
  file: string | string[],
  nodeHeader?: string[],
  parseLinks = false,
  strictHeader = true
): Result<NodeType> {
  // First order
  // # path flow name node_id
  // 1:1 0.166667 "i" 1
  // 1:2 0.166667 "j" 2
  // 1:3 0.166667 "k" 3

  // States
  // # path flow name state_id node_id
  // 1:1 0.166667 "i" 1 1
  // 1:2 0.166667 "j" 2 2
  // 1:3 0.166667 "k" 3 3

  // Multilayer
  // # path flow name state_id node_id layer_id
  // 1:1 0.166667 "i" 3 1 2
  // 1:2 0.166667 "j" 4 2 2
  // 1:3 0.166667 "k" 5 3 2

  if (typeof file === "string") {
    file = lines(file);
  }

  const header = parseHeader(file, strictHeader);

  if (!nodeHeader) {
    nodeHeader = parseNodeHeader(file);
  }

  const nodes: NodeType[] = [];

  let lineNo = 0;

  for (const [i, line] of nodeSection(file)) {
    lineNo = i;

    const match = line.match(/[^\s"']+|"([^"]*)"/g);

    if (match === null || match.length < nodeHeader.length) {
      continue;
    }

    const node: Partial<NodeType> = {};

    for (let j = 0; j < nodeHeader.length; ++j) {
      const field = nodeHeader[j] as keyof typeof map;
      const key = map[field] as keyof NodeType;

      switch (field) {
        case "path":
          node.path = match[j]; //.split(":").map(Number);
          break;
        case "name":
          node.name = match[j].slice(1, -1);
          break;
        case "node_id": // FIXME can this be merged with the next case?
          node.id = Number(match[j]);
          break;
        case "state_id":
        case "layer_id":
        case "flow":
          // @ts-ignore
          node[key] = Number(match[j]);
          break;
        default:
          console.warn(`Unknown field ${field}`);
          break;
      }
    }

    if (
      node.path === undefined ||
      node.flow === undefined ||
      node.name === undefined ||
      node.id === undefined
    ) {
      throw new Error(`Invalid node at line ${i + 1}: ${line}`);
    }

    nodes.push(node as NodeType);
  }

  // lineNo is the index of the last line of the node section
  ++lineNo;

  if (file.length <= lineNo + 1) {
    return {
      ...header,
      nodes,
    };
  }

  const modules = [];
  let module: Module | null = null;

  const linkHeader = file[lineNo];

  let directed = false;

  if (linkHeader?.startsWith("*Links")) {
    const parts = linkHeader.trim().split(" ");
    directed = parts.length > 1 && parts[1] === "directed";
  }

  for (const [, line] of linkSection(file, lineNo)) {
    if (line.startsWith("*Links")) {
      const [, moduleId, ...rest] = line.split(" ");

      const path = moduleId === "root" ? [0] : moduleId.split(":").map(Number);
      const [enterFlow, exitFlow, numEdges, numChildren] = rest.map(Number);

      module = {
        path,
        enterFlow,
        exitFlow,
        numEdges,
        numChildren,
        links: parseLinks ? [] : undefined,
      };

      modules.push(module);
      continue;
    }

    if (!parseLinks || !module) {
      continue;
    }

    const [source, target, flow] = line.split(" ").map(Number);

    module.links?.push({ source, target, flow });
  }

  return {
    ...header,
    directed,
    nodes,
    modules,
  };
}

function* nodeSection(lines: string[]) {
  for (let i = 0; i < lines.length; ++i) {
    const line = lines[i];

    if (line.startsWith("#") || line.length === 0) {
      continue;
    }

    if (line.startsWith("*")) {
      break;
    }

    yield [i, line] as [number, string];
  }
}

function* linkSection(lines: string[], start: number = 0) {
  // links section start at lines[start]
  if (lines.length <= start || !lines[start].startsWith("*Links")) {
    return;
  }

  // skip the section header
  start++;

  for (let i = start; i < lines.length; ++i) {
    const line = lines[i];

    if (line.length === 0) {
      continue;
    }

    yield [i, line] as [number, string];
  }
}

function parseNodeHeader(file: string | string[]): string[] {
  // First order tree
  // # path flow name node_id
  // 1:2 0.166667 "j" 2

  // States tree
  // # path flow name state_id node_id
  // 1:2 0.166667 "j" 1 2

  // Multilayer tree states
  // # path flow name state_id node_id layer_id
  // 1:2 0.166667 "j" 1 2 3

  // First order clu
  // # node_id module flow
  // 1 1 0.166667

  // States clu
  // # state_id module flow node_id
  // 1 1 0.166667 2

  // Multilayer clu states
  // # state_id module flow node_id layer_id
  // 1 1 0.166667 2 3

  const prefix = /^#\s*(path|node_id|state_id)/;

  if (typeof file === "string") {
    file = lines(file);
  }

  // Look for the commented line with the header
  for (const line of file) {
    if (!line.startsWith("#")) {
      break;
    }

    if (prefix.test(line)) {
      return line.slice(1).trim().split(" ");
    }
  }

  function isNumeric(arr: string[]) {
    return arr.every((x) => !isNaN(Number(x)));
  }

  function bail(line: string): never {
    throw new Error(`Invalid node line: ${line}`);
  }

  // Try parsing a single node line
  for (const line of file) {
    if (line.startsWith("#")) {
      continue;
    }

    // Only tree files include names
    if (line.includes('"')) {
      // Don't split inside quotes
      const match = line.match(/[^\s"']+|"([^"]*)"/g);
      if (!match) bail(line);

      const pathIndex = 0;
      const nameIndex = 2;
      const numericFields = match.slice();
      numericFields.splice(nameIndex, 1);
      numericFields.splice(pathIndex, 1);

      const path = match[pathIndex].split(":");

      if (!(isNumeric(numericFields) && isNumeric(path))) {
        bail(line);
      }

      switch (match.length) {
        case 4:
          return ["path", "flow", "name", "node_id"];
        case 5:
          return ["path", "flow", "name", "state_id", "node_id"];
        case 6:
          return ["path", "flow", "name", "state_id", "node_id", "layer_id"];
        default:
          bail(line);
      }
    } else {
      const fields = line.split(" ");

      if (!isNumeric(fields)) bail(line);

      switch (fields.length) {
        case 3:
          return ["node_id", "module", "flow"];
        case 4:
          return ["state_id", "module", "flow", "node_id"];
        case 5:
          return ["state_id", "module", "flow", "node_id", "layer_id"];
        default:
          bail(line);
      }
    }

    break;
  }

  throw new Error("Could not parse node header");
}

export function parseHeader(file: string | string[], strict = true): Header {
  // # v1.7.3
  // # ./Infomap examples/networks/ninetriangles.net .
  // # started at 2021-11-02 13:56:33
  // # completed in 0.011717 s
  // # partitioned into 3 levels with 7 top modules
  // # codelength 3.49842 bits
  // # relative codelength savings 26.2781%
  // # flow model undirected
  // # bipartite start id 123
  // # path flow name node_id

  const result: Partial<Header> = {};

  if (typeof file === "string") {
    file = lines(file);
  }

  if (file.length < 8 && strict) {
    throw new Error("Expected file header to have at least 8 lines");
  }

  const version = file[0] && file[0].match(/^# (v\d+\.\d+\.\d+)/)?.[1];

  if (version) {
    result.version = version;
  } else if (strict) {
    throw new Error("Could not parse version");
  }

  const args = file[1] && file[1].match(/^# (.+)/)?.[1];

  if (args) {
    result.args = args;
  } else if (strict) {
    throw new Error("Could not parse args");
  }

  for (let i = 2; i < file.length; ++i) {
    const line = file[i];

    if (!line.startsWith("#")) {
      break;
    }

    const startedAt = line.match(/^# started at (.+)/)?.[1];
    if (startedAt) {
      result.startedAt = startedAt;
      continue;
    }

    const completedIn = line.match(/^# completed in (.+) s/)?.[1];
    if (completedIn) {
      result.completedIn = Number(completedIn);
      continue;
    }

    const partitionedInto = line.match(
      /^# partitioned into (\d+) levels with (\d+) top modules/
    );
    if (partitionedInto) {
      result.numLevels = Number(partitionedInto[1]);
      result.numTopModules = Number(partitionedInto[2]);
      continue;
    }

    const codelength = line.match(/^# codelength (.+) bits/)?.[1];
    if (codelength) {
      result.codelength = Number(codelength);
      continue;
    }

    const relativeCodelengthSavings = line.match(
      /^# relative codelength savings (.+)%/
    )?.[1];
    if (relativeCodelengthSavings) {
      result.relativeCodelengthSavings =
        Number(relativeCodelengthSavings) / 100;
      continue;
    }

    const flowModel = line.match(/^# flow model (.+)/)?.[1];
    if (flowModel) {
      result.flowModel = flowModel;
      continue;
    }

    const cluLevel = line.match(/^# module level (\d+)/)?.[1];
    if (cluLevel) {
      result.cluLevel = Number(cluLevel);
      continue;
    }

    const higherOrder = line.match(/^# higher order/)?.[1];
    if (higherOrder) {
      result.higherOrder = true;
      continue;
    }

    const stateLevel = line.match(/^# (state|physical) level/)?.[1];
    if (stateLevel) {
      result.stateLevel = stateLevel === "state";
      continue;
    }

    const bipartiteStartId = line.match(/^# bipartite start id (\d+)/)?.[1];
    if (bipartiteStartId) {
      result.bipartiteStartId = Number(bipartiteStartId);
      //continue;
    }
  }

  if (!result.higherOrder) {
    result.higherOrder = false;
  }

  if (
    strict &&
    !(
      result.startedAt &&
      result.completedIn &&
      result.numLevels &&
      result.numTopModules &&
      result.codelength &&
      result.relativeCodelengthSavings
    )
  ) {
    throw new Error("Could not parse file header");
  }

  return result as Header;
}
