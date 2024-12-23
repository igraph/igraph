export type NodeBase = {
  id: number;
  flow?: number;
};

type StateNodeBase = {
  stateId: number;
  layerId?: number;
};

export type CluNode<T extends string = "moduleId"> = {
  [key in T]: number;
} & NodeBase;

export type CluStateNode = CluNode & StateNodeBase;

export type Clu = CluNode[];

export type MetaClu = CluNode<"meta">[];

export type TreeNode = {
  path: number[];
  modules?: number[];
  mec?: number;
  name?: string;
} & NodeBase;

export type TreeStateNode = TreeNode & StateNodeBase;

export type Tree = TreeNode[];

export type StateTree = TreeStateNode[];

export type FileTypes = Clu | MetaClu | Tree | StateTree;

const isClu = (file: FileTypes): file is Clu | MetaClu => {
  const clu = file as Clu | MetaClu;
  return clu.length > 0 && ("moduleId" in file[0] || "meta" in file[0]);
};

const isTree = (file: FileTypes): file is Tree | StateTree => {
  const clu = file as Tree | StateTree;
  return clu.length > 0 && "path" in file[0];
};

export default function fileToString(file: FileTypes): string {
  if (file.length === 0) {
    return "";
  }

  if (isClu(file)) {
    return cluToString(file);
  } else if (isTree(file)) {
    return treeToString(file);
  }

  return "" as never;
}

function cluToString(clu: Clu | MetaClu): string {
  let result = "";

  for (let node of clu) {
    result += node.id;
    if ("moduleId" in node) result += ` ${node.moduleId}`;
    else if ("meta" in node) result += ` ${node.meta}`;
    if (node.flow != null) result += ` ${node.flow}`;
    result += "\n";
  }

  return result;
}

function treeToString(tree: Tree | StateTree): string {
  let result = "";

  for (let node of tree) {
    const flow = node.flow ?? 0;
    const name = node.name ?? node.id;

    result += `${node.path.join(":")} ${flow} "${name}"`;
    if ("stateId" in node) result += ` ${node.stateId}`;
    result += ` ${node.id}`;
    if ("layerId" in node) result += ` ${node.layerId}`;
    result += "\n";
  }

  return result;
}
