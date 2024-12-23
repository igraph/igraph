export interface Node {
  id: number;
  name?: string;
  weight?: number;
}

export interface StateNode extends Omit<Node, "weight"> {
  stateId: number;
}

export interface Link {
  source: number;
  target: number;
  weight?: number;
}

export interface Network<LinkType = Link> {
  nodes?: Node[];
  links: LinkType[];
}

export interface BipartiteNetwork extends Network {
  bipartiteStartId: number;
}

export interface StateNetwork extends Required<Network> {
  states: StateNode[];
}

export type BipartiteStateNetwork = Required<BipartiteNetwork> & StateNetwork;

export interface MultilayerLink extends Link {
  sourceLayer: number;
  targetLayer: number;
}

export interface IntraLink extends Link {
  layerId: number;
}

export interface InterLink extends Omit<MultilayerLink, "source" | "target"> {
  id: number;
}

export type MultilayerNetwork = Network<MultilayerLink>;

export interface MultilayerIntraInterNetwork extends Omit<Network, "links"> {
  intra: IntraLink[];
  inter?: InterLink[];
}

export type NetworkTypes =
  | Network
  | BipartiteNetwork
  | StateNetwork
  | BipartiteStateNetwork
  | MultilayerNetwork
  | MultilayerIntraInterNetwork;

const isBipartite = (network: NetworkTypes): network is BipartiteNetwork => {
  const net = network as BipartiteNetwork;
  return "bipartiteStartId" in net && !("states" in net);
};

const isState = (network: NetworkTypes): network is StateNetwork => {
  const net = network as StateNetwork;
  return "states" in net && !("bipartiteStartId" in net);
};

const isBipartiteState = (
  network: NetworkTypes
): network is BipartiteStateNetwork => {
  const net = network as BipartiteStateNetwork;
  return "bipartiteStartId" in net && "states" in net;
};

const isMultilayer = (network: NetworkTypes): network is MultilayerNetwork => {
  const net = network as MultilayerNetwork;
  return (
    "links" in net &&
    net.links.length > 0 &&
    "sourceLayer" in net.links[0] &&
    "targetLayer" in net.links[0]
  );
};

const isMultilayerIntraInter = (
  network: NetworkTypes
): network is MultilayerIntraInterNetwork => {
  const net = network as MultilayerIntraInterNetwork;
  return "intra" in net;
};

export default function toString(network: NetworkTypes) {
  if (isMultilayerIntraInter(network)) {
    return multilayerIntraInterToString(network);
  } else if (isMultilayer(network)) {
    return multilayerToString(network);
  } else if (isBipartiteState(network)) {
    return bipartiteStateToString(network);
  } else if (isBipartite(network)) {
    return bipartiteToString(network);
  } else if (isState(network)) {
    return stateToString(network);
  } else if ("links" in network) {
    return networkToString(network);
  }
}

function multilayerIntraInterToString(network: MultilayerIntraInterNetwork) {
  let result = "";

  if (network.nodes != null) {
    result += nodesToString(network.nodes);
  }

  if (network.intra.length > 0) {
    result += "*Intra\n";

    for (let { layerId, source, target, weight } of network.intra) {
      result += `${layerId} ${source} ${target}`;
      if (weight != null) result += ` ${weight}`;
      result += "\n";
    }
  }

  if (network.inter != null && network.inter.length > 0) {
    result += "*Inter\n";

    for (let { sourceLayer, id, targetLayer, weight } of network.inter) {
      result += `${sourceLayer} ${id} ${targetLayer}`;
      if (weight != null) result += ` ${weight}`;
      result += "\n";
    }
  }

  return result;
}

function multilayerToString(network: MultilayerNetwork) {
  let result = "";

  if (network.nodes != null) {
    result += nodesToString(network.nodes);
  }

  if (network.links.length > 0) {
    result += "*Multilayer\n";

    for (let {
      sourceLayer,
      source,
      targetLayer,
      target,
      weight,
    } of network.links) {
      result += `${sourceLayer} ${source} ${targetLayer} ${target}`;
      if (weight != null) result += ` ${weight}`;
      result += "\n";
    }
  }

  return result;
}

function bipartiteStateToString(network: BipartiteStateNetwork) {
  let result = "";

  result += nodesToString(network.nodes);

  result += statesToString(network.states);

  result += linksToString(
    network.links,
    `*Bipartite ${network.bipartiteStartId}\n`
  );

  return result;
}

function bipartiteToString(network: BipartiteNetwork) {
  let result = "";

  if (network.nodes != null) {
    result += nodesToString(network.nodes);
  }

  result += linksToString(
    network.links,
    `*Bipartite ${network.bipartiteStartId}\n`
  );

  return result;
}

function stateToString(network: StateNetwork) {
  let result = "";

  result += nodesToString(network.nodes);

  result += statesToString(network.states);

  result += linksToString(network.links);

  return result;
}

function networkToString(network: Network) {
  let result = "";

  if (network.nodes != null) {
    result += nodesToString(network.nodes);
  }

  result += linksToString(network.links);

  return result;
}

function nodesToString(nodes: Node[]) {
  if (nodes.length === 0) {
    return "";
  }

  let result = "*Vertices\n";

  for (let { id, name, weight } of nodes) {
    result += id;
    if (name != null) result += ` "${name}"`;
    else result += ` "${id}"`;
    if (weight != null) result += ` ${weight}`;
    result += "\n";
  }

  return result;
}

function statesToString(nodes: StateNode[]) {
  if (nodes.length === 0) {
    return "";
  }

  let result = "*States\n";

  for (let { stateId, id, name } of nodes) {
    result += `${stateId} ${id}`;
    if (name != null) result += ` "${name}"`;
    result += "\n";
  }

  return result;
}

function linksToString(links: Link[], header = "*Links\n") {
  if (links.length === 0) {
    return "";
  }

  let result = header;

  for (let { source, target, weight } of links) {
    result += `${source} ${target}`;
    if (weight != null) result += ` ${weight}`;
    result += "\n";
  }

  return result;
}
