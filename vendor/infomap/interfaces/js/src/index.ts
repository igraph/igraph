import argumentsToString, { Arguments } from "./arguments";
import fileToString, {
  FileTypes,
  TreeNode as Node,
  TreeStateNode as StateNode,
} from "./filetypes";
import networkToString, { NetworkTypes } from "./network";
// @ts-ignore
import InfomapWorker from "./worker/infomap.worker.js";
// @ts-ignore
import MemFile from "./worker/infomap.worker.js.mem";

export interface Changelog {
  body: string | null;
  date: string;
  footer: string | null;
  header: string;
  mentions: string[];
  merge: string | null;
  notes: string[];
  references: string[];
  revert: string | null;
  scope: string | null;
  subject: string;
  type: string | null;
}

export interface Parameter<Required = false> {
  long: string;
  short: string;
  description: string;
  group: string;
  required: Required;
  advanced: boolean;
  incremental: boolean;
  default: boolean | string | number;
}

export interface RequiredParameter extends Parameter<true> {
  longType: string;
  shortType: string;
  default: string;
}

// @ts-ignore
const changelog: Changelog[] = CHANGELOG;
// @ts-ignore
const parameters: (Parameter | RequiredParameter)[] = PARAMETERS;

export type Module = {
  path: number[];
  enterFlow: number;
  exitFlow: number;
  numEdges: number;
  numChildren: number;
  codelength: number;
  links?: {
    source: number;
    target: number;
    flow: number;
  }[];
};

export interface Header {
  version: string;
  args: string;
  startedAt: string;
  completedIn: number;
  codelength: number;
  numLevels: number;
  numTopModules: number;
  relativeCodelengthSavings: number;
  directed: boolean;
  flowModel: string;
  higherOrder: boolean;
  stateLevel?: boolean;
  bipartiteStartId?: number;
}

export type Tree<NodeType = Required<Node>> = Header & {
  nodes: NodeType[];
  modules: Module[];
};

export type StateTree = Tree<Required<StateNode>>;

export interface Result {
  clu?: string;
  clu_states?: string;
  tree?: string;
  tree_states?: string;
  ftree?: string;
  ftree_states?: string;
  newick?: string;
  newick_states?: string;
  json?: Tree;
  json_states?: StateTree;
  csv?: string;
  csv_states?: string;
  net?: string;
  states_as_physical?: string;
  states?: string;
  flow?: string;
  flow_as_physical?: string;
}

export interface EventCallbacks {
  data?: (output: string, id: number) => void;
  progress?: (progress: number, id: number) => void;
  error?: (message: string, id: number) => void;
  finished?: (result: Result, id: number) => void;
}

interface Event<Type extends keyof EventCallbacks> {
  type: Type;
  content: Parameters<Required<EventCallbacks>[Type]>[0];
}

type EventData =
  | Event<"data">
  | Event<"progress">
  | Event<"error">
  | Event<"finished">;

const workerUrl = URL.createObjectURL(
  new Blob([InfomapWorker], { type: "application/javascript" })
);

class Infomap {
  // @ts-ignore
  static __version__: string = VERSION;

  protected events: EventCallbacks = {};
  protected workerId = 0;
  protected workers: { [id: number]: Worker } = {};

  run(...args: Parameters<Infomap["createWorker"]>) {
    const id = this.createWorker(...args);
    this.setHandlers(id);
    return id;
  }

  runAsync(...args: Parameters<Infomap["createWorker"]>) {
    const id = this.createWorker(...args);
    return new Promise<Result>((finished, error) =>
      this.setHandlers(id, { finished, error })
    );
  }

  on<E extends keyof EventCallbacks>(event: E, callback: EventCallbacks[E]) {
    this.events[event] = callback;
    return this;
  }

  protected createWorker({
    network,
    filename,
    args,
    files,
  }: {
    network?: string | NetworkTypes;
    filename?: string;
    args?: string | Arguments;
    files?: { [filename: string]: string | FileTypes };
  }) {
    network = network ?? "";
    filename = filename ?? "network.net";
    args = args ?? "";
    files = files ?? {};

    if (typeof network !== "string") {
      network = networkToString(network);
    }

    if (typeof args !== "string") {
      args = argumentsToString(args);
    }

    for (let key of Object.keys(files)) {
      const file = files[key];
      if (typeof file !== "string") {
        files[key] = fileToString(file);
      }
    }

    const index = filename.lastIndexOf(".");
    const networkName = index > 0 ? filename.slice(0, index) : filename;
    const outNameMatch = args.match(/--out-name\s(\S+)/);
    const outName =
      outNameMatch && outNameMatch[1] ? outNameMatch[1] : networkName;

    const worker = new Worker(workerUrl);
    const id = this.workerId++;
    this.workers[id] = worker;

    worker.postMessage({
      memBuffer: new Uint8Array(MemFile),
      arguments: args.split(" "),
      filename,
      network,
      outName,
      files,
    });

    return id;
  }

  protected setHandlers(id: number, events = this.events) {
    const worker = this.workers[id];
    const { data, progress, error, finished } = { ...this.events, ...events };

    worker.onmessage = (event: MessageEvent<EventData>) => {
      if (event.data.type === "data") {
        if (data) {
          data(event.data.content, id);
        }
        if (progress) {
          const match = event.data.content.match(/^Trial (\d+)\/(\d+)/);
          if (match) {
            const trial = Number(match[1]);
            const totTrials = Number(match[2]);
            const currentProgress = (100 * trial) / (totTrials + 1);
            progress(currentProgress, id);
          } else {
            const summary = event.data.content.match(/^Summary after/);
            if (summary) {
              progress(100, id);
            }
          }
        }
      } else if (error && event.data.type === "error") {
        this._terminate(id);
        error(event.data.content, id);
      } else if (finished && event.data.type === "finished") {
        this._terminate(id);
        finished(event.data.content, id);
      }
    };

    worker.onerror = (err: ErrorEvent) => {
      err.preventDefault();
      if (error) error(err.message, id);
    };
  }

  protected _terminate(id: number, timeout = 1000): Promise<boolean> {
    if (!this.workers[id]) return Promise.resolve(false);

    const worker = this.workers[id];

    return new Promise<boolean>((resolve) => {
      const terminate = () => {
        worker.terminate();
        resolve(true);
      };

      if (timeout <= 0) {
        terminate();
      } else {
        setTimeout(terminate, timeout);
      }

      delete this.workers[id];
    });
  }

  // Returns true if the worker with `id` existed and was terminated.
  async terminate(id: number, timeout = 1000) {
    return await this._terminate(id, timeout);
  }
}

export {
  Infomap as default,
  Infomap,
  changelog,
  parameters,
  networkToString,
  argumentsToString,
};
