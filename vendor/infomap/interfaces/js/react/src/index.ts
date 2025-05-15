import Infomap from "@mapequation/infomap";
import type { Arguments } from "@mapequation/infomap/arguments";
import { useCallback, useDebugValue, useState } from "react";

export function useInfomap(args?: Arguments) {
  const [progress, setProgress] = useState(0);
  const [running, setRunning] = useState(false);
  const [infomap] = useState(() => new Infomap().on("progress", setProgress));

  const run = useCallback(
    (...params: Parameters<Infomap["run"]>) => {
      setRunning(true);
      setProgress(0);
      applyArgs(params, args);
      infomap.run(...params);
      setRunning(false);
    },
    [infomap, args]
  );

  const runAsync = useCallback(
    (...params: Parameters<Infomap["runAsync"]>) => {
      setRunning(true);
      setProgress(0);
      applyArgs(params, args);
      return infomap
        .runAsync(...params)
        .then((result) => {
          setRunning(false);
          return result;
        })
        .catch((e) => {
          setRunning(false);
          throw e;
        });
    },
    [infomap, args]
  );

  useDebugValue(running ? "Running" : "Stopped");

  return {
    infomap,
    progress,
    running,
    run,
    runAsync,
    on(...params: Parameters<Infomap["on"]>) {
      infomap.on(...params);
      return this;
    },
  };
}

function applyArgs(params: Parameters<Infomap["run"]>, args?: Arguments) {
  if (!params.length) params[0] = {};
  const param = params[0];

  if (typeof param.args === "object" && typeof args === "object") {
    param.args = { ...args, ...param.args };
  } else {
    param.args = args;
  }
}
