
source("load-infomap.R")

im <- InfomapWrapper("-2 -d")

im$addLink(0, 1);
im$addLink(0, 2);
im$addLink(0, 3);
im$addLink(1, 0);
im$addLink(1, 2);
im$addLink(2, 1);
im$addLink(2, 0);
im$addLink(3, 0);
im$addLink(3, 4);
im$addLink(3, 5);
im$addLink(4, 3);
im$addLink(4, 5);
im$addLink(5, 4);
im$addLink(5, 3);

im$run();

cat("Partitioned network in", im$numTopModules(), "modules with codelength", im$codelength(), "bits:\n")

it <- im$iterLeafNodes()

while (!it$isEnd()) {
	cat("Node:", it$stateId, "module:", it$moduleId(), '\n')
	it$stepForward()
}
