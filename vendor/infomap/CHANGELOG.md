# Changelog

All notable changes to this project will be documented in this file. See [standard-version](https://github.com/conventional-changelog/standard-version) for commit guidelines.

### [2.6.1](https://github.com/mapequation/infomap/compare/v2.6.0...v2.6.1) (2022-10-31)


### Bug Fixes

* Scale variable Markov time logarithmic as default ([#320](https://github.com/mapequation/infomap/issues/320)) ([3f40134](https://github.com/mapequation/infomap/commit/3f40134b0346586b9d051c3a2996a77500ba10d6))
* Set num trials to 1 if no infomap ([be60a09](https://github.com/mapequation/infomap/commit/be60a09c785d8bee9ac885ebc92473494db368c1))

## [2.6.0](https://github.com/mapequation/infomap/compare/v2.5.0...v2.6.0) (2022-08-09)


### Features

* **python:** Add entropy_rate property ([a243dfc](https://github.com/mapequation/infomap/commit/a243dfcdf4629aff85b425d11aa7317b1fc6b80c))
* Variable Markov time ([#315](https://github.com/mapequation/infomap/issues/315)) ([cd7b25b](https://github.com/mapequation/infomap/commit/cd7b25b8aded24f85034ddc0d2d3bbacb881b24c))


### Bug Fixes

* Adjust entropy correction to Miller-Madow ([324c9f7](https://github.com/mapequation/infomap/commit/324c9f757d94b2d045ec77d14b53491ef30aba17))
* Don't write multilevel json output for aggregated higher-order networks ([1f5e7ea](https://github.com/mapequation/infomap/commit/1f5e7eab31b865fd2208cecbd0ed50642b9ca199))
* Fix calc entropy rate for undirected networks ([ccd6b6c](https://github.com/mapequation/infomap/commit/ccd6b6cbdb8b322a5e457355a214801301351d4a))
* **py:** Fix ignored variable_markov_time argument ([eedc777](https://github.com/mapequation/infomap/commit/eedc77705e44255c580d973c9403d1062a6e1712))
* **py:** Fix test for entropy rate ([78772b5](https://github.com/mapequation/infomap/commit/78772b5b81323d65456220e23456210c48fc0871))
* **python:** Fix physical leaf node iterator on state networks ([#312](https://github.com/mapequation/infomap/issues/312)) ([4e2918c](https://github.com/mapequation/infomap/commit/4e2918c131beeda4187012e8c9672e3635f77df4)), closes [#313](https://github.com/mapequation/infomap/issues/313)
* **python:** Only aggregate states into physical nodes when we have a higher-order network ([7e3277c](https://github.com/mapequation/infomap/commit/7e3277c97881057b36756f6328089ed7252234eb))
* **python:** Segfault when accessing physical_nodes ([cc5c95e](https://github.com/mapequation/infomap/commit/cc5c95e302e5093ee168a9952bbb70dc1e5b74ce)), closes [#300](https://github.com/mapequation/infomap/issues/300)

## [2.5.0](https://github.com/mapequation/infomap/compare/v2.4.1...v2.5.0) (2022-06-07)


### Features

* Output per-module codelength in json ([eae6650](https://github.com/mapequation/infomap/commit/eae6650c48dbf90a938d011f1c736794a786d7ef))
* Print all trials using --print-all-trials ([e65bea8](https://github.com/mapequation/infomap/commit/e65bea86c21ebcabafb053e79306049e92429ca1)), closes [#298](https://github.com/mapequation/infomap/issues/298)
* Write node metadata to json ([472b791](https://github.com/mapequation/infomap/commit/472b791742f0f692a0a4b06262b33165245f6f36))


### Bug Fixes

* Constrain numerical arguments to sensible values ([174c99f](https://github.com/mapequation/infomap/commit/174c99fb56c4a50a2e11c37208bcec3e98ba4ab4))
* Don't allow less than 1 trial ([18e97e2](https://github.com/mapequation/infomap/commit/18e97e21d76e7d27520eb2625c281031d72376d3))
* Json multilevel modules for higher-order networks ([d020c7c](https://github.com/mapequation/infomap/commit/d020c7c444d1227cb26aa249f7cae5e7de5ecaa4)), closes [#266](https://github.com/mapequation/infomap/issues/266)
* Only print all trials when running with more than 1 trial ([5f9d8b5](https://github.com/mapequation/infomap/commit/5f9d8b5e56ae348b142404f3f775ad526d396c55))
* **python:** Fix elapsed running time ([1131e84](https://github.com/mapequation/infomap/commit/1131e845243d3053cceae5bdf069569b2856bc4a))
* **python:** Replace MersenneTwister implementation with stdlib ([be3effd](https://github.com/mapequation/infomap/commit/be3effd9813a00892437c806b62cb8ee15387d3d))

### [2.4.1](https://github.com/mapequation/infomap/compare/v2.4.0...v2.4.1) (2022-05-27)


### Bug Fixes

* Json output causes segfault ([78e80c5](https://github.com/mapequation/infomap/commit/78e80c581edd92001818abdf50066b798251c604))

## [2.4.0](https://github.com/mapequation/infomap/compare/v2.3.0...v2.4.0) (2022-05-24)


### Features

* **js:** React hook ([c19a819](https://github.com/mapequation/infomap/commit/c19a81917849e27298d00bc1b7c6bde2d66f46c5))


### Bug Fixes

* **js-parser:** Parse module level from clu files ([bb91066](https://github.com/mapequation/infomap/commit/bb9106625761716dfabbb240ba39702c110f011a))
* **js:** Correctly parse clu-files without header ([04c39ea](https://github.com/mapequation/infomap/commit/04c39ea46b716fb10af987a568d5ab226e6ad69d))
* **python:** Don't allow adding multilayer intra-inter links and ordinary nodes. ([8c3d09c](https://github.com/mapequation/infomap/commit/8c3d09c5c5188c1c94d691beff2c519100b939f9)), closes [#287](https://github.com/mapequation/infomap/issues/287)
* **python:** Improve constructor and run signatures ([6852e14](https://github.com/mapequation/infomap/commit/6852e142343c91e61a9affe9bcff2897375cf70c)), closes [#286](https://github.com/mapequation/infomap/issues/286)
* **python:** no_self_links was always set unless include_self_links=True ([9afe5fa](https://github.com/mapequation/infomap/commit/9afe5fa7afbf3f143dd2b8d293be59f572a14c79))
* **python:** Prevent infinite loop in getMultilevelModules ([1174b66](https://github.com/mapequation/infomap/commit/1174b6658ed5c0bbc2c976223524f5f58e4d6c35)), closes [#284](https://github.com/mapequation/infomap/issues/284)
* Using deprecated option --include-self-links references missing option --no-loops ([1497bf7](https://github.com/mapequation/infomap/commit/1497bf71fe7a6f5dad725742db8c36a49171eb4d))

## [2.3.0](https://github.com/mapequation/infomap/compare/v2.2.0...v2.3.0) (2022-04-13)


### Features

* Matchable multilayer ids ([c9d28b6](https://github.com/mapequation/infomap/commit/c9d28b6e739731043e0dd6b82d0ae9898f686618)), closes [#279](https://github.com/mapequation/infomap/issues/279)


### Bug Fixes

* **js:** Don't call error on worker terminate ([cc57c8c](https://github.com/mapequation/infomap/commit/cc57c8cbafa3ec313426e1356a8140fa6d09c9b5))
* **js:** More robust heuristic to parse output format type ([33bde00](https://github.com/mapequation/infomap/commit/33bde00775c2fd41fb03b7592ae80633f09e8263))
* Match multilayer nodes by node and layer ids when using --cluster-data ([b5f6eb3](https://github.com/mapequation/infomap/commit/b5f6eb38d1740559ae93591e659ff8f2f53b7927))
* **python:** Add no_self_links to api ([b3d903d](https://github.com/mapequation/infomap/commit/b3d903d16e7619edebabe58a9004f40c88ea3614)), closes [#285](https://github.com/mapequation/infomap/issues/285)
* **python:** Add recorded teleportation keyword argument ([839e1e0](https://github.com/mapequation/infomap/commit/839e1e057b2f7181a87d77758ae00c20aa368891)), closes [#270](https://github.com/mapequation/infomap/issues/270)
* **python:** Fix matchable multilayer id option ([8377e07](https://github.com/mapequation/infomap/commit/8377e077832eafba1cb989dcad6e24b7a8f87edf))

## [2.2.0](https://github.com/mapequation/infomap/compare/v2.1.0...v2.2.0) (2022-03-07)


### Features

* Add flow-model, higher-order, and state-level fields to json, tree and clu output ([24125e5](https://github.com/mapequation/infomap/commit/24125e5874cd0599e7274e08a2691855f212186b))


### Bug Fixes

* **js:** More forgiving clu and tree parsing ([4764946](https://github.com/mapequation/infomap/commit/47649466a0ef2a200dc766b1442c7f33a8fec20f))
* Write json output with --no-infomap ([c0cec87](https://github.com/mapequation/infomap/commit/c0cec8715981f53dd53c3134c8f16e5b4d644717)), closes [#268](https://github.com/mapequation/infomap/issues/268)

## [2.1.0](https://github.com/mapequation/infomap/compare/v2.0.2...v2.1.0) (2022-02-11)


### Features

* **js:** Add progress event callback ([0fab95b](https://github.com/mapequation/infomap/commit/0fab95bd59717aaf2f36792f2be5a54223e38fab))


### Bug Fixes

* **js:** Change Typescript types to use mec instead of modularCentrality ([298eed3](https://github.com/mapequation/infomap/commit/298eed32d977758b07af058327775c234340b7b5))
* **js:** Parse node path as string instead of array of number ([51feb73](https://github.com/mapequation/infomap/commit/51feb73254b17857e4feb27197d4bb8c8157db05))
* **js:** Validate parsed results and throw instead of returning null on failure ([5f56b94](https://github.com/mapequation/infomap/commit/5f56b94fe0f956483f61a03ef07e94d8def2205e))

### [2.0.2](https://github.com/mapequation/infomap/compare/v2.0.1...v2.0.2) (2022-01-31)


### Bug Fixes

* Change json key modularCentrality to mec ([78c7a57](https://github.com/mapequation/infomap/commit/78c7a573914e127bda03a2b8001df321e5bffd04))
* **js:** Parse directed field in ftree files ([916f54a](https://github.com/mapequation/infomap/commit/916f54a05fb2ffd86ddf0ba297952e85ccfaa5b8))
* **js:** Pass parseLinks to main parse function and correctly detect multilayer tree files ([0a40ca6](https://github.com/mapequation/infomap/commit/0a40ca6c1606c189f02e691a2fdb74a4494ffbe7))
* Minify json output ([123c2b1](https://github.com/mapequation/infomap/commit/123c2b10f1e6486135346f8cd546ca0e72fb4d2e)), closes [#260](https://github.com/mapequation/infomap/issues/260)

### [2.0.1](https://github.com/mapequation/infomap/compare/v2.0.0...v2.0.1) (2022-01-26)


### Bug Fixes

* Bug in write state network to json ([50f8de2](https://github.com/mapequation/infomap/commit/50f8de2008cf6422693a74279b347d043735200d)), closes [#262](https://github.com/mapequation/infomap/issues/262)
* **js:** Broken parser due to duplicate variable name ([f437e13](https://github.com/mapequation/infomap/commit/f437e138f7b24140732831e63fac8c28e44e44bf))
* Use min mass for dangling nodes when regularized ([f741d07](https://github.com/mapequation/infomap/commit/f741d079b04558102130cbbbd742a0ddc3c1c366))

## [2.0.0](https://github.com/mapequation/infomap/compare/v1.9.0...v2.0.0) (2022-01-13)


### ⚠ BREAKING CHANGES

* Removes --input-format.
* Removes --include-self-links.

### Features

* Entropy bias correction ([#258](https://github.com/mapequation/infomap/issues/258)) ([8ea235c](https://github.com/mapequation/infomap/commit/8ea235cf9159fc1a672a0470798ab7e443590fde))
* Regularized map equation ([#181](https://github.com/mapequation/infomap/issues/181)) ([0673000](https://github.com/mapequation/infomap/commit/06730005329acb4ba781d647e7009d4ff2510ae7)), closes [#256](https://github.com/mapequation/infomap/issues/256)


### Bug Fixes

* Add new feature flags to python and js api ([48245a7](https://github.com/mapequation/infomap/commit/48245a795f87675b62b94b2f8bd8251b71e5f62e))
* Include self links by default ([#255](https://github.com/mapequation/infomap/issues/255)) ([1f68940](https://github.com/mapequation/infomap/commit/1f68940a90633405f696efa93471ca8cbb25f14c))
* **example:** Fix building c++ examples ([0d5c063](https://github.com/mapequation/infomap/commit/0d5c0639cb76b9b275388c1a45873496f675aeea)), closes [#252](https://github.com/mapequation/infomap/issues/252)


* Remove superfluous --input-format flag ([2dd77c4](https://github.com/mapequation/infomap/commit/2dd77c4689c67485454478d739c22c8deb44eda0))

## [1.9.0](https://github.com/mapequation/infomap/compare/v1.8.0...v1.9.0) (2021-11-17)


### Features

* **python:** Add API for intra/inter-layer links ([#245](https://github.com/mapequation/infomap/issues/245)) ([aa31ab2](https://github.com/mapequation/infomap/commit/aa31ab22ff8a706f7b1e5bbf75740d35ef5716bf)), closes [#244](https://github.com/mapequation/infomap/issues/244)


### Bug Fixes

* Incorrect one-level codelength for higher-order input ([5a6d1a1](https://github.com/mapequation/infomap/commit/5a6d1a10657c568f0fddf3afd6729a63e3543081)), closes [#182](https://github.com/mapequation/infomap/issues/182)
* Write "directed", not "1", in cli output notice ([f0a9dc7](https://github.com/mapequation/infomap/commit/f0a9dc722f8e8d5b4ca7e9e30d201f129d210df4))

## [1.8.0](https://github.com/mapequation/infomap/compare/v1.7.3...v1.8.0) (2021-11-05)


### Features

* **js:** Add clu and tree javascript parser ([f0c738c](https://github.com/mapequation/infomap/commit/f0c738c2ebbd4d1f43d14b75f82767b30dbfb02f))
* Add modular centrality score ([31cdd7f](https://github.com/mapequation/infomap/commit/31cdd7fa7a12add63f5b8818e503fe81b3f48555))

### [1.7.3](https://github.com/mapequation/infomap/compare/v1.7.2...v1.7.3) (2021-10-26)


### Bug Fixes

* **js:** Don't run error callback on finish ([13a21f0](https://github.com/mapequation/infomap/commit/13a21f00f445e8a378ca287e080085b47df4653f))

### [1.7.2](https://github.com/mapequation/infomap/compare/v1.7.1...v1.7.2) (2021-10-26)


### Bug Fixes

* **js:** Run error callback on terminate ([e83e897](https://github.com/mapequation/infomap/commit/e83e8970656e654aee353705a85bd7f97aae3f71))

### [1.7.1](https://github.com/mapequation/infomap/compare/v1.7.0...v1.7.1) (2021-10-14)

## [1.7.0](https://github.com/mapequation/infomap/compare/v1.6.0...v1.7.0) (2021-10-04)


### Features

* **python:** Add multilevel_modules property ([0d4cc37](https://github.com/mapequation/infomap/commit/0d4cc375cebe8477ea2e0073f9e7056411c3091e))
* **python:** Add num_levels property ([abd4f7f](https://github.com/mapequation/infomap/commit/abd4f7f8cc272cc20e4996b03a10dd2411a46260))
* **python:** Pass dict to add_nodes ([e30870f](https://github.com/mapequation/infomap/commit/e30870f28366676e105b189980a2c5eb921b9481))
* **python:** Pass dict to add_state_nodes ([5cd9ba2](https://github.com/mapequation/infomap/commit/5cd9ba25db758d2cb1a928c32adfdea21dda87d0))
* **python:** Pass dict to set_names ([0ed79e8](https://github.com/mapequation/infomap/commit/0ed79e8b08ab1ae410a3fdbb4dd59ddb2ecf149a))
* Add directed field to json output ([7db19f2](https://github.com/mapequation/infomap/commit/7db19f216b3d7f825a7eb90c89e0a0696c3a3e95))
* Write json multilevel modules and module information ([e6796c8](https://github.com/mapequation/infomap/commit/e6796c8d91499c1287278ab818c7e62d7426c303))


### Bug Fixes

* **python:** Adding a networkx DiGraph should set the directed flag ([a60a404](https://github.com/mapequation/infomap/commit/a60a404eac88d04a3bb6beb1b130452b28a99d83))
* **python:** Using directed=False should use undirected flow ([be0f6e2](https://github.com/mapequation/infomap/commit/be0f6e2f9ffc80ee58305709ebc6b801187e21bb))
* Metadata codelength works after multiple runs ([ec00ccd](https://github.com/mapequation/infomap/commit/ec00ccd7d14bdff54e6757b4cd1fd84e93dd4e58)), closes [#167](https://github.com/mapequation/infomap/issues/167)
* Module ids wrong for uneven depth and chosen level > 1 ([c3a0c90](https://github.com/mapequation/infomap/commit/c3a0c9025ae82a84f05e929f105b88d1d28f8be1)), closes [#221](https://github.com/mapequation/infomap/issues/221)

## [1.6.0](https://github.com/mapequation/infomap/compare/v1.5.1...v1.6.0) (2021-07-08)


### Features

* **python:** Add get_dataframe ([13066b3](https://github.com/mapequation/infomap/commit/13066b3300c921dcb18b46be2a633153ecf0ffa8))
* **python:** Add write method ([fdb3ce1](https://github.com/mapequation/infomap/commit/fdb3ce13a96c77b46b064fcda5fdcaf0f0a3c753))
* **python:** Add write_state_network and write_pajek ([da715a9](https://github.com/mapequation/infomap/commit/da715a9ccfc27b6dd2a079bd9f4e3ef6501a95a1))


### Bug Fixes

* Intra-module links missing from ftree output ([b347e24](https://github.com/mapequation/infomap/commit/b347e2491f343f2a541e6a3e4b7ac8837d0ae507))
* **js:** Add type to files parameter ([269b426](https://github.com/mapequation/infomap/commit/269b426326e7eca02b53dcbddf3ae0925a4b5a23))
* **js:** Output arguments should not be required to be an array ([a292318](https://github.com/mapequation/infomap/commit/a29231853e4d900e5f99dee962b1f44a20fc4e93))
* **js:** Worker should accept extra files as objects ([524a0e7](https://github.com/mapequation/infomap/commit/524a0e7e80c2bef897dab85cd60473a3e06408c9))
* **python:** get_multilevel_modules should return a dict ([f050fb5](https://github.com/mapequation/infomap/commit/f050fb568df4fb678e6c4689bef85cf18a6d93f2))

### [1.5.1](https://github.com/mapequation/infomap/compare/v1.5.0...v1.5.1) (2021-06-23)


### Bug Fixes

* **js:** Add missing flow network file output ([6b8ddb8](https://github.com/mapequation/infomap/commit/6b8ddb820bc48eaed27a4e66d2810151b65b8bb9))
* **js:** Handle json parse exceptions in worker ([a8a127b](https://github.com/mapequation/infomap/commit/a8a127b02e81b1bbd2a919edbc7952724b2e8682))
* **js:** Pass arguments as object ([13e71a0](https://github.com/mapequation/infomap/commit/13e71a0b1a7465e5f59bb3d05ef1bbdef7b5f65f))
* **js:** Worker should accept network as object ([d7184b8](https://github.com/mapequation/infomap/commit/d7184b840e7aae3b7d59b9b8bb2fa432d8d4c616))
* Avoid nan values for relative codelength savings ([e6a429e](https://github.com/mapequation/infomap/commit/e6a429e6bd463fc405161f6b4e981439979fc5c3))

## [1.5.0](https://github.com/mapequation/infomap/compare/v1.4.4...v1.5.0) (2021-06-22)


### Features

* **js:** Add Typescript types to JS worker ([f2fecf4](https://github.com/mapequation/infomap/commit/f2fecf4898dd0200709def0da11bea1883e70eaa))
* **python:** Add get_links, links and flow_links ([a4c2611](https://github.com/mapequation/infomap/commit/a4c26112821abc05185442b1eb27fd1bf83fdbc9))
* Write tree to CSV file ([5d0d36a](https://github.com/mapequation/infomap/commit/5d0d36a58413c0755f97d93dd00d8f6db7b2b785))

### [1.4.4](https://github.com/mapequation/infomap/compare/v1.4.3...v1.4.4) (2021-06-14)


### Bug Fixes

* **python:** Fix OpenMP detection ([27f9070](https://github.com/mapequation/infomap/commit/27f9070d4676aee859f14bbf0ba2e1de6a417636))

### [1.4.3](https://github.com/mapequation/infomap/compare/v1.4.2...v1.4.3) (2021-06-14)

### [1.4.2](https://github.com/mapequation/infomap/compare/v1.4.1...v1.4.2) (2021-06-14)


### Bug Fixes

* **build:** Fix building Python package with OpenMP if available ([cfaade5](https://github.com/mapequation/infomap/commit/cfaade5bb699980954fc6d618ef1eb211bb2a63e))

### [1.4.1](https://github.com/mapequation/infomap/compare/v1.4.0...v1.4.1) (2021-06-10)


### Bug Fixes

* **js:** Add newick and json to output from JS worker ([549ee84](https://github.com/mapequation/infomap/commit/549ee84858860680d2b9ed93e2d79a419c24e7d1))
* **python:** Add missing write_newick and write_json methods ([62c2736](https://github.com/mapequation/infomap/commit/62c2736dfb048de432093313e19b3d67869d3138))
* Change Newick output extention to .nwk ([f0e2e8e](https://github.com/mapequation/infomap/commit/f0e2e8ef8b6f5da6985c548b2f2b00521177cad8))

## [1.4.0](https://github.com/mapequation/infomap/compare/v1.3.0...v1.4.0) (2021-06-10)


### Features

* Support printing JSON tree ([22116d6](https://github.com/mapequation/infomap/commit/22116d6a9077909c2e5eb2160feea913022752e8))
* **python:** Accept named arguments ([5dabf4a](https://github.com/mapequation/infomap/commit/5dabf4abbaa90512d5a1ffd11c8986f8623ccb13))
* **python:** Add 'add_networkx_graph' method to API ([4067d40](https://github.com/mapequation/infomap/commit/4067d40268c4b2726652b5ef47358ca049c58edd))
* **python:** Add effective number of modules ([0f9a516](https://github.com/mapequation/infomap/commit/0f9a5169cacde5c9106bb6cf7dfa3b726a1c28d5))
* **python:** Add shorthand for accessing number of nodes and links ([8808045](https://github.com/mapequation/infomap/commit/88080454e054c51516ee7047975baf480ce1bd42))
* Add 'flow' to output options for node/edge flow ([d74a773](https://github.com/mapequation/infomap/commit/d74a773dd33102b85ff444158c958a8048c8cb02))


### Bug Fixes

* **python:** Also add nodes in networkx example ([3a88878](https://github.com/mapequation/infomap/commit/3a888789cc229521cffdf6da62fd0c239a5607eb))
* **python:** Fix effective_num_modules ([f41827f](https://github.com/mapequation/infomap/commit/f41827fde6f30e914e73d8a141992c4363f14adc))
* **R:** Expose base class methods to avoid c stack usage error ([47bdec0](https://github.com/mapequation/infomap/commit/47bdec0ab9a3d0e80cd7a77a7629ad9d9b400183))
* **R:** Fix swig build by hiding python specifics ([6b8a67e](https://github.com/mapequation/infomap/commit/6b8a67e6f09a02daad8e85c8b9fcec0f632526e3))
* **R:** Replace enum class with struct for R ([73c702f](https://github.com/mapequation/infomap/commit/73c702fd23e7a02ae8ccc3630f5713181cb842a9))
* **R:** Update R example and build to Infomap v1 ([17e1cc0](https://github.com/mapequation/infomap/commit/17e1cc00bd789c8798bd0633a39285a4160854a4))
* **R:** Update swig version in docker build ([74c4e86](https://github.com/mapequation/infomap/commit/74c4e86cfdd395ec1b335fc34dd5e0b34d7f6c1e))

## [1.3.0](https://github.com/mapequation/infomap/compare/v1.2.1...v1.3.0) (2020-12-03)


### Features

* Add --hide-bipartite-nodes ([#156](https://github.com/mapequation/infomap/issues/156)) ([f5e673c](https://github.com/mapequation/infomap/commit/f5e673c64efaf5ae825ee6b303e78389b8d594f9))
* New optional flow model for bipartite input ([91f35fe](https://github.com/mapequation/infomap/commit/91f35fef9b1246b87cde46d8e36ad7a992291a28))
* Support Newick tree format in output ([9749a05](https://github.com/mapequation/infomap/commit/9749a0535b9db47826c198b6f98b8fe2510eadc4))


### Bug Fixes

* Continue searching from input tree ([2199921](https://github.com/mapequation/infomap/commit/21999211ed822a8e46a8f1cd158185ce79b3075e))
* Fine-tune bottom modules in input tree with -F ([c187718](https://github.com/mapequation/infomap/commit/c187718e49618e2dac6c060233c844ff9a40893f))
* Implement missing multilayer relax by JSD ([95141bf](https://github.com/mapequation/infomap/commit/95141bf89452053a5f8f7bb2afbfc4048f066ee3)), closes [#154](https://github.com/mapequation/infomap/issues/154)
* Use correct objective in coarse tune for memory and multilevel networks ([97db558](https://github.com/mapequation/infomap/commit/97db558f27a81994a09aa6834bff5d8810a93130)), closes [#174](https://github.com/mapequation/infomap/issues/174)

### [1.2.1](https://github.com/mapequation/infomap/compare/v1.2.0...v1.2.1) (2020-11-10)


### Bug Fixes

* --clu-level was not respected in command line usage ([b82a0b1](https://github.com/mapequation/infomap/commit/b82a0b1a0582583b4552887650b0a6b92aba14fc)), closes [#153](https://github.com/mapequation/infomap/issues/153)
* Allow ftree cluster data ([bfe3601](https://github.com/mapequation/infomap/commit/bfe3601a50675cf353379efdb940b184773fb2b2))
* Apply weight-threshold also to multilayer nodes ([db2d781](https://github.com/mapequation/infomap/commit/db2d781c5257932ad5a4e5fe1976498b01985a39))
* Remove ClusterReader/BipartiteClusterReader ([8e08309](https://github.com/mapequation/infomap/commit/8e0830915f47bd4460a1b1b4f3794d60919a429a)), closes [#162](https://github.com/mapequation/infomap/issues/162)
* **windows:** Include algorithm header ([4f9ef88](https://github.com/mapequation/infomap/commit/4f9ef88409012ff00457291fc0ad8d0f33063c59)), closes [#150](https://github.com/mapequation/infomap/issues/150)

## [1.2.0](https://github.com/mapequation/infomap/compare/v1.1.4...v1.2.0) (2020-09-30)


### Features

* Add option to --use-node-weights-as-flow ([0d75174](https://github.com/mapequation/infomap/commit/0d751745ba565e06f93dfad6feb79ac9c6794e1e))
* Show bipartite start id in bipartite output ([f6e9233](https://github.com/mapequation/infomap/commit/f6e9233932f39f2dbe72b842fe2ddb259b725e08))


### Bug Fixes

* Fix flow calculator for bipartite networks ([b725b27](https://github.com/mapequation/infomap/commit/b725b27563e84092438839606aba294d7b993de5))
* **python:** Fix bipartite_start_id getter ([d326298](https://github.com/mapequation/infomap/commit/d326298c22570ded7002a8ee9cc361a448d8e4c5))

### [1.1.4](https://github.com/mapequation/infomap/compare/v1.1.3...v1.1.4) (2020-09-13)


### Bug Fixes

* Clarify that the algorithm searches for multi-level solutions by default ([47be4de](https://github.com/mapequation/infomap/commit/47be4defa349e3be99851105aff290cd488b9282))
* **python:** Fix color indices in networkx example ([dcf5888](https://github.com/mapequation/infomap/commit/dcf588800faf7f5d1f9be13d0bcc5528a2c2cf18))

### [1.1.3](https://github.com/mapequation/infomap/compare/v1.1.2...v1.1.3) (2020-05-09)


### Bug Fixes

* **python:** Correct codelength terms after N > 1 ([146cee6](https://github.com/mapequation/infomap/commit/146cee6246e55fb9c4a711c1278bdc251b54fced))
* **python:** Define module_codelength as L - L_index ([fa6eb94](https://github.com/mapequation/infomap/commit/fa6eb946acb22d808cd286bcfad6d76d3f2d4293))
* **python:** Don't create empty node names if they didn't already exist ([#126](https://github.com/mapequation/infomap/issues/126)) ([4d1cd56](https://github.com/mapequation/infomap/commit/4d1cd565e356f5c4ac4e8ce1620ae128a0f9ad4c))
* **python:** Expose num_non_trivial_top_modules ([6b91184](https://github.com/mapequation/infomap/commit/6b911842693f46546fdc89beb145394ed74dacdf))

### [1.1.2](https://github.com/mapequation/infomap/compare/v1.1.1...v1.1.2) (2020-04-08)


### Bug Fixes

* **python:** Expose node names with get_name and get_names ([#125](https://github.com/mapequation/infomap/issues/125)) ([b8cb409](https://github.com/mapequation/infomap/commit/b8cb409eec171d63f60634ce3a19cc01bc3a6eef)), closes [#52](https://github.com/mapequation/infomap/issues/52)

### [1.1.1](https://github.com/mapequation/infomap/compare/v1.1.0...v1.1.1) (2020-04-08)


### Bug Fixes

* Remove implied link direction from --skip-adjust-bipartite-flow ([05c3eba](https://github.com/mapequation/infomap/commit/05c3eba891521e6833292802c9c9671dee1303f7))

## [1.1.0](https://github.com/mapequation/infomap/compare/v1.0.11...v1.1.0) (2020-03-31)


### Features

* Show num levels and modules in file output ([0b70cd8](https://github.com/mapequation/infomap/commit/0b70cd82030693faaf24a318d45b7e929dc61c77))
* **python:** Expose codelengths for all trials ([eab7ae8](https://github.com/mapequation/infomap/commit/eab7ae861b2d93a9aea640994b53a80f4a52030c))


### Bug Fixes

* Don't use physical names in _states_as_physical.net ([5889fb1](https://github.com/mapequation/infomap/commit/5889fb17b13a38b1ce7a3b5ae714810bba1406c2))

### [1.0.11](https://github.com/mapequation/infomap/compare/v1.0.10...v1.0.11) (2020-03-31)


### Bug Fixes

* Handle zero weights in intra-layer links ([ddc5a3d](https://github.com/mapequation/infomap/commit/ddc5a3d827013cc006bdfa4b3518c32aef2ec722))

### [1.0.10](https://github.com/mapequation/infomap/compare/v1.0.9...v1.0.10) (2020-03-17)


### Bug Fixes

* Handle unassigned nodes after input tree ([d337535](https://github.com/mapequation/infomap/commit/d33753500939912035e30db4b4dfaf38beeb3eb1)), closes [#119](https://github.com/mapequation/infomap/issues/119)

### [1.0.9](https://github.com/mapequation/infomap/compare/v1.0.8...v1.0.9) (2020-03-14)


### Bug Fixes

* **windows:** compilation error on std::min ([b6c1a41](https://github.com/mapequation/infomap/commit/b6c1a419ad7f732520723a5df621eaed4bb5d2da))

### [1.0.8](https://github.com/mapequation/infomap/compare/v1.0.7...v1.0.8) (2020-03-12)


### Bug Fixes

* Fix reconstruct physically merged state nodes ([dc91e1a](https://github.com/mapequation/infomap/commit/dc91e1ad93fe150c58994c9dd2ab3aebefd84ef8)), closes [#118](https://github.com/mapequation/infomap/issues/118)
* Remove debug output ([#117](https://github.com/mapequation/infomap/issues/117)) ([60e7ab3](https://github.com/mapequation/infomap/commit/60e7ab346f50d3cedff06a6b7a12ef106b296f5a))

### [1.0.7](https://github.com/mapequation/infomap/compare/v1.0.6...v1.0.7) (2020-03-09)


### Bug Fixes

* Add header in network and states output ([c415a4b](https://github.com/mapequation/infomap/commit/c415a4bf05625eaa2330e6e9e3c9b1e2797895c5))
* Use *Edges/*Arcs instead of *Links in pajek output ([ff3a98e](https://github.com/mapequation/infomap/commit/ff3a98e18c119a64b8324018920452b6611b77d2))
* Use node_id instead of id in clu output ([39cd4e4](https://github.com/mapequation/infomap/commit/39cd4e4d866ee03c1d701946225a2bd76a0a9790))
* **js:** Fix support for all file io in Infomap.js ([ea9d899](https://github.com/mapequation/infomap/commit/ea9d899bdebdd361f0ff90e673471f92bc46cfaa))

### [1.0.6](https://github.com/mapequation/infomap/compare/v1.0.5...v1.0.6) (2020-03-03)


### Bug Fixes

* **python:** Enable numpy.int64 in link weights ([8d91206](https://github.com/mapequation/infomap/commit/8d9120681f932bd0cfad37cc4740f719f2acb358)), closes [#107](https://github.com/mapequation/infomap/issues/107)
* **python:** Print tree by default with python cli ([1565cab](https://github.com/mapequation/infomap/commit/1565cab1f3afec8c555180dd620379557709d3d0)), closes [#106](https://github.com/mapequation/infomap/issues/106)

### [1.0.5](https://github.com/mapequation/infomap/compare/v1.0.4...v1.0.5) (2020-03-02)


### Bug Fixes

* Common parameters should not be advanced ([#101](https://github.com/mapequation/infomap/issues/101)) ([2907c86](https://github.com/mapequation/infomap/commit/2907c8605313d60fc7c74e36d7a841b2f96c7b1d)), closes [#100](https://github.com/mapequation/infomap/issues/100)
* consistently use 1-based indexing for paths and 0 for indexes ([f283818](https://github.com/mapequation/infomap/commit/f2838189771cfc118ee20152d96694b438168ca8)), closes [#103](https://github.com/mapequation/infomap/issues/103)
* Fix ftree links since remapping path from 1 ([f447b48](https://github.com/mapequation/infomap/commit/f447b4872cb29b651747c08cd919e837353b4354)), closes [#102](https://github.com/mapequation/infomap/issues/102)

### [1.0.4](https://github.com/mapequation/infomap/compare/v1.0.3...v1.0.4) (2020-02-28)


### Bug Fixes

* **js:** Revert attempted worker blob optimization ([c521373](https://github.com/mapequation/infomap/commit/c521373fadaf741a34792ea27cd2efba6813cff9))

### [1.0.3](https://github.com/mapequation/infomap/compare/v1.0.2...v1.0.3) (2020-02-28)


### Bug Fixes

* **js:** Handle Infomap exceptions in js ([7a641f9](https://github.com/mapequation/infomap/commit/7a641f94ddae3df99c2065e77e3a1bc78922e174)), closes [#99](https://github.com/mapequation/infomap/issues/99)
* **python:** Start module id from 1 in get_[multilevel_]modules ([9419798](https://github.com/mapequation/infomap/commit/94197982597df86ae315dd9520ae24c27f67c552))

### [1.0.2](https://github.com/mapequation/infomap/compare/v1.0.1...v1.0.2) (2020-02-28)


### Bug Fixes

* **windows:** Fix compilation error on std::min ([3b000dc](https://github.com/mapequation/infomap/commit/3b000dcea45875e222d69b63ae0c289bcc6efcdc)), closes [#97](https://github.com/mapequation/infomap/issues/97)

### [1.0.1](https://github.com/mapequation/infomap/compare/v1.0.0...v1.0.1) (2020-02-27)


### Bug Fixes

* **python:** Fix missing package_meta on pip install ([dd24b1d](https://github.com/mapequation/infomap/commit/dd24b1d8fa7aa0c0b276dd611c57ef8a110002b8)), closes [#95](https://github.com/mapequation/infomap/issues/95)

## 1.0.0 (2020-02-26)


### ⚠ BREAKING CHANGES

* **python:** Drop support for python 2
* **python:** Rename physicalId to node_id.
Remove initial partition as first argument to run, use initial_partition argument.
Start InfomapIterator.path indexing from one.
* Output header format has changed
* Ftree output format has changed to include enter flow
* Clu output now contains top modules instead of leaf modules
* Tree output is sorted on flow
* Clu module ids starts from 1 instead of 0
* Output contains name of the physical node instead of state id
* Full multilayer format require the *multilayer heading
* Remove undirdir and outdirdir, use --flow-model. Rename min-improment to --core-loop-codelength-thresh. Remove random-loop-limit. Rename tune-iteration-threshold to --tune-iteration-relative-threshold. Rename skip-replace-to-one-module to --prefer-modular-solution.
* Removed --print-state-network, use -o states
* Renamed --set-unidentified-nodes-to-closest-module to
--assing-to-neighbouring-module
* No support for *path input

### Features

* **python:** Add meta info to infomap package ([a2caaf4](https://github.com/mapequation/infomap/commit/a2caaf4ff7f368c237c2752207b097028bacb900)), closes [#63](https://github.com/mapequation/infomap/issues/63)
* Add enter flow to ftree output ([bd3255e](https://github.com/mapequation/infomap/commit/bd3255e61977ce8f1d9e0babb95ec8158710e3a8)), closes [#82](https://github.com/mapequation/infomap/issues/82)
* Show entropy rate in the beginning ([f61692c](https://github.com/mapequation/infomap/commit/f61692c60c681d532f7a45b980f5041199b088b5))
* **python:** Improve Python API ([#87](https://github.com/mapequation/infomap/issues/87)) ([c39cd9f](https://github.com/mapequation/infomap/commit/c39cd9f2a310a10a30e0de83dab95902e8c6aa91))
* Add -o/--output option for comma-separated formats ([a255a18](https://github.com/mapequation/infomap/commit/a255a181f109b51eed3e1cbb82bb5f73f1894dd0))
* Add meta data in output file header ([66ccc64](https://github.com/mapequation/infomap/commit/66ccc64fec74a47a92b053642d7ad5fb63b74df2))
* Add option to print parameters in json ([3a2509d](https://github.com/mapequation/infomap/commit/3a2509d28dec022d5bdd2c8b8a4a5fd87448edab))
* Allow multilayer formats without specifying -i. ([08e09f1](https://github.com/mapequation/infomap/commit/08e09f136321dcf917d0d131d57041d32c923e64))
* Allow one-sided multilayer relax limit ([a116b83](https://github.com/mapequation/infomap/commit/a116b83bcf1240a17aca12a576fe634fd3db91dd)), closes [#68](https://github.com/mapequation/infomap/issues/68)
* Begin module id from 1 instead of 0 in clu output ([2bbea5d](https://github.com/mapequation/infomap/commit/2bbea5d537ef0e68ac167d1f84acbf0888d22fb8))
* Create npm package with Infomap worker, changelog and parameters ([44dae36](https://github.com/mapequation/infomap/commit/44dae36d65e9dccd2462db330298887318bc2463))
* Output name of physical node instead of state id ([dda3277](https://github.com/mapequation/infomap/commit/dda3277a2f97152a5b5011c94e440f9bb5ee9b9b))
* Remove --skip-complete-dangling-memory-nodes ([6909f21](https://github.com/mapequation/infomap/commit/6909f21a3ebc6eb84c5bc10f6f43c1839020ba41))
* Remove support for *path input ([e7be175](https://github.com/mapequation/infomap/commit/e7be175378838b294d64390bd6f59cdc6b8a42f1))
* Rename --set-unidentified-nodes-to-closest-module ([16fdc15](https://github.com/mapequation/infomap/commit/16fdc155d2d1fefbb22af1305f7ae7efcd1f06f4))
* Simplify command line interface ([f537132](https://github.com/mapequation/infomap/commit/f5371328e63d761ccaf6ee1cdce543ced80fe25c))
* Sort tree on flow ([6e27858](https://github.com/mapequation/infomap/commit/6e278584c45367145f57ce6323f6a8a9c2d64ffe)), closes [#71](https://github.com/mapequation/infomap/issues/71)
* Write top modules to .clu, not bottom level ([1e16a3c](https://github.com/mapequation/infomap/commit/1e16a3ccdf15ee6155a22d8d0cd4bc6ebfdeb94f))


### Bug Fixes

* **js:** Add missing this in error handler ([7345cdf](https://github.com/mapequation/infomap/commit/7345cdff8aa2b362df9d6522960fb992f5634f6e))
