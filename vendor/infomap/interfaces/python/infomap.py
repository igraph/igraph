import os
import warnings
from collections import namedtuple
from contextlib import contextmanager
from math import log2

if __name__ != "__main__":
    try:
        import pandas
    except ImportError:
        pandas = None


MultilayerNode = namedtuple("MultilayerNode", "layer_id, node_id")


def plogp(p):
    return (x * log2(x) if x > 0 else 0 for x in p)


def entropy(p):
    return -sum(plogp(p))


def perplexity(p):
    return 2 ** entropy(p)


_DEFAULT_META_DATA_RATE = 1.0
_DEFAULT_VERBOSITY_LEVEL = 1
_DEFAULT_TELEPORTATION_PROB = 0.15
_DEFAULT_MULTILAYER_RELAX_RATE = 0.15
_DEFAULT_SEED = 123
_DEFAULT_CORE_LOOP_CODELENGTH_THRESHOLD = 1e-10
_DEFAULT_TUNE_ITER_RELATIVE_THRESHOLD = 1e-5


def _construct_args(
    args=None,
    # input
    cluster_data=None,
    no_infomap=False,
    skip_adjust_bipartite_flow=False,
    bipartite_teleportation=False,
    weight_threshold=None,
    include_self_links=None,
    no_self_links=False,
    node_limit=None,
    matchable_multilayer_ids=None,
    assign_to_neighbouring_module=False,
    meta_data=None,
    meta_data_rate=_DEFAULT_META_DATA_RATE,
    meta_data_unweighted=False,
    # output
    tree=False,
    ftree=False,
    clu=False,
    verbosity_level=_DEFAULT_VERBOSITY_LEVEL,
    silent=False,
    out_name=None,
    no_file_output=False,
    clu_level=None,
    output=None,
    hide_bipartite_nodes=False,
    print_all_trials=False,
    # algorithm
    two_level=False,
    flow_model=None,
    directed=None,
    recorded_teleportation=False,
    use_node_weights_as_flow=False,
    to_nodes=False,
    teleportation_probability=_DEFAULT_TELEPORTATION_PROB,
    regularized=False,
    regularization_strength=1.0,
    entropy_corrected=False,
    entropy_correction_strength=1.0,
    markov_time=1.0,
    variable_markov_time=False,
    variable_markov_damping=1.0,
    preferred_number_of_modules=None,
    multilayer_relax_rate=_DEFAULT_MULTILAYER_RELAX_RATE,
    multilayer_relax_limit=-1,
    multilayer_relax_limit_up=-1,
    multilayer_relax_limit_down=-1,
    multilayer_relax_by_jsd=False,
    # accuracy
    seed=_DEFAULT_SEED,
    num_trials=1,
    core_loop_limit=10,
    core_level_limit=None,
    tune_iteration_limit=None,
    core_loop_codelength_threshold=_DEFAULT_CORE_LOOP_CODELENGTH_THRESHOLD,
    tune_iteration_relative_threshold=_DEFAULT_TUNE_ITER_RELATIVE_THRESHOLD,
    fast_hierarchical_solution=None,
    prefer_modular_solution=False,
    inner_parallelization=False,
):
    if args is None:
        args = ""

    # input
    if cluster_data is not None:
        args += " --cluster-data {}".format(cluster_data)

    if no_infomap:
        args += " --no-infomap"

    if skip_adjust_bipartite_flow:
        args += " --skip-adjust-bipartite-flow"

    if bipartite_teleportation:
        args += " --bipartite-teleportation"

    if weight_threshold is not None:
        args += " --weight-threshold {}".format(weight_threshold)

    if include_self_links is not None:
        warnings.warn(
            "include_self_links is deprecated, use no_self_links to exclude self-links",
            DeprecationWarning,
        )

    if include_self_links is not None and not include_self_links:
        args += " --no-self-links"

    if no_self_links:
        args += " --no-self-links"

    if node_limit is not None:
        args += " --node-limit {}".format(node_limit)

    if matchable_multilayer_ids is not None:
        args += " --matchable-multilayer-ids {}".format(matchable_multilayer_ids)

    if assign_to_neighbouring_module:
        args += " --assign-to-neightbouring-module"

    if meta_data is not None:
        args += " --meta-data {}".format(meta_data)

    if meta_data_rate != _DEFAULT_META_DATA_RATE:
        args += " --meta-data-rate {}".format(meta_data_rate)

    if meta_data_unweighted:
        args += " --meta-data-unweighted"

    # output
    if tree:
        args += " --tree"

    if ftree:
        args += " --ftree"

    if clu:
        args += " --clu"

    if verbosity_level > 1:
        args += " -{}".format("v" * verbosity_level)

    if silent:
        args += " --silent"

    if out_name is not None:
        args += " --out-name {}".format(out_name)

    if no_file_output:
        args += " --no-file-output"

    if clu_level is not None:
        args += " --clu-level {}".format(clu_level)

    if output is not None:
        args += " --output {}".format(",".join(output))

    if hide_bipartite_nodes:
        args += " --hide-bipartite-nodes"

    if print_all_trials:
        args += " --print-all-trials"

    # algorithm
    if two_level:
        args += " --two-level"

    if flow_model is not None:
        args += " --flow-model {}".format(flow_model)

    if directed is not None:
        args += " --directed" if directed else " --flow-model undirected"

    if recorded_teleportation:
        args += " --recorded-teleportation"

    if use_node_weights_as_flow:
        args += " --use-node-weights-as-flow"

    if to_nodes:
        args += " --to-nodes"

    if teleportation_probability != _DEFAULT_TELEPORTATION_PROB:
        args += " --teleportation-probability {}".format(teleportation_probability)

    if regularized:
        args += " --regularized"
    if regularization_strength != 1.0:
        args += " --regularization-strength {}".format(regularization_strength)
    if entropy_corrected:
        args += " --entropy-corrected"
    if entropy_correction_strength != 1.0:
        args += " --entropy-correction-strength {}".format(entropy_correction_strength)

    if markov_time != 1.0:
        args += " --markov-time {}".format(markov_time)
    if variable_markov_time:
        args += " --variable-markov-time"
    if variable_markov_damping != 1.0:
        args += " --variable-markov-damping {}".format(
            variable_markov_damping
        )

    if preferred_number_of_modules is not None:
        args += " --preferred-number-of-modules {}".format(preferred_number_of_modules)

    if multilayer_relax_rate != _DEFAULT_MULTILAYER_RELAX_RATE:
        args += " --multilayer-relax-rate {}".format(multilayer_relax_rate)

    if multilayer_relax_limit != -1:
        args += " --multilayer-relax-limit {}".format(multilayer_relax_limit)

    if multilayer_relax_limit_up != -1:
        args += " --multilayer-relax-limit-up {}".format(multilayer_relax_limit_up)

    if multilayer_relax_limit_down != -1:
        args += " --multilayer-relax-limit-down {}".format(multilayer_relax_limit_down)

    if multilayer_relax_by_jsd:
        args += " --multilayer-relax-by-jsd"

    # accuracy
    if seed != _DEFAULT_SEED:
        args += " --seed {}".format(seed)

    if num_trials != 1:
        args += " --num-trials {}".format(num_trials)

    if core_loop_limit != 10:
        args += " --core-loop-limit {}".format(core_loop_limit)

    if core_level_limit is not None:
        args += " --core-level-limit {}".format(core_level_limit)

    if tune_iteration_limit is not None:
        args += " --tune-iteration-limit {}".format(tune_iteration_limit)

    if core_loop_codelength_threshold != _DEFAULT_CORE_LOOP_CODELENGTH_THRESHOLD:
        args += " --core-loop-codelength-threshold {}".format(
            core_loop_codelength_threshold
        )

    if tune_iteration_relative_threshold != _DEFAULT_TUNE_ITER_RELATIVE_THRESHOLD:
        args += " --tune-iteration-relative-threshold {}".format(
            tune_iteration_relative_threshold
        )

    if fast_hierarchical_solution is not None:
        args += " -{}".format("F" * fast_hierarchical_solution)

    if prefer_modular_solution:
        args += " --prefer-modular-solution"

    if inner_parallelization:
        args += " --inner-parallelization"

    return args


class Infomap(InfomapWrapper):
    """Infomap

    This class is a thin wrapper around Infomap C++ compiled with SWIG.

    Examples
    --------
    Create an instance and add nodes and links:

    >>> from infomap import Infomap
    >>> im = Infomap(silent=True)
    >>> im.add_node(1)
    >>> im.add_node(2)
    >>> im.add_link(1, 2)
    >>> im.run()
    >>> im.codelength
    1.0


    Create an instance and read a network file:

    >>> from infomap import Infomap
    >>> im = Infomap(silent=True, num_trials=10)
    >>> im.read_file("ninetriangles.net")
    >>> im.run()
    >>> tol = 1e-4
    >>> abs(im.codelength - 3.4622731375264144) < tol
    True


    For more examples, see the examples directory.
    """

    def __init__(
        self,
        args=None,
        # input
        cluster_data=None,
        no_infomap=False,
        skip_adjust_bipartite_flow=False,
        bipartite_teleportation=False,
        weight_threshold=None,
        include_self_links=None,
        no_self_links=False,
        node_limit=None,
        matchable_multilayer_ids=None,
        assign_to_neighbouring_module=False,
        meta_data=None,
        meta_data_rate=_DEFAULT_META_DATA_RATE,
        meta_data_unweighted=False,
        # output
        tree=False,
        ftree=False,
        clu=False,
        verbosity_level=_DEFAULT_VERBOSITY_LEVEL,
        silent=False,
        out_name=None,
        no_file_output=False,
        clu_level=None,
        output=None,
        hide_bipartite_nodes=False,
        print_all_trials=False,
        # algorithm
        two_level=False,
        flow_model=None,
        directed=None,
        recorded_teleportation=False,
        use_node_weights_as_flow=False,
        to_nodes=False,
        teleportation_probability=_DEFAULT_TELEPORTATION_PROB,
        regularized=False,
        regularization_strength=1.0,
        entropy_corrected=False,
        entropy_correction_strength=1.0,
        markov_time=1.0,
        variable_markov_time=False,
        variable_markov_damping=1.0,
        preferred_number_of_modules=None,
        multilayer_relax_rate=_DEFAULT_MULTILAYER_RELAX_RATE,
        multilayer_relax_limit=-1,
        multilayer_relax_limit_up=-1,
        multilayer_relax_limit_down=-1,
        multilayer_relax_by_jsd=False,
        # accuracy
        seed=_DEFAULT_SEED,
        num_trials=1,
        core_loop_limit=10,
        core_level_limit=None,
        tune_iteration_limit=None,
        core_loop_codelength_threshold=_DEFAULT_CORE_LOOP_CODELENGTH_THRESHOLD,
        tune_iteration_relative_threshold=_DEFAULT_TUNE_ITER_RELATIVE_THRESHOLD,
        fast_hierarchical_solution=None,
        prefer_modular_solution=False,
        inner_parallelization=False,
    ):
        """Create a new Infomap instance.

        Parameters
        ----------
        args : str, optional
            String of Infomap arguments.
        cluster_data : str, optional
            Provide an initial two-level (clu format) or multi-layer (tree
            format) solution.
        no_infomap : bool, optional
            Don't run the optimizer. Useful to calculate codelength of provided
            cluster data or to print non-modular statistics.
        skip_adjust_bipartite_flow : bool, optional
            Skip distributing all flow from the bipartite nodes to the primary
            nodes.
        bipartite_teleportation : bool, optional
            Teleport like the bipartite flow instead of two-step (unipartite)
            teleportation.
        weight_threshold : float, optional
            Limit the number of links to read from the network. Ignore links
            with less weight than the threshold.
        include_self_links : bool, optional
            Include links with the same source and target node.
            DEPRECATED. Include self links by default now, exclude with no_self_links.
        no_self_links : bool, optional
            Exclude self links in the input network.
        node_limit : int, optional
            Limit the number of nodes to read from the network. Ignore links
            connected to ignored nodes.
        matchable_multilayer_ids : int, optional
            Construct state ids from node and layer ids that are consistent
            across networks for the same max number of layers.
            Set to at least the largest layer id among networks to match.
        assign_to_neighbouring_module : bool, optional
            Assign nodes without module assignments (from ``cluster_data``) to
            the module assignment of a neighbouring node if possible.
        meta_data : str, optional
            Provide meta data (clu format) that should be encoded.
        meta_data_rate : float, optional
            Metadata encoding rate. Default is to encode each step.
        meta_data_unweighted : bool, optional
            Don't weight meta data by node flow.
        tree : bool, optional
            Write a tree file with the modular hierarchy. Automatically enabled
            if no other output is specified.
        ftree : bool, optional
            Write a ftree file with the modular hierarchy including aggregated
            links between (nested) modules.
        clu : bool, optional
            Write a clu file with the top cluster ids for each node.
        verbosity_level : int, optional
            Verbose output on the console.
        silent : bool, optional
            No output on the console.
        out_name : str, optional
            Name for the output files, e.g.
            ``"[output_directory]/[out-name].tree"``
        no_file_output : bool, optional
            Don't write output to file.
        clu_level : int, optional
            For clu output, print modules at specified depth from root.
            Use ``-1`` for bottom level modules.
        output : list(str), optional
            List of output formats.
            Options: clu, tree, ftree, newick, json, csv, network, states.
        hide_bipartite_nodes : bool, optional
            Project bipartite solution to unipartite.
        print_all_trials : bool, optional
            Print all trials to separate files.
        two_level : bool, optional
            Optimize a two-level partition of the network.
            Default (``false``) is multi-level.
        flow_model : str, optional
            Specify flow model.
            Options: undirected, directed, undirdir, outdirdir, rawdir.
        directed : bool, optional
            Assume directed links. Shorthand for ``flow_model="directed"``.
        recorded_teleportation : bool, optional
            If teleportation is used to calculate the flow, also record it
            when minimizing codelength. Default False.
        use_node_weights_as_flow : bool, optional
            Use node weights (from api or after names in Pajek format) as flow,
            normalized to sum to 1.
        to_nodes : bool, optional
            Teleport to nodes instead of to links, assuming uniform node
            weights if no such input data.
        teleportation_probability : float, optional
            Probability of teleporting to a random node or link.
        regularized : bool, optional
            Effectively add a fully connected Bayesian prior network to not overfit
            due to missing links. Implies recorded teleportation.
        regularization_strength : float, optional
            Adjust relative strength of Bayesian prior network with this multiplier.
        entropy_corrected : bool, optional
            Correct for negative entropy bias in small samples (many modules).
        entropy_correction_strength : float, optional
            Increase or decrease the default entropy correction with this multiplier.
        markov_time : float, optional
            Scales link flow to change the cost of moving between modules.
            Higher values results in fewer modules.
        variable_markov_time : bool, optional
            Increase Markov time locally to level out link flow. Reduces risk of
            overpartitioning sparse areas while keeping high resolution in dense areas.
        variable_markov_damping : float, optional
            Damping parameter for variable Markov time, to scale with local effective 
            degree (0) or local entropy (1).
        preferred_number_of_modules : int, optional
            Penalize solutions the more they differ from this number.
        multilayer_relax_rate : float, optional
            Probability to relax the constraint to move only in the current
            layer.
        multilayer_relax_limit : int, optional
            Number of neighboring layers in each direction to relax to.
            If negative, relax to any layer.
        multilayer_relax_limit_up : int, optional
            Number of neighboring layers with higher id to relax to.
            If negative, relax to any layer.
        multilayer_relax_limit_down : int, optional
            Number of neighboring layers with lower id to relax to.
            If negative, relax to any layer.
        multilayer_relax_by_jsd : bool, optional
            Relax proportional to the out-link similarity measured by the
            Jensen-Shannon divergence.
        seed : int, optional
            A seed (integer) to the random number generator for reproducible
            results.
        num_trials : int, optional
            Number of outer-most loops to run before picking the best solution.
        core_loop_limit : int, optional
            Limit the number of loops that tries to move each node into the
            best possible module.
        core_level_limit : int, optional
            Limit the number of times the core loops are reapplied on existing
            modular network to search bigger structures.
        tune_iteration_limit : int, optional
            Limit the number of main iterations in the two-level partition
            algorithm. 0 means no limit.
        core_loop_codelength_threshold : float, optional
            Minimum codelength threshold for accepting a new solution in core
            loop.
        tune_iteration_relative_threshold : float, optional
            Set codelength improvement threshold of each new tune iteration to
            ``f`` times the initial two-level codelength.
        fast_hierarchical_solution : int, optional
            Find top modules fast. Use ``2`` to keep all fast levels.
            Use ``3`` to skip recursive part.
        prefer_modular_solution : bool, optional
            Prefer modular solutions even if they are worse than putting all
            nodes in one module.
        inner_parallelization : bool, optional
            Parallelize the inner-most loop for greater speed.
            This may give some accuracy tradeoff.
        """
        super().__init__(
            _construct_args(
                args,
                cluster_data=cluster_data,
                no_infomap=no_infomap,
                skip_adjust_bipartite_flow=skip_adjust_bipartite_flow,
                bipartite_teleportation=bipartite_teleportation,
                weight_threshold=weight_threshold,
                include_self_links=include_self_links,
                no_self_links=no_self_links,
                node_limit=node_limit,
                matchable_multilayer_ids=matchable_multilayer_ids,
                assign_to_neighbouring_module=assign_to_neighbouring_module,
                meta_data=meta_data,
                meta_data_rate=meta_data_rate,
                meta_data_unweighted=meta_data_unweighted,
                tree=tree,
                ftree=ftree,
                clu=clu,
                verbosity_level=verbosity_level,
                silent=silent,
                out_name=out_name,
                no_file_output=no_file_output,
                clu_level=clu_level,
                output=output,
                hide_bipartite_nodes=hide_bipartite_nodes,
                print_all_trials=print_all_trials,
                two_level=two_level,
                flow_model=flow_model,
                directed=directed,
                recorded_teleportation=recorded_teleportation,
                use_node_weights_as_flow=use_node_weights_as_flow,
                to_nodes=to_nodes,
                teleportation_probability=teleportation_probability,
                regularized=regularized,
                regularization_strength=regularization_strength,
                entropy_corrected=entropy_corrected,
                entropy_correction_strength=entropy_correction_strength,
                markov_time=markov_time,
                variable_markov_time=variable_markov_time,
                variable_markov_damping=variable_markov_damping,
                preferred_number_of_modules=preferred_number_of_modules,
                multilayer_relax_rate=multilayer_relax_rate,
                multilayer_relax_limit=multilayer_relax_limit,
                multilayer_relax_limit_up=multilayer_relax_limit_up,
                multilayer_relax_limit_down=multilayer_relax_limit_down,
                multilayer_relax_by_jsd=multilayer_relax_by_jsd,
                seed=seed,
                num_trials=num_trials,
                core_loop_limit=core_loop_limit,
                core_level_limit=core_level_limit,
                tune_iteration_limit=tune_iteration_limit,
                core_loop_codelength_threshold=core_loop_codelength_threshold,
                tune_iteration_relative_threshold=tune_iteration_relative_threshold,
                fast_hierarchical_solution=fast_hierarchical_solution,
                prefer_modular_solution=prefer_modular_solution,
                inner_parallelization=inner_parallelization,
            )
        )

    # ----------------------------------------
    # Input
    # ----------------------------------------

    def read_file(self, filename, accumulate=True):
        """Read network data from file.

        Parameters
        ----------
        filename : str
        accumulate : bool, optional
            If the network data should be accumulated to already added
            nodes and links. Default ``True``.
        """
        return super().readInputData(filename, accumulate)

    def add_node(self, node_id, name=None, teleportation_weight=None):
        """Add a node.

        See Also
        --------
        set_name
        add_nodes

        Parameters
        ----------
        node_id : int
        name : str, optional
        teleportation_weight : float, optional
            Used for teleporting between layers in multilayer networks.
        """
        if name is None:
            name = ""

        if len(name) and teleportation_weight is not None:
            return super().addNode(node_id, name, teleportation_weight)
        elif len(name) and teleportation_weight is None:
            return super().addNode(node_id, name)
        elif not len(name) and teleportation_weight is not None:
            return super().addNode(node_id, teleportation_weight)

        return super().addNode(node_id)

    def add_nodes(self, nodes):
        """Add nodes.

        See Also
        --------
        add_node

        Examples
        --------
        Add nodes

        >>> from infomap import Infomap
        >>> im = Infomap()
        >>> im.add_nodes(range(4))


        Add named nodes

        >>> from infomap import Infomap
        >>> im = Infomap()
        >>> nodes = (
        ...     (1, "Node 1"),
        ...     (2, "Node 2"),
        ...     (3, "Node 3")
        ... )
        >>> im.add_nodes(nodes)
        >>> im.names
        {1: 'Node 1', 2: 'Node 2', 3: 'Node 3'}


        Add named nodes with teleportation weights

        >>> from infomap import Infomap
        >>> im = Infomap()
        >>> nodes = (
        ...     (1, "Node 1", 0.5),
        ...     (2, "Node 2", 0.2),
        ...     (3, "Node 3", 0.8)
        ... )
        >>> im.add_nodes(nodes)
        >>> im.names
        {1: 'Node 1', 2: 'Node 2', 3: 'Node 3'}

        Add named nodes using dict

        >>> from infomap import Infomap
        >>> im = Infomap()
        >>> nodes = {
        ...     1: "Node 1",
        ...     2: "Node 2",
        ...     3: "Node 3"
        ... }
        >>> im.add_nodes(nodes)
        >>> im.names
        {1: 'Node 1', 2: 'Node 2', 3: 'Node 3'}


        Add named nodes with teleporation weights using dict

        >>> from infomap import Infomap
        >>> im = Infomap()
        >>> nodes = {
        ...     1: ("Node 1", 0.5),
        ...     2: ("Node 2", 0.2),
        ...     3: ("Node 3", 0.8)
        ... }
        >>> im.add_nodes(nodes)
        >>> im.names
        {1: 'Node 1', 2: 'Node 2', 3: 'Node 3'}


        Parameters
        ----------
        nodes : iterable of tuples or iterable of int or dict of int: str or dict of int: tuple of (str, float)
            Iterable of tuples on the form
            ``(node_id, [name], [teleportation_weight])``
        """
        try:
            for node, attr in nodes.items():
                if isinstance(attr, str):
                    self.add_node(node, attr)
                else:
                    self.add_node(node, *attr)
        except AttributeError:
            for node in nodes:
                if isinstance(node, int):
                    self.add_node(node)
                else:
                    self.add_node(*node)

    def add_state_node(self, state_id, node_id):
        """Add a state node.

        Notes
        -----
        If a physical node with id node_id does not exist, it will be created.
        If you want to name the physical node, use ``set_name``.

        Parameters
        ----------
        state_id : int
        node_id : int
            Id of the physical node the state node should be added to.
        """
        return super().addStateNode(state_id, node_id)

    def add_state_nodes(self, state_nodes):
        """Add state nodes.

        See Also
        --------
        add_state_node

        Examples
        --------
        With tuples

        >>> from infomap import Infomap
        >>> im = Infomap()
        >>> states = (
        ...     (1, 1),
        ...     (2, 1),
        ...     (3, 2)
        ... )
        >>> im.add_state_nodes(states)


        With dict

        >>> from infomap import Infomap
        >>> im = Infomap()
        >>> states = {
        ...     1: 1,
        ...     2: 1,
        ...     3: 2
        ... }
        >>> im.add_state_nodes(states)


        Parameters
        ----------
        state_nodes : iterable of tuples or dict of int: int
            Iterable of tuples of the form ``(state_id, node_id)``
            or dict of the form ``{state_id: node_id}``.
        """
        try:
            for node in state_nodes.items():
                self.add_state_node(*node)
        except AttributeError:
            for node in state_nodes:
                self.add_state_node(*node)

    def set_name(self, node_id, name):
        """Set the name of a node.

        Parameters
        ----------
        node_id : int
        name : str
        """
        if name is None:
            name = ""
        return super().addName(node_id, name)

    def set_names(self, names):
        """Set names to several nodes at once.

        Examples
        --------
        With tuples

        >>> from infomap import Infomap
        >>> im = Infomap()
        >>> names = (
        ...     (1, "Node 1"),
        ...     (2, "Node 2")
        ... )
        >>> im.set_names(names)
        >>> im.names
        {1: 'Node 1', 2: 'Node 2'}


        With dict

        >>> from infomap import Infomap
        >>> im = Infomap()
        >>> names = {
        ...     1: "Node 1",
        ...     2: "Node 2"
        ... }
        >>> im.set_names(names)
        >>> im.names
        {1: 'Node 1', 2: 'Node 2'}


        See Also
        --------
        set_name, names

        Parameters
        ----------
        names : iterable of tuples or dict of int: str
            Iterable of tuples on the form ``(node_id, name)``
            or dict of the form ``{node_id: name}``.
        """
        try:
            for name in names.items():
                self.set_name(*name)
        except AttributeError:
            for name in names:
                self.set_name(*name)

    def add_link(self, source_id, target_id, weight=1.0):
        """Add a link.

        Notes
        -----
        If the source or target nodes does not exist, they will be created.

        See Also
        --------
        remove_link

        Parameters
        ----------
        source_id : int
        target_id : int
        weight : float, optional
        """
        return super().addLink(source_id, target_id, weight)

    def add_links(self, links):
        """Add several links.

        Examples
        --------

        >>> from infomap import Infomap
        >>> im = Infomap()
        >>> links = (
        ...     (1, 2),
        ...     (1, 3)
        ... )
        >>> im.add_links(links)


        See Also
        --------
        add_link
        remove_link

        Parameters
        ----------
        links : iterable of tuples
            Iterable of tuples of int of the form
            ``(source_id, target_id, [weight])``
        """
        for link in links:
            self.add_link(*link)

    def remove_link(self, source_id, target_id):
        """Remove a link.

        Notes
        -----
        Removing links will not remove nodes if they become disconnected.

        See Also
        --------
        add_link

        Parameters
        ----------
        source_id : int
        target_id : int
        """
        return self.network.removeLink(source_id, target_id)

    def remove_links(self, links):
        """Remove several links.

        Examples
        --------

        >>> from infomap import Infomap
        >>> im = Infomap()
        >>> links = (
        ...     (1, 2),
        ...     (1, 3)
        ... )
        >>> im.add_links(links)
        >>> im.remove_links(links)
        >>> im.num_links
        0


        See Also
        --------
        remove_link

        Parameters
        ----------
        links : iterable of tuples
            Iterable of tuples of the form ``(source_id, target_id)``
        """
        for link in links:
            self.remove_link(*link)

    def add_multilayer_link(
        self, source_multilayer_node, target_multilayer_node, weight=1.0
    ):
        """Add a multilayer link.

        Adds a link between layers in a multilayer network.

        Examples
        --------
        Usage with tuples:

        >>> from infomap import Infomap
        >>> im = Infomap()
        >>> source_multilayer_node = (0, 1) # layer_id, node_id
        >>> target_multilayer_node = (1, 2) # layer_id, node_id
        >>> im.add_multilayer_link(source_multilayer_node, target_multilayer_node)


        Usage with MultilayerNode

        >>> from infomap import Infomap, MultilayerNode
        >>> im = Infomap()
        >>> source_multilayer_node = MultilayerNode(layer_id=0, node_id=1)
        >>> target_multilayer_node = MultilayerNode(layer_id=1, node_id=2)
        >>> im.add_multilayer_link(source_multilayer_node, target_multilayer_node)

        Notes
        -----
        This is the full multilayer format that supports both undirected
        and directed links. Infomap will not make any changes to the network.

        Parameters
        ----------
        source_multilayer_node : tuple of int, or MultilayerNode
            If passed a tuple, it should be of the format
            ``(layer_id, node_id)``.
        target_multilayer_node : tuple of int, or MultilayerNode
            If passed a tuple, it should be of the format
            ``(layer_id, node_id)``.
        weight : float, optional

        """
        source_layer_id, source_node_id = source_multilayer_node
        target_layer_id, target_node_id = target_multilayer_node
        return super().addMultilayerLink(
            source_layer_id, source_node_id, target_layer_id, target_node_id, weight
        )

    def add_multilayer_intra_link(
        self, layer_id, source_node_id, target_node_id, weight=1.0
    ):
        """Add an intra-layer link.

        Adds a link within a layer in a multilayer network.

        Examples
        --------

        >>> from infomap import Infomap
        >>> im = Infomap()
        >>> im.add_multilayer_intra_link(1, 1, 2)
        >>> im.add_multilayer_intra_link(1, 2, 3)
        >>> im.add_multilayer_intra_link(2, 1, 3)
        >>> im.add_multilayer_intra_link(2, 3, 4)

        Notes
        -----
        This multilayer format requires a directed network, so if
        the directed flag is not present, it will add all links
        also in their opposite direction to transform the undirected
        input to directed. If no inter-layer links are added, Infomap
        will simulate those by relaxing the random walker's constraint
        to its current layer. The final state network will be generated
        on run, which will clear the temporary data structure that holds
        the provided intra-layer links.

        Parameters
        ----------
        layer_id : int
        source_node_id : int
        target_node_id : int
        weight : float, optional


        """
        if self.num_nodes != 0:
            raise RuntimeError(
                """Using the multilayer intra/inter api and explicitly adding nodes using add_node is unsupported.
If you want to set node names, use set_name."""
            )

        return super().addMultilayerIntraLink(
            layer_id, source_node_id, target_node_id, weight
        )

    def add_multilayer_inter_link(
        self, source_layer_id, node_id, target_layer_id, weight=1.0
    ):
        """Add an inter-layer link.

        Adds a link between two layers in a multilayer network.
        The link is specified through a shared physical node, but
        that jump will not be recorded so Infomap will spread out
        this link to the next possible steps for the random walker
        in the target layer.


        Notes
        -----
        This multilayer format requires a directed network, so if
        the directed flag is not present, it will add all links
        also in their opposite direction to transform the undirected
        input to directed. If no inter-layer links are added, Infomap
        will simulate these by relaxing the random walker's constraint
        to its current layer. The final state network will be generated
        on run, which will clear the temporary data structure that holds
        the provided inter-layer links.

        Examples
        --------

        >>> from infomap import Infomap
        >>> im = Infomap()
        >>> im.add_multilayer_inter_link(1, 1, 2)
        >>> im.add_multilayer_inter_link(1, 2, 2)
        >>> im.add_multilayer_inter_link(2, 1, 1)
        >>> im.add_multilayer_inter_link(2, 3, 1)

        Parameters
        ----------
        source_layer_id : int
        node_id : int
        target_layer_id : int
        weight : float, optional

        """
        if self.num_nodes != 0:
            raise RuntimeError(
                """Using the multilayer intra/inter api and explicitly adding nodes using add_node is unsupported.
If you want to set node names, use set_name."""
            )

        return super().addMultilayerInterLink(
            source_layer_id, node_id, target_layer_id, weight
        )

    def add_multilayer_links(self, links):
        """Add several multilayer links.

        Examples
        --------

        >>> from infomap import Infomap
        >>> im = Infomap()
        >>> links = (
        ...     ((0, 1), (1, 2)),
        ...     ((0, 3), (1, 2))
        ... )
        >>> im.add_multilayer_links(links)


        See Also
        --------
        add_multilayer_link

        Parameters
        ----------
        links : iterable of tuples
            Iterable of tuples of the form
            ``(source_node, target_node, [weight])``
        """
        for link in links:
            self.add_multilayer_link(*link)

    def remove_multilayer_link(self):
        raise NotImplementedError

    def set_meta_data(self, node_id, meta_category):
        """Set meta data to a node.

        Examples
        --------

        >>> from infomap import Infomap
        >>> im = Infomap(silent=True, num_trials=10)
        >>> im.add_links((
        ...     (1, 2), (1, 3), (2, 3),
        ...     (3, 4),
        ...     (4, 5), (4, 6), (5, 6)
        ... ))
        >>> im.set_meta_data(1, 0)
        >>> im.set_meta_data(2, 0)
        >>> im.set_meta_data(3, 1)
        >>> im.set_meta_data(4, 1)
        >>> im.set_meta_data(5, 0)
        >>> im.set_meta_data(6, 0)
        >>> im.run(meta_data_rate=0)
        >>> im.num_top_modules
        2
        >>> im.run(meta_data_rate=2)
        >>> im.num_top_modules
        3


        Parameters
        ----------
        node_id : int
        meta_category : int
        """
        return self.network.addMetaData(node_id, meta_category)

    def add_networkx_graph(self, g, weight="weight"):
        """Add NetworkX graph

        Examples
        --------

        >>> import networkx as nx
        >>> from infomap import Infomap
        >>> G = nx.Graph([("a", "b"), ("a", "c")])
        >>> im = Infomap(silent=True)
        >>> mapping = im.add_networkx_graph(G)
        >>> mapping
        {0: 'a', 1: 'b', 2: 'c'}
        >>> im.run()
        >>> for node in im.nodes:
        ...     print(node.node_id, node.module_id, node.flow, mapping[node.node_id])
        0 1 0.5 a
        1 1 0.25 b
        2 1 0.25 c

        Notes
        -----
        Transforms non-int labels to unique int ids.
        Assumes that all nodes are of the same type.
        If node type is string, they are added as names
        to Infomap.
        If the NetworkX graph is directed (``nx.DiGraph``), and no flow
        model has been specified in the constructor, this method
        sets the ``directed`` flag to ``True``.

        Parameters
        ----------
        g : nx.Graph
            A NetworkX graph
        weight : str, optional
            Key to lookup link weight in edge data if present.
            Default ``"weight"``.

        Returns
        -------
        dict
            Dict with the node ids as keys and labels as values.
        """
        try:
            first = next(iter(g.nodes))
        except StopIteration:
            return dict()

        # If no flow model has been set, and the graph is directed,
        # set the flow model to directed.
        if not super().flowModelIsSet and g.is_directed():
            super().setDirected(True)

        if isinstance(first, int):
            node_map = {node: node for node in g.nodes}
        else:
            node_map = {label: node for node, label in enumerate(g.nodes)}

        is_string_id = isinstance(first, str)

        for label, node in node_map.items():
            self.add_node(node, name=label if is_string_id else None)

        for source, target, data in g.edges(data=True):
            u, v = node_map[source], node_map[target]
            w = data[weight] if weight is not None and weight in data else 1.0
            self.add_link(u, v, w)

        return {node: label for label, node in node_map.items()}

    # ----------------------------------------
    # Run
    # ----------------------------------------

    @property
    def bipartite_start_id(self):
        """Get or set the bipartite start id.

        Examples
        --------

        >>> from infomap import Infomap
        >>> im = Infomap(silent=True, num_trials=10)
        >>> im.add_node(1, "Left 1")
        >>> im.add_node(2, "Left 2")
        >>> im.bipartite_start_id = 3
        >>> im.add_node(3, "Right 3")
        >>> im.add_node(4, "Right 4")
        >>> im.add_link(1, 3)
        >>> im.add_link(1, 4)
        >>> im.add_link(2, 4)
        >>> im.run()
        >>> tol = 1e-4
        >>> abs(im.codelength - 0.9182958340544896) < tol
        True


        Parameters
        ----------
        start_id : int
            The node id where the second node type starts.

        Returns
        -------
        int
            The node id where the second node type starts.
        """
        return self.network.bipartiteStartId()

    @bipartite_start_id.setter
    def bipartite_start_id(self, start_id):
        super().setBipartiteStartId(start_id)

    @property
    def initial_partition(self):
        """Get or set the initial partition.

        This is a initial configuration of nodes into modules where Infomap
        will start the optimizer.

        Examples
        --------

        >>> from infomap import Infomap
        >>> im = Infomap(silent=True)
        >>> im.add_node(1)
        >>> im.add_node(2)
        >>> im.add_node(3)
        >>> im.add_node(4)
        >>> im.add_link(1, 2)
        >>> im.add_link(1, 3)
        >>> im.add_link(2, 3)
        >>> im.add_link(2, 4)
        >>> im.initial_partition = {
        ...     1: 0,
        ...     2: 0,
        ...     3: 1,
        ...     4: 1
        ... }
        >>> im.run(no_infomap=True)
        >>> tol = 1e-4
        >>> abs(im.codelength - 3.4056390622295662) < tol
        True


        Notes
        -----
        The initial partition is saved between runs.
        If you want to use an initial partition for one run only,
        use ``run(initial_partition=partition)``.

        Parameters
        ----------
        module_ids : dict of int, or None
            Dict with node ids as keys and module ids as values.

        Returns
        -------
        dict of int
            Dict with node ids as keys and module ids as values.
        """
        return super().getInitialPartition()

    @initial_partition.setter
    def initial_partition(self, module_ids):
        if module_ids is None:
            module_ids = {}
        super().setInitialPartition(module_ids)

    @contextmanager
    def _initial_partition(self, partition):
        # Internal use only
        # Store the possiply set initial partition
        # and use the new initial partition for one run only.
        old_partition = self.initial_partition
        try:
            self.initial_partition = partition
            yield
        finally:
            self.initial_partition = old_partition

    def run(
        self,
        args=None,
        initial_partition=None,
        # input
        cluster_data=None,
        no_infomap=False,
        skip_adjust_bipartite_flow=False,
        bipartite_teleportation=False,
        weight_threshold=None,
        include_self_links=None,
        no_self_links=False,
        node_limit=None,
        matchable_multilayer_ids=None,
        assign_to_neighbouring_module=False,
        meta_data=None,
        meta_data_rate=_DEFAULT_META_DATA_RATE,
        meta_data_unweighted=False,
        # output
        tree=False,
        ftree=False,
        clu=False,
        verbosity_level=_DEFAULT_VERBOSITY_LEVEL,
        silent=False,
        out_name=None,
        no_file_output=False,
        clu_level=None,
        output=None,
        hide_bipartite_nodes=False,
        print_all_trials=False,
        # algorithm
        two_level=False,
        flow_model=None,
        directed=None,
        recorded_teleportation=False,
        use_node_weights_as_flow=False,
        to_nodes=False,
        teleportation_probability=_DEFAULT_TELEPORTATION_PROB,
        regularized=False,
        regularization_strength=1.0,
        entropy_corrected=False,
        entropy_correction_strength=1.0,
        markov_time=1.0,
        variable_markov_time=False,
        variable_markov_damping=1.0,
        preferred_number_of_modules=None,
        multilayer_relax_rate=_DEFAULT_MULTILAYER_RELAX_RATE,
        multilayer_relax_limit=-1,
        multilayer_relax_limit_up=-1,
        multilayer_relax_limit_down=-1,
        multilayer_relax_by_jsd=False,
        # accuracy
        seed=_DEFAULT_SEED,
        num_trials=1,
        core_loop_limit=10,
        core_level_limit=None,
        tune_iteration_limit=None,
        core_loop_codelength_threshold=_DEFAULT_CORE_LOOP_CODELENGTH_THRESHOLD,
        tune_iteration_relative_threshold=_DEFAULT_TUNE_ITER_RELATIVE_THRESHOLD,
        fast_hierarchical_solution=None,
        prefer_modular_solution=False,
        inner_parallelization=False,
    ):
        """Run Infomap.

        Parameters
        ----------
        args : str, optional
            String of Infomap arguments.
        initial_partition : dict, optional
            Initial partition to start optimizer from.
            (See ``initial_partition``).
        cluster_data : str, optional
            Provide an initial two-level (clu format) or multi-layer (tree
            format) solution.
        no_infomap : bool, optional
            Don't run the optimizer. Useful to calculate codelength of provided
            cluster data or to print non-modular statistics.
        skip_adjust_bipartite_flow : bool, optional
            Skip distributing all flow from the bipartite nodes to the primary
            nodes.
        bipartite_teleportation : bool, optional
            Teleport like the bipartite flow instead of two-step (unipartite)
            teleportation.
        weight_threshold : float, optional
            Limit the number of links to read from the network. Ignore links
            with less weight than the threshold.
        include_self_links : bool, optional
            Include links with the same source and target node.
            DEPRECATED. Include self links by default now, exclude with no_self_links.
        no_self_links : bool, optional
            Exclude self links in the input network.
        node_limit : int, optional
            Limit the number of nodes to read from the network. Ignore links
            connected to ignored nodes.
        matchable_multilayer_ids : int, optional
            Construct state ids from node and layer ids that are consistent
            across networks for the same max number of layers.
            Set to at least the largest layer id among networks to match.
        assign_to_neighbouring_module : bool, optional
            Assign nodes without module assignments (from ``cluster_data``) to
            the module assignment of a neighbouring node if possible.
        meta_data : str, optional
            Provide meta data (clu format) that should be encoded.
        meta_data_rate : float, optional
            Metadata encoding rate. Default is to encode each step.
        meta_data_unweighted : bool, optional
            Don't weight meta data by node flow.
        tree : bool, optional
            Write a tree file with the modular hierarchy. Automatically enabled
            if no other output is specified.
        ftree : bool, optional
            Write a ftree file with the modular hierarchy including aggregated
            links between (nested) modules.
        clu : bool, optional
            Write a clu file with the top cluster ids for each node.
        verbosity_level : int, optional
            Verbose output on the console.
        silent : bool, optional
            No output on the console.
        out_name : str, optional
            Name for the output files, e.g.
            ``"[output_directory]/[out-name].tree"``
        no_file_output : bool, optional
            Don't write output to file.
        clu_level : int, optional
            For clu output, print modules at specified depth from root.
            Use ``-1`` for bottom level modules.
        output : list(str), optional
            List of output formats.
            Options: clu, tree, ftree, newick, json, csv, network, states.
        hide_bipartite_nodes : bool, optional
            Project bipartite solution to unipartite.
        print_all_trials : bool, optional
            Print all trials to separate files.
        two_level : bool, optional
            Optimize a two-level partition of the network.
            Default (``false``) is multi-level.
        flow_model : str, optional
            Specify flow model.
            Options: undirected, directed, undirdir, outdirdir, rawdir.
        directed : bool, optional
            Assume directed links. Shorthand for ``flow_model="directed"``.
        recorded_teleportation : bool, optional
            If teleportation is used to calculate the flow, also record it
            when minimizing codelength. Default False.
        use_node_weights_as_flow : bool, optional
            Use node weights (from api or after names in Pajek format) as flow,
            normalized to sum to 1.
        to_nodes : bool, optional
            Teleport to nodes instead of to links, assuming uniform node
            weights if no such input data.
        teleportation_probability : float, optional
            Probability of teleporting to a random node or link.
        regularized : bool, optional
            Effectively add a fully connected Bayesian prior network to not overfit
            due to missing links. Implies recorded teleportation.
        regularization_strength : float, optional
            Adjust relative strength of Bayesian prior network with this multiplier.
        entropy_corrected : bool, optional
            Correct for negative entropy bias in small samples (many modules).
        entropy_correction_strength : float, optional
            Increase or decrease the default entropy correction with this multiplier.
        markov_time : float, optional
            Scales link flow to change the cost of moving between modules.
            Higher values results in fewer modules.
        variable_markov_time : bool, optional
            Increase Markov time locally to level out link flow. Reduces risk of
            overpartitioning sparse areas while keeping high resolution in dense areas.
        variable_markov_damping : float, optional
            Exponent for variable Markov time scale. 0 means no rescaling and 1 means
            full rescaling to constant transition flow rate.
        preferred_number_of_modules : int, optional
            Penalize solutions the more they differ from this number.
        multilayer_relax_rate : float, optional
            Probability to relax the constraint to move only in the current
            layer.
        multilayer_relax_limit : int, optional
            Number of neighboring layers in each direction to relax to.
            If negative, relax to any layer.
        multilayer_relax_limit_up : int, optional
            Number of neighboring layers with higher id to relax to.
            If negative, relax to any layer.
        multilayer_relax_limit_down : int, optional
            Number of neighboring layers with lower id to relax to.
            If negative, relax to any layer.
        multilayer_relax_by_jsd : bool, optional
            Relax proportional to the out-link similarity measured by the
            Jensen-Shannon divergence.
        seed : int, optional
            A seed (integer) to the random number generator for reproducible
            results.
        num_trials : int, optional
            Number of outer-most loops to run before picking the best solution.
        core_loop_limit : int, optional
            Limit the number of loops that tries to move each node into the
            best possible module.
        core_level_limit : int, optional
            Limit the number of times the core loops are reapplied on existing
            modular network to search bigger structures.
        tune_iteration_limit : int, optional
            Limit the number of main iterations in the two-level partition
            algorithm. 0 means no limit.
        core_loop_codelength_threshold : float, optional
            Minimum codelength threshold for accepting a new solution in core
            loop.
        tune_iteration_relative_threshold : float, optional
            Set codelength improvement threshold of each new tune iteration to
            ``f`` times the initial two-level codelength.
        fast_hierarchical_solution : int, optional
            Find top modules fast. Use ``2`` to keep all fast levels.
            Use ``3`` to skip recursive part.
        prefer_modular_solution : bool, optional
            Prefer modular solutions even if they are worse than putting all
            nodes in one module.
        inner_parallelization : bool, optional
            Parallelize the inner-most loop for greater speed.
            This may give some accuracy tradeoff.

        See Also
        --------
        initial_partition
        """
        args = _construct_args(
            args,
            cluster_data=cluster_data,
            no_infomap=no_infomap,
            skip_adjust_bipartite_flow=skip_adjust_bipartite_flow,
            bipartite_teleportation=bipartite_teleportation,
            weight_threshold=weight_threshold,
            include_self_links=include_self_links,
            no_self_links=no_self_links,
            node_limit=node_limit,
            matchable_multilayer_ids=matchable_multilayer_ids,
            assign_to_neighbouring_module=assign_to_neighbouring_module,
            meta_data=meta_data,
            meta_data_rate=meta_data_rate,
            meta_data_unweighted=meta_data_unweighted,
            tree=tree,
            ftree=ftree,
            clu=clu,
            verbosity_level=verbosity_level,
            silent=silent,
            out_name=out_name,
            no_file_output=no_file_output,
            clu_level=clu_level,
            output=output,
            hide_bipartite_nodes=hide_bipartite_nodes,
            print_all_trials=print_all_trials,
            two_level=two_level,
            flow_model=flow_model,
            directed=directed,
            recorded_teleportation=recorded_teleportation,
            use_node_weights_as_flow=use_node_weights_as_flow,
            to_nodes=to_nodes,
            teleportation_probability=teleportation_probability,
            regularized=regularized,
            regularization_strength=regularization_strength,
            entropy_corrected=entropy_corrected,
            entropy_correction_strength=entropy_correction_strength,
            markov_time=markov_time,
            variable_markov_time=False,
            variable_markov_damping=1.0,
            preferred_number_of_modules=preferred_number_of_modules,
            multilayer_relax_rate=multilayer_relax_rate,
            multilayer_relax_limit=multilayer_relax_limit,
            multilayer_relax_limit_up=multilayer_relax_limit_up,
            multilayer_relax_limit_down=multilayer_relax_limit_down,
            multilayer_relax_by_jsd=multilayer_relax_by_jsd,
            seed=seed,
            num_trials=num_trials,
            core_loop_limit=core_loop_limit,
            core_level_limit=core_level_limit,
            tune_iteration_limit=tune_iteration_limit,
            core_loop_codelength_threshold=core_loop_codelength_threshold,
            tune_iteration_relative_threshold=tune_iteration_relative_threshold,
            fast_hierarchical_solution=fast_hierarchical_solution,
            prefer_modular_solution=prefer_modular_solution,
            inner_parallelization=inner_parallelization,
        )

        if initial_partition:
            with self._initial_partition(initial_partition):
                return super().run(args)

        return super().run(args)

    # ----------------------------------------
    # Iterators/Output
    # ----------------------------------------

    def get_modules(self, depth_level=1, states=False):
        """Get a dict with node ids as keys and module ids as values for a given depth in the hierarchical tree.

        ::

            Level                            Root

              0                               ┌─┐
                                    ┌─────────┴─┴────────┐
                                    │                    │
                                    │                    │
                                    │                    │
                              Path  │  Module      Path  │  Module
              1                  1 ┌┼┐ 1              2 ┌┼┐ 2
                               ┌───┴─┴───┐          ┌───┴─┴───┐
                               │         │          │         │
                               │         │          │         │
                               │         │          │         │
                               │         │          │         │
              2             1 ┌┼┐ 1   2 ┌┼┐ 2    1 ┌┼┐ 3   2 ┌┼┐ 3
                          ┌───┴─┴───┐   └─┴────┐   └─┘       └─┘
                          │         │          │
                          │         │          │    ▲         ▲
                          │         │          │    └────┬────┘
                          │         │          │         │
              3        1 ┌┼┐     2 ┌┼┐      1 ┌┼┐
                         └─┘       └─┘        └─┘  ◄───  Leaf-nodes


        Path to the left of the nodes. Depth dependent module ids to the right.
        The five leaf-nodes are network-nodes. All other tree-nodes are modules.

        For example:

        The left-most node on level 3 has path 1:1:1 and belong to module 1 on level 1.

        The right-most node on level 2 has path 2:2 and belong to module 2 on level 1
        which is renamed to module 3 on level 2 as we have more modules in total on this level.

        Assuming the nodes are labelled 1-5 from left to right, then the first three nodes
        are in module 1, and the last two nodes are in module 2::

            > im.get_modules(depth_level=1)
            {1: 1, 2: 1, 3: 1, 4: 2, 5: 2}

        However, at level 2, the first two nodes are in module 1, the third node in module 2,
        and the last two nodes are in module 3::

            > im.get_modules(depth_level=2)
            {1: 1, 2: 1, 3: 2, 4: 3, 5: 3}


        Examples
        --------

        >>> from infomap import Infomap
        >>> im = Infomap(silent=True)
        >>> im.read_file("twotriangles.net")
        >>> im.run()
        >>> im.get_modules()
        {1: 1, 2: 1, 3: 1, 4: 2, 5: 2, 6: 2}

        >>> from infomap import Infomap
        >>> im = Infomap(silent=True)
        >>> im.read_file("states.net")
        >>> im.run()
        >>> im.get_modules(states=True)
        {1: 1, 2: 1, 3: 1, 4: 2, 5: 2, 6: 2}


        Notes
        -----
        In a higher-order network, a physical node (defined by ``node_id``)
        may partially exist in multiple modules. However, the ``node_id``
        can not exist multiple times as a key in the node-to-module map,
        so only one occurrence of a physical node will be retrieved.
        To get all states, use ``get_modules(states=True)``.

        Parameters
        ----------
        depth_level : int, optional
            The level in the hierarchical tree. Set to ``1`` (default) to
            return the top modules (coarsest level). Set to ``2`` for second
            coarsest level etc. Set to ``-1`` to return the bottom level
            modules (finest level). Default ``1``.
        states : bool, optional
            For higher-order networks, if ``states`` is True, it will return
            state node ids. Otherwise it will return physical node ids,
            merging state nodes with same ``node_id`` if they are in the same
            module. Note that the same physical node may end up on different
            paths in the tree.
            Default ``false``.

        Returns
        -------
        dict of int
            Dict with node ids as keys and module ids as values.
        """
        return super().getModules(depth_level, states)

    def get_multilevel_modules(self, states=False):
        """Get a dict with node ids as keys and a tuple of module ids as values.
        Each position in the tuple corresponds to a depth in the hierarchical tree,
        with the first level being the top level.

        See Also
        --------
        get_modules


        Examples
        --------

        >>> from infomap import Infomap
        >>> im = Infomap(silent=True, num_trials=10)
        >>> im.read_file("ninetriangles.net")
        >>> im.run()
        >>> for modules in sorted(im.get_multilevel_modules().values()):
        ...     print(modules)
        (1, 1)
        (1, 1)
        (1, 1)
        (1, 2)
        (1, 2)
        (1, 2)
        (1, 3)
        (1, 3)
        (1, 3)
        (2, 4)
        (2, 4)
        (2, 4)
        (2, 5)
        (2, 5)
        (2, 5)
        (2, 6)
        (2, 6)
        (2, 6)
        (3, 7)
        (3, 7)
        (3, 7)
        (4, 8)
        (4, 8)
        (4, 8)
        (5, 9)
        (5, 9)
        (5, 9)

        >>> from infomap import Infomap
        >>> im = Infomap(silent=True)
        >>> im.read_file("states.net")
        >>> im.run()
        >>> for node, modules in im.get_multilevel_modules(states=True).items():
        ...     print(node, modules)
        1 (1,)
        2 (1,)
        3 (1,)
        4 (2,)
        5 (2,)
        6 (2,)


        Notes
        -----
        In a higher-order network, a physical node (defined by ``node_id``)
        may partially exist in multiple modules. However, the ``node_id``
        can not exist multiple times as a key in the node-to-module map,
        so only one occurrence of a physical node will be retrieved.
        To get all states, use ``get_multilevel_modules(states=True)``.

        Parameters
        ----------
        states : bool, optional
            For higher-order networks, if ``states`` is True, it will return
            state node ids. Otherwise it will return physical node ids,
            merging state nodes with same ``node_id`` if they are in the same
            module. Note that the same physical node may end up on different
            paths in the tree.
            Default ``false``.

        Returns
        -------
        dict of list of int
            Dict with node ids as keys and tuple of module ids as values.
        """
        return super().getMultilevelModules(states)

    @property
    def modules(self):
        """A view of the top-level modules, mapping
        ``node_id`` to ``module_id``.

        Notes
        -----
        In a higher-order network, a physical node (defined by ``node_id``)
        may partially exist in multiple modules. However, the ``node_id``
        can not exist multiple times as a key in the node-to-module map,
        so only one occurrence of a physical node will be retrieved.
        To get all states, use ``get_modules(states=True)``.

        Examples
        --------
        >>> from infomap import Infomap
        >>> im = Infomap(silent=True, num_trials=5)
        >>> im.read_file("twotriangles.net")
        >>> im.run()
        >>> for node_id, module_id in im.modules:
        ...     print(node_id, module_id)
        ...
        1 1
        2 1
        3 1
        4 2
        5 2
        6 2

        See Also
        --------
        get_modules

        Yields
        -------
        tuple of int, int
            An iterator of ``(node_id, module_id)`` pairs.
        """
        return self.get_modules(depth_level=1, states=False).items()

    @property
    def multilevel_modules(self):
        """A view of the multilevel modules, mapping
        ``node_id`` to a tuple of ``module_id``.

        Notes
        -----
        In a higher-order network, a physical node (defined by ``node_id``)
        may partially exist in multiple modules. However, the ``node_id``
        can not exist multiple times as a key in the node-to-module map,
        so only one occurrence of a physical node will be retrieved.
        To get all states, use ``get_multilevel_modules(states=True)``.

        See Also
        --------
        get_multilevel_modules

        Yields
        -------
        tuple of (int, tuple of int)
            An iterator of ``(node_id, (module_ids...)`` pairs.
        """
        return self.get_multilevel_modules().items()

    def get_tree(self, depth_level=1, states=False):
        """A view of the hierarchical tree, iterating
        over the modules as well as the leaf-nodes.

        Parameters
        ----------
        depth_level : int, optional
            The module level returned by ``iterator.module_id``. Set to ``1``
            (default) to return the top modules (coarsest level). Set to ``2``
            for second coarsest level etc. Set to ``-1`` to return the bottom
            level modules (finest level).
        states : bool, optional
            For higher-order networks, if ``states`` is True, it will iterate
            over state nodes. Otherwise it will iterate over physical nodes,
            merging state nodes with same ``node_id`` if they are in the same
            module. Note that the same physical node may end up on different
            paths in the tree.
            Default ``false``.

        Notes
        ----
        For higher-order networks, each node is represented by a set of state
        nodes with the same ``node_id``, where each state node represents a
        different constraint on the random walker. This enables
        overlapping modules, where state nodes with the same ``node_id`` end up
        in different modules. However, the state nodes with the same
        ``node_id`` within each module are only visible as one (partial)
        physical node (if ``states = False``).

        Returns
        -------
        InfomapIterator or InfomapIteratorPhysical
            An iterator over each node in the tree, depth first from the root
        """
        if self.have_memory and not states:
            return super().iterTreePhysical(depth_level)
        return super().iterTree(depth_level)

    def get_nodes(self, depth_level=1, states=False):
        """A view of the nodes in the hierarchical tree, iterating depth first
        from the root.

        Parameters
        ----------
        depth_level : int, optional
            The module level returned by ``iterator.module_id``. Set to ``1``
            (default) to return the top modules (coarsest level). Set to ``2``
            for second coarsest level etc. Set to ``-1`` to return the bottom
            level modules (finest level). Default ``1``.
        states : bool, optional
            For higher-order networks, if ``states`` is True, it will iterate
            over state nodes. Otherwise it will iterate over physical nodes,
            merging state nodes with same ``node_id`` if they are in the same
            module. Note that the same physical node may end up on different
            paths in the tree.
            See notes on ``physical_tree``. Default ``false``.

        Notes
        -----
        For higher-order networks, each node is represented by a set of state
        nodes with the same ``node_id``, where each state node represents a
        different constraint on the random walker. This enables
        overlapping modules, where state nodes with the same ``node_id`` end up
        in different modules. However, the state nodes with the same
        ``node_id`` within each module are only visible as one (partial)
        physical node (if ``states = False``).

        Returns
        -------
        InfomapIterator or InfomapIteratorPhysical
            An iterator over each node in the tree, depth first from the root
        """

        class LeafIterWrapper:
            def __init__(self, treeIterator):
                self.it = treeIterator

            def __iter__(self):
                return self

            def __next__(self):
                self.it.stepForward()
                while not self.it.isEnd() and not self.it.isLeaf():
                    self.it.stepForward()
                if not self.it.isEnd():
                    return self.it
                raise StopIteration

        if self.have_memory and not states:
            # super().iterLeafNodesPhysical(depth_level) is unreliable in python
            return LeafIterWrapper(super().iterTreePhysical(depth_level))
        return super().iterLeafNodes(depth_level)

    @property
    def tree(self):
        """A view of the hierarchical tree, iterating
        over the modules as well as the leaf-nodes.

        Convinience method for ``get_tree(depth_level=1, states=True)``.

        See Also
        --------
        get_tree
        InfomapIterator

        Returns
        -------
        InfomapIterator
            An iterator over each node in the tree, depth first from the root
        """
        return self.get_tree(depth_level=1, states=True)

    @property
    def physical_tree(self):
        """A view of the hierarchical tree, iterating
        over the modules as well as the leaf-nodes.
        All state nodes with the same ``node_id`` are merged to one physical node.

        Convinience method for ``get_tree(depth_level=1, states=False)``.

        See Also
        --------
        get_tree
        InfomapIteratorPhysical

        Returns
        -------
        InfomapIteratorPhysical
            An iterator over each physical node in the tree, depth first from
            the root
        """
        return self.get_tree(depth_level=1, states=False)

    @property
    def leaf_modules(self):
        """A view of the leaf modules, i.e. the bottom modules containing leaf nodes.

        See Also
        --------
        get_modules
        InfomapLeafModuleIterator

        Returns
        -------
        InfomapLeafModuleIterator
            An iterator over each leaf module in the tree, depth first from the
            root
        """
        return super().iterLeafModules()

    @property
    def nodes(self):
        """A view of the nodes in the hierarchical tree, iterating depth first
        from the root.

        Convinience method for ``get_nodes(depth_level=1, states=True)``.

        See Also
        --------
        get_nodes
        InfomapLeafIterator

        Returns
        -------
        InfomapLeafIterator
            An iterator over each leaf node in the tree, depth first from the
            root
        """
        return self.get_nodes(depth_level=1, states=True)

    @property
    def physical_nodes(self):
        """A view of the nodes in the hierarchical tree, iterating depth first
        from the root.
        All state nodes with the same ``node_id`` are merged to one physical node.

        Convinience method for ``get_nodes(depth_level=1, states=False)``.

        See Also
        --------
        get_nodes
        InfomapLeafIteratorPhysical

        Returns
        -------
        InfomapLeafIteratorPhysical
            An iterator over each physical leaf node in the tree, depth first
            from the root
        """
        return self.get_nodes(depth_level=1, states=False)

    def get_dataframe(self, columns=["path", "flow", "name", "node_id"]):
        """Get a Pandas DataFrame with the selected columns.

        Examples
        --------

        >>> from infomap import Infomap
        >>> im = Infomap(silent=True)
        >>> im.read_file("twotriangles.net")
        >>> im.run()
        >>> im.get_dataframe()
             path      flow name  node_id
        0  (1, 1)  0.214286    C        3
        1  (1, 2)  0.142857    A        1
        2  (1, 3)  0.142857    B        2
        3  (2, 1)  0.214286    D        4
        4  (2, 2)  0.142857    E        5
        5  (2, 3)  0.142857    F        6
        >>> im.get_dataframe(columns=["node_id", "module_id"])
           node_id  module_id
        0        3          1
        1        1          1
        2        2          1
        3        4          2
        4        5          2
        5        6          2


        See Also
        --------
        InfoNode
        InfomapLeafIterator
        InfomapLeafIteratorPhysical

        Parameters
        ---------
        columns : list(str)
            A list of columns that should be extracted from each node.
            Must be available as an attribute of ``InfoNode``,
            ``InfomapLeafIterator`` (for state nodes),
            or ``InfomapLeafIteratorPhysical``.
            One exception to this is ``"name"`` which is looked up internally.
            Default ``["path", "flow", "name", "node_id"]``.

        Raises
        ------
        ImportError
            If the pandas package is not available.
        AttributeError
            If a column name is not available as an ``InfoNode`` attribute.

        Returns
        -------
        pandas.DataFrame
            A DataFrame containing the selected columns.
        """
        if pandas is None:
            raise ImportError("Cannot import package 'pandas'")

        return pandas.DataFrame(
            [
                [
                    getattr(node, attr)
                    if attr != "name"
                    else self.get_name(node.node_id, default=node.node_id)
                    for attr in columns
                ]
                for node in self.nodes
            ],
            columns=columns,
        )

    def get_name(self, node_id, default=None):
        """Get the name of a node.

        Notes
        -----
        If the node name is an empty string,
        the ``default`` will be returned.

        See Also
        --------
        set_name
        names

        Parameters
        ----------
        node_id : int
        default : str, optional
            The return value if the node name is missing, default ``None``

        Returns
        -------
        str
            The node name if it exists, else the ``default``.
        """
        name = super().getName(node_id)
        if name == "":
            return default
        return name

    def get_names(self):
        """Get all node names.

        See Also
        --------
        names
        get_name

        Returns
        -------
        dict of string
            A dict with node ids as keys and node names as values.
        """
        return super().getNames()

    @property
    def names(self):
        """Get all node names.

        Short-hand for ``get_names``.

        See Also
        --------
        get_names
        get_name

        Returns
        -------
        dict of string
            A dict with node ids as keys and node names as values.
        """
        return super().getNames()

    def get_links(self, data="weight"):
        """A view of the currently assigned links and their weights or flow.

        The sources and targets are state ids when we have a
        state or multilayer network.

        Examples
        --------

        >>> from infomap import Infomap
        >>> im = Infomap(silent=True)
        >>> im.read_file("twotriangles.net")
        >>> im.run()
        >>> for link in im.get_links():
        ...     print(link)
        (1, 2, 1.0)
        (1, 3, 1.0)
        (2, 3, 1.0)
        (3, 4, 1.0)
        (4, 5, 1.0)
        (4, 6, 1.0)
        (5, 6, 1.0)
        >>> for link in im.get_links(data="flow"):
        ...     print(link)
        (1, 2, 0.14285714285714285)
        (1, 3, 0.14285714285714285)
        (2, 3, 0.14285714285714285)
        (3, 4, 0.14285714285714285)
        (4, 5, 0.14285714285714285)
        (4, 6, 0.14285714285714285)
        (5, 6, 0.14285714285714285)


        See Also
        --------
        links
        flow_links

        Parameters
        ----------
        data : str
            The kind of data to return, one of ``"weight"`` or ``"flow"``.
            Default ``"weight"``.

        Returns
        -------
        tuple of int, int, float
            An iterator of source, target, weight/flow tuples.
        """
        if data not in ("weight", "flow"):
            raise RuntimeError('data must one of "weight" or "flow"')

        return (
            (source, target, value)
            for (source, target), value in self.getLinks(data != "weight").items()
        )

    @property
    def links(self):
        """A view of the currently assigned links and their weights.

        The sources and targets are state ids when we have a
        state or multilayer network.

        Examples
        --------

        >>> from infomap import Infomap
        >>> im = Infomap(silent=True)
        >>> im.read_file("twotriangles.net")
        >>> im.run()
        >>> for link in im.links:
        ...     print(link)
        (1, 2, 1.0)
        (1, 3, 1.0)
        (2, 3, 1.0)
        (3, 4, 1.0)
        (4, 5, 1.0)
        (4, 6, 1.0)
        (5, 6, 1.0)


        See Also
        --------
        flow_links

        Returns
        -------
        tuple of int, int, float
            An iterator of source, target, weight tuples.
        """
        return self.get_links()

    @property
    def flow_links(self):
        """A view of the currently assigned links and their flow.

        The sources and targets are state ids when we have a
        state or multilayer network.

        Examples
        --------

        >>> from infomap import Infomap
        >>> im = Infomap(silent=True)
        >>> im.read_file("twotriangles.net")
        >>> im.run()
        >>> for link in im.flow_links:
        ...     print(link)
        (1, 2, 0.14285714285714285)
        (1, 3, 0.14285714285714285)
        (2, 3, 0.14285714285714285)
        (3, 4, 0.14285714285714285)
        (4, 5, 0.14285714285714285)
        (4, 6, 0.14285714285714285)
        (5, 6, 0.14285714285714285)


        See Also
        --------
        links

        Returns
        -------
        tuple of int, int, float
            An iterator of source, target, flow tuples.
        """
        return self.get_links(data="flow")

    @property
    def network(self):
        """Get the internal network."""
        return super().network()

    # ----------------------------------------
    # num_* getters
    # ----------------------------------------

    @property
    def num_nodes(self):
        """The number of state nodes if we have a higher order network, or the
        number of physical nodes.

        See Also
        --------
        num_physical_nodes

        Returns
        -------
        int
            The number of nodes
        """
        return self.network.numNodes()

    @property
    def num_links(self):
        """The number of links.

        Returns
        -------
        int
            The number of links
        """
        return self.network.numLinks()

    @property
    def num_physical_nodes(self):
        """The number of physical nodes.

        See Also
        --------
        num_nodes

        Returns
        -------
        int
            The number of nodes
        """
        return self.network.numPhysicalNodes()

    @property
    def num_top_modules(self):
        """Get the number of top modules in the tree

        Returns
        -------
        int
            The number of top modules
        """
        return super().numTopModules()

    @property
    def num_non_trivial_top_modules(self):
        """Get the number of non-trivial top modules in the tree

        A trivial module is a module with either one or all nodes within.

        Returns
        -------
        int
            The number of non-trivial top modules
        """
        return super().numNonTrivialTopModules()

    @property
    def num_leaf_modules(self):
        """Get the number of leaf modules in the tree

        Returns
        -------
        int
            The number of leaf modules
        """
        num_leaf_modules = 0
        for _ in self.leaf_modules:
            num_leaf_modules += 1
        return num_leaf_modules

    def get_effective_num_modules(self, depth_level=1):
        """The flow weighted effective number of modules.

        Measured as the perplexity of the module flow distribution.

        Parameters
        ----------
        depth_level : int, optional
            The module level returned by ``iterator.depth``.
            Set to ``1`` (default) to return the top modules (coarsest level).
            Set to ``2`` for second coarsest level etc. Set to ``-1`` to return
            the bottom level modules (finest level).

        Returns
        -------
        float
            The effective number of modules
        """
        return perplexity(
            [
                module.flow
                for module in self.get_tree(depth_level)
                if depth_level == -1
                and module.is_leaf_module
                or module.depth == depth_level
            ]
        )

    @property
    def effective_num_top_modules(self):
        """The flow weighted effective number of top modules.

        Measured as the perplexity of the module flow distribution.

        Returns
        -------
        float
            The effective number of top modules
        """
        return self.get_effective_num_modules(depth_level=1)

    @property
    def effective_num_leaf_modules(self):
        """The flow weighted effective number of leaf modules.

        Measured as the perplexity of the module flow distribution.

        Returns
        -------
        float
            The effective number of top modules
        """
        return self.get_effective_num_modules(depth_level=-1)

    @property
    def max_depth(self):
        """Get the max depth of the hierarchical tree.

        Returns
        -------
        int
            The max depth
        """
        return super().maxTreeDepth()

    @property
    def num_levels(self):
        """Get the max depth of the hierarchical tree.
        Alias of ``max_depth``.

        See Also
        --------
        max_depth

        Returns
        -------
        int
            The max depth
        """
        return self.max_depth

    @property
    def have_memory(self):
        """Returns true for multilayer and memory networks.

        Returns
        -------
        bool
            True if the network is a multilayer or memory network.
        """
        return super().haveMemory()

    # ----------------------------------------
    # Codelength
    # ----------------------------------------

    @property
    def codelength(self):
        """Get the total (hierarchical) codelength.

        See Also
        --------
        index_codelength
        module_codelength

        Returns
        -------
        float
            The codelength
        """
        return super().codelength()

    @property
    def codelengths(self):
        """Get the total (hierarchical) codelength for each trial.

        See Also
        --------
        codelength

        Returns
        -------
        tuple of float
            The codelengths for each trial
        """
        return super().codelengths()

    @property
    def index_codelength(self):
        """Get the two-level index codelength.

        See Also
        --------
        codelength
        module_codelength

        Returns
        -------
        float
            The two-level index codelength
        """
        return super().getIndexCodelength()

    @property
    def module_codelength(self):
        """Get the total codelength of the modules.

        The module codelength is defined such that
        ``codelength = index_codelength + module_codelength``

        For a hierarchical solution, the module codelength
        is the sum of codelengths for each top module.

        See Also
        --------
        codelength
        index_codelength

        Returns
        -------
        float
            The module codelength
        """
        return super().getModuleCodelength()

    @property
    def one_level_codelength(self):
        """Get the one-level codelength.

        See Also
        --------
        codelength

        Returns
        -------
        float
            The one-level codelength
        """
        return super().getOneLevelCodelength()

    @property
    def relative_codelength_savings(self):
        """Get the relative codelength savings.

        This is defined as the reduction in codelength
        relative to the non-modular one-level solution::

            S_L = 1 - L / L_1

        where ``L`` is the ``codelength`` and ``L_1``
        the ``one_level_codelength``.

        See Also
        --------
        codelength
        one_level_codelength

        Returns
        -------
        float
            The relative codelength savings
        """
        return super().getRelativeCodelengthSavings()

    @property
    def entropy_rate(self):
        """Get the entropy rate of the network.

        The entropy rate is an indication of the sparsity of a network.
        A higher entropy rate corresponds to a densely connected network.

        Notes
        -----
        This value is only accessable after running the optimizer (``im.run()``).

        Examples
        --------

        >>> from infomap import Infomap
        >>> im = Infomap(silent=True, no_infomap=True)
        >>> im.read_file("twotriangles.net")
        >>> im.run()
        >>> f"{im.entropy_rate:.5f}"
        '1.25070'


        Returns
        -------
        float
            The entropy rate
        """
        return super().getEntropyRate()

    @property
    def meta_codelength(self):
        """Get the meta codelength.

        This is the meta entropy times the meta data rate.

        See Also
        --------
        meta_entropy

        Returns
        -------
        float
            The meta codelength
        """
        return super().getMetaCodelength()

    @property
    def meta_entropy(self):
        """Get the meta entropy (unweighted by meta data rate).

        See Also
        --------
        meta_codelength

        Returns
        -------
        float
            The meta entropy
        """
        return super().getMetaCodelength(True)

    # ----------------------------------------
    # Write Results
    # ----------------------------------------

    def write(self, filename, *args, **kwargs):
        """Write results to file.

        Raises
        ------
        NotImplementedError
            If the file format is not supported.

        Parameters
        ----------
        filename : str
            The filename.
        """
        _, ext = os.path.splitext(filename)

        # remove the dot
        ext = ext[1:]

        if ext == "ftree":
            ext = "flow_tree"
        elif ext == "nwk":
            ext = "newick"
        elif ext == "net":
            ext = "pajek"

        writer = "write_{}".format(ext)

        if hasattr(self, writer):
            return getattr(self, writer)(filename, *args, **kwargs)

        raise NotImplementedError("No method found for writing {} files".format(ext))

    def write_clu(self, filename, states=False, depth_level=1):
        """Write result to a clu file.

        See Also
        --------
        write_tree
        write_flow_tree

        Parameters
        ----------
        filename : str
        states : bool, optional
            If the state nodes should be included. Default ``False``.
        depth_level : int, optional
            The depth in the hierarchical tree to write.
        """
        return self.writeClu(filename, states, depth_level)

    def write_tree(self, filename, states=False):
        """Write result to a tree file.

        See Also
        --------
        write_clu
        write_flow_tree

        Parameters
        ----------
        filename : str
        states : bool, optional
            If the state nodes should be included. Default ``False``.
        """
        return self.writeTree(filename, states)

    def write_flow_tree(self, filename, states=False):
        """Write result to a ftree file.

        See Also
        --------
        write_clu
        write_tree

        Parameters
        ----------
        filename : str
        states : bool, optional
            If the state nodes should be included. Default ``False``.
        """
        return self.writeFlowTree(filename, states)

    def write_newick(self, filename, states=False):
        """Write result to a Newick file.

        See Also
        --------
        write_clu
        write_tree

        Parameters
        ----------
        filename : str
        states : bool, optional
            If the state nodes should be included. Default ``False``.
        """
        return self.writeNewickTree(filename, states)

    def write_json(self, filename, states=False):
        """Write result to a JSON file.

        See Also
        --------
        write_clu
        write_tree

        Parameters
        ----------
        filename : str
        states : bool, optional
            If the state nodes should be included. Default ``False``.
        """
        return self.writeJsonTree(filename, states)

    def write_csv(self, filename, states=False):
        """Write result to a CSV file.

        See Also
        --------
        write_clu
        write_tree

        Parameters
        ----------
        filename : str
        states : bool, optional
            If the state nodes should be included. Default ``False``.
        """
        return self.writeCsvTree(filename, states)

    def write_state_network(self, filename):
        """Write internal state network to file.

        See Also
        --------
        write_pajek

        Parameters
        ----------
        filename : str
        """
        return self.network.writeStateNetwork(filename)

    def write_pajek(self, filename, flow=False):
        """Write network to a Pajek file.

        See Also
        --------
        write_state_network

        Parameters
        ----------
        filename : str
        flow : bool, optional
            If the flow should be included. Default ``False``.
        """
        return self.network.writePajekNetwork(filename, flow)


def main():
    import sys

    args = " ".join(sys.argv[1:])
    conf = Config(args, True)
    im = Infomap(conf)
    im.run()


if __name__ == "__main__":
    main()
