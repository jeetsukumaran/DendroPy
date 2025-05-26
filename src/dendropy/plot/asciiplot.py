#! /usr/bin/env python
# -*- coding: utf-8 -*-

from dendropy.utility import terminal

class AsciiTreePlot(object):
    class NullEdgeLengthError(ValueError):
        def __init__(self, *args, **kwargs):
            ValueError.__init__(self, *args, **kwargs)

    def __init__(self, **kwargs):
        """

        Keyword Arguments
        -----------------

        plot_metric : str
            A string which specifies how branches should be scaled, one of:
            'age' (distance from tips), 'depth' (distance from root),
            'level' (number of branches from root) or 'length' (edge
            length/weights).
        show_internal_node_labels : bool
            Whether or not to write out internal node labels.
        leaf_spacing_factor : int
            Positive integer: number of rows between each leaf.
        width : int
            Force a particular display width, in terms of number of columns.
        node_label_compose_fn : function object
            A function that takes a Node object as an argument and returns
            the string to be used to display it.

        """
        self.plot_metric = kwargs.pop("plot_metric", "depth")
        self.show_external_node_labels = kwargs.pop("show_external_node_labels", True)
        self.show_internal_node_labels = kwargs.pop("show_internal_node_labels", False)
        self.leaf_spacing_factor = kwargs.pop("leaf_spacing_factor", 2)
        self.width = kwargs.pop("width", None)
        self.display_width = kwargs.pop("display_width", self.width)  # legacy
        self.compose_node = kwargs.pop("node_label_compose_fn", None)
        if self.compose_node is None:
            self.compose_node = self.default_compose_node
        if kwargs:
            raise TypeError("Unrecognized or unsupported arguments: {}".format(kwargs))

    def default_compose_node(self, node):
        if node.taxon is not None and node.taxon.label is not None:
            return node.taxon.label
        elif node.label is not None:
            return node.label
        else:
            return "@"

    def reset(self):
        self.grid = []
        self.node_row = {}
        self.node_col = {}
        self.node_offset = {}
        self.current_leaf_row = 0
        self.node_label_map = {}

    def _calc_node_offsets(self, tree):
        if self.plot_metric == "age" or self.plot_metric == "depth":

            for nd in tree.postorder_node_iter():
                cnds = nd.child_nodes()
                if self.plot_metric == "depth":  # 'number of branchings from tip'
                    if len(cnds) == 0:
                        curr_node_offset = 0.0
                    else:
                        depths = [self.node_offset[v] for v in cnds]
                        curr_node_offset = max(depths) + 1
                elif self.plot_metric == "age":  # 'sum of edge weights from tip'
                    # note: no enforcement of ultrametricity!
                    if len(cnds) == 0:
                        curr_node_offset = 0.0
                    else:
                        if cnds[0].edge.length is not None:
                            curr_node_offset = (
                                self.node_offset[cnds[0]] + cnds[0].edge.length
                            )
                else:
                    raise ValueError(
                        "Unrecognized plot metric '%s' (must be one of: 'age', 'depth',"
                        " 'level', or 'length')"
                        % self.plot_metric
                    )
                self.node_offset[nd] = curr_node_offset
            flipped_origin = max(self.node_offset.values())
            for nd in self.node_offset:
                self.node_offset[nd] = flipped_origin - self.node_offset[nd]
        else:
            for nd in tree.preorder_node_iter():
                if self.plot_metric == "level":  # 'number of branchings from root'
                    curr_edge_len = 1
                elif self.plot_metric == "length":  # 'sum of edge weights from root'
                    if nd.edge.length is not None:
                        curr_edge_len = nd.edge.length
                    else:
                        curr_edge_len = 0
                else:
                    raise ValueError(
                        "Unrecognized plot metric '%s' (must be one of: 'age', 'depth',"
                        " 'level', or 'length')"
                        % self.plot_metric
                    )
                if nd._parent_node is None:
                    self.node_offset[nd] = curr_edge_len
                else:
                    self.node_offset[nd] = (
                        curr_edge_len + self.node_offset[nd._parent_node]
                    )

    def draw(self, tree, dest):
        dest.write(self.compose(tree))

    def get_label_for_node(self, node):
        try:
            return self.node_label_map[node]
        except KeyError:
            if node._child_nodes and self.show_internal_node_labels:
                label = self.compose_node(node)
            elif not node._child_nodes and self.show_external_node_labels:
                label = self.compose_node(node)
            else:
                label = ""
            self.node_label_map[node] = label
            return label

    def compose(self, tree):
        self.reset()
        if self.display_width is None:
            display_width = terminal.terminal_width() - 1
        else:
            display_width = self.display_width
        max_label_len = max(
            [len(self.get_label_for_node(i)) for i in tree.leaf_node_iter()]
        )
        if max_label_len <= 0:
            max_label_len = 0
        # effective_display_width = display_width - max_label_len - len(tree.internal_nodes) - 1
        effective_display_width = display_width - max_label_len - 1
        self._calc_node_offsets(tree)
        widths = [
            self.node_offset[i]
            for i in tree.leaf_node_iter()
            if self.node_offset[i] is not None
        ]
        max_width = float(max(widths))
        if max_width == 0:
            raise AsciiTreePlot.NullEdgeLengthError(
                "Tree cannot be plotted under metric '%s' due to zero or null edge"
                " lengths: '%s'" % (self.plot_metric, tree._as_newick_string())
            )
        edge_scale_factor = float(effective_display_width) / max_width
        self.calc_plot(tree.seed_node, edge_scale_factor=edge_scale_factor)
        for i in range(len(tree.leaf_nodes()) * self.leaf_spacing_factor + 1):
            self.grid.append([" " for i in range(0, display_width)])
        self.draw_node(tree.seed_node)
        display = "\n".join(["".join(i) for i in self.grid])
        return display

    def calc_plot(self, node, edge_scale_factor):
        """
        First pass through tree, post-order traversal to calculate
        coordinates of each node.
        """
        child_nodes = node.child_nodes()
        if child_nodes:
            for n in child_nodes:
                self.calc_plot(n, edge_scale_factor)
            ys = [self.node_row[n] for n in child_nodes]
            self.node_row[node] = int(float((max(ys) - min(ys)) / 2) + min(ys))
        else:
            self.node_row[node] = self.current_leaf_row
            self.current_leaf_row = self.current_leaf_row + self.leaf_spacing_factor
        if node.edge.length is None:
            self.node_col[node] = 1
        else:
            self.node_col[node] = int(float(self.node_offset[node]) * edge_scale_factor)
        self.node_col[node] = int(float(self.node_offset[node]) * edge_scale_factor)

    def draw_label(self, label, row, start_col):
        if label:
            for i in range(len(label)):
                if start_col + i < len(self.grid[row]):
                    self.grid[row][start_col + i] = label[i]

    def draw_node(self, node):
        """
        Second pass through tree, plotting nodes onto given self.grid.
        """
        child_nodes = node.child_nodes()
        if child_nodes:
            for i, child_node in enumerate(child_nodes):
                start_row = min([self.node_row[node], self.node_row[child_node]])
                end_row = max([self.node_row[node], self.node_row[child_node]])
                if i == 0:
                    self.grid[self.node_row[child_node]][self.node_col[node]] = "/"
                    start_row = start_row + 1
                    edge_row = self.node_row[child_node]
                elif i == len(child_nodes) - 1:
                    self.grid[self.node_row[child_node]][self.node_col[node]] = "\\"
                    edge_row = self.node_row[child_node]
                else:
                    self.grid[self.node_row[child_node]][self.node_col[node]] = "+"
                    edge_row = self.node_row[child_node]
                self.draw_node(child_node)
                for x in range(self.node_col[node] + 1, self.node_col[child_node]):
                    self.grid[edge_row][x] = "-"
                for y in range(start_row, end_row):
                    self.grid[y][self.node_col[node]] = "|"
            label = []
            if self.show_internal_node_labels:
                label = self.get_label_for_node(node)
                self.draw_internal_text(label, self.node_row[node], self.node_col[node])
            else:
                self.grid[self.node_row[node]][self.node_col[node]] = "+"
        else:
            label = self.get_label_for_node(node)
            self.draw_label(label, self.node_row[node], self.node_col[node] + 1)

    def draw_internal_text(self, label, r, c):
        row = self.grid[r]
        try:
            for n, letter in enumerate(label):
                row[c + n] = letter
        except:
            pass

