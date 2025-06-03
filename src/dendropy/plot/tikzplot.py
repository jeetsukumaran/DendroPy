#! /usr/bin/env python
# -*- coding: utf-8 -*-

from dendropy.utility import error

class TikzTreePlot(object):

    @classmethod
    def _display(cls, tikz_code):
        try:
            from jupyter_tikz import TexFragment
        except ImportError:
            raise error.LibraryDependencyError("This method requires `jupyter_tikz` to be installed in the Python environment: `$ python3 -m pip install jupyter_tikz`")
        tikz_picture = TexFragment(tikz_code)
        return tikz_picture.run_latex()

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
        show_external_node_labels : bool
            Whether or not to write out external node labels.
        show_edge_labels : bool
            Whether or not to write out edge labels.
        node_label_compose_fn : function object
            A function that takes a Node object as an argument and returns
            the string to be used to display it.
        edge_label_compose_fn : function object
            A function that takes a Node object as an argument and returns
            the string to be used to display its incoming edge label.
        scale : float
            Scale factor for the entire tree plot.
        """
        self.plot_metric = kwargs.pop("plot_metric", "depth")
        self.show_external_node_labels = kwargs.pop("show_external_node_labels", True)
        self.show_internal_node_labels = kwargs.pop("show_internal_node_labels", False)
        self.show_edge_labels = kwargs.pop("show_edge_labels", False)
        self.compose_node = kwargs.pop("node_label_compose_fn", None)
        self.compose_edge = kwargs.pop("edge_label_compose_fn", None)
        self.scale = kwargs.pop("scale", 1.0)
        if self.compose_node is None:
            self.compose_node = self.default_compose_node
        if self.compose_edge is None:
            self.compose_edge = self.default_compose_edge

        self.tikzpicture_options = kwargs.pop("tikzpicture_options", None)
        self.scope_options = kwargs.pop("scope_options", None)
        self.extend_drawing_area_x = kwargs.pop("extend_drawing_area_x", None)
        self.extend_drawing_area_y = kwargs.pop("extend_drawing_area_y", None)

        if kwargs:
            raise TypeError("Unrecognized or unsupported arguments: {}".format(kwargs))

    def default_compose_node(self, node):
        if node.taxon is not None and node.taxon.label is not None:
            return node.taxon.label
        elif node.label is not None:
            return node.label
        else:
            return "@"

    def default_compose_edge(self, node):
        """Default edge label composer - returns edge length if available."""
        if node.edge is not None and node.edge.length is not None:
            return str(node.edge.length)
        else:
            return ""

    def reset(self):
        self.node_coords = {}
        self.node_label_map = {}
        self.edge_label_map = {}

    def _calc_node_offsets(self, tree):
        if self.plot_metric == "age" or self.plot_metric == "depth":
            for nd in tree.postorder_node_iter():
                cnds = nd.child_nodes()
                if self.plot_metric == "depth":  # 'number of branchings from tip'
                    if len(cnds) == 0:
                        curr_node_offset = 0.0
                    else:
                        depths = [self.node_coords[v][0] for v in cnds]
                        curr_node_offset = max(depths) + 1
                elif self.plot_metric == "age":  # 'sum of edge weights from tip'
                    if len(cnds) == 0:
                        curr_node_offset = 0.0
                    else:
                        if cnds[0].edge.length is not None:
                            curr_node_offset = self.node_coords[cnds[0]][0] + cnds[0].edge.length
                else:
                    raise ValueError(
                        "Unrecognized plot metric '%s' (must be one of: 'age', 'depth',"
                        " 'level', or 'length')"
                        % self.plot_metric
                    )
                self.node_coords[nd] = [curr_node_offset, 0]  # [x, y] coordinates
            flipped_origin = max(coord[0] for coord in self.node_coords.values())
            for nd in self.node_coords:
                self.node_coords[nd][0] = flipped_origin - self.node_coords[nd][0]
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
                    self.node_coords[nd] = [curr_edge_len, 0]
                else:
                    self.node_coords[nd] = [
                        curr_edge_len + self.node_coords[nd._parent_node][0],
                        0
                    ]

    def _assign_y_coordinates(self, tree):
        """Assign y-coordinates to nodes based on their position in the tree."""
        leaf_nodes = list(tree.leaf_node_iter())
        y_spacing = 1.0
        for i, node in enumerate(leaf_nodes):
            self.node_coords[node][1] = i * y_spacing

        def assign_internal_y_coords(node):
            if node.is_leaf():
                return
            children = node.child_nodes()
            for child in children:
                assign_internal_y_coords(child)
            y_coords = [self.node_coords[child][1] for child in children]
            self.node_coords[node][1] = sum(y_coords) / len(y_coords)

        assign_internal_y_coords(tree.seed_node)

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

    def get_label_for_edge(self, node):
        """Get label for the incoming edge of a node."""
        try:
            return self.edge_label_map[node]
        except KeyError:
            if self.show_edge_labels and node._parent_node is not None:
                label = self.compose_edge(node)
            else:
                label = ""
            self.edge_label_map[node] = label
            return label

    def display(self, tree):
        return self.__class__._display(self.compose(tree))

    def compose(self, tree):
        """Generate TikZ code for the tree."""
        self.reset()
        self._calc_node_offsets(tree)
        self._assign_y_coordinates(tree)

        # Scale coordinates
        for node in self.node_coords:
            self.node_coords[node] = [
                coord * self.scale for coord in self.node_coords[node]
            ]

        # Generate TikZ code
        tikz_code = []
        tikzpicture_options = f"[{self.tikzpicture_options}]" if self.tikzpicture_options else ""
        tikz_code.append("\\begin{tikzpicture}" + tikzpicture_options)

        # Draw edges
        for node in tree.preorder_node_iter():
            if node._parent_node is not None:
                parent_coords = self.node_coords[node._parent_node]
                node_coords = self.node_coords[node]
                tikz_code.append(
                    f"\\draw ({parent_coords[0]},{parent_coords[1]}) -- "
                    f"({node_coords[0]},{node_coords[1]});"
                )

        # Draw edge labels
        if self.show_edge_labels:
            for node in tree.preorder_node_iter():
                if node._parent_node is not None:
                    edge_label = self.get_label_for_edge(node)
                    if edge_label:
                        parent_coords = self.node_coords[node._parent_node]
                        node_coords = self.node_coords[node]
                        # Calculate midpoint of edge
                        mid_x = (parent_coords[0] + node_coords[0]) / 2
                        mid_y = (parent_coords[1] + node_coords[1]) / 2
                        tikz_code.append(
                            f"\\node[font=\\tiny,fill=white,inner sep=1pt] at "
                            f"({mid_x},{mid_y}) {{{edge_label}}};"
                        )

        # Draw nodes and labels
        for node in tree.preorder_node_iter():
            coords = self.node_coords[node]
            label = self.get_label_for_node(node)

            # Draw node
            if node.is_leaf():
                tikz_code.append(
                    f"\\node[anchor=west] at ({coords[0]},{coords[1]}) {{{label}}};"
                )
            else:
                tikz_code.append(
                    f"\\node[circle,fill,inner sep=1pt] at ({coords[0]},{coords[1]}) {{}};"
                )
                if label:
                    tikz_code.append(
                        f"\\node[anchor=east] at ({coords[0]-0.1},{coords[1]}) {{{label}}};"
                    )

        if self.extend_drawing_area_x or self.extend_drawing_area_y:
            xs = [coords[0] for coords in self.node_coords.values()]
            ys = [coords[1] for coords in self.node_coords.values()]
            min_x = min(xs)
            max_x = max(xs)
            min_y = min(ys)
            max_y = max(ys)
            extend_x = self.extend_drawing_area_x if self.extend_drawing_area_x is not None else 0
            extend_y = self.extend_drawing_area_y if self.extend_drawing_area_y is not None else 0
            tikz_code.append(
                f"\\path[draw=none] ({min_x-extend_x},{min_y-extend_y}) rectangle ({max_x+extend_x},{max_y+extend_y});"
            )

        tikz_code.append("\\end{tikzpicture}")
        return "\n".join(tikz_code)

    def draw(self, tree, dest):
        """Write TikZ code to the destination stream."""
        dest.write(self.compose(tree))
