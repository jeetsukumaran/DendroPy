#! /usr/bin/env python
# -*- coding: utf-8 -*-

from dendropy.utility import terminal

class TreeDrawComposer(object):
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
        leaf_spacing_factor : int
            Positive integer: number of rows between each leaf.
        width : int
            Force a particular display width, in terms of number of columns.
        """
        self.plot_metric = kwargs.pop("plot_metric", "depth")
        self.leaf_spacing_factor = kwargs.pop("leaf_spacing_factor", 2)
        self.width = kwargs.pop("width", None)
        self.display_width = kwargs.pop("display_width", self.width)  # legacy
        if kwargs:
            raise TypeError("Unrecognized or unsupported arguments: {}".format(kwargs))

    def reset(self):
        self.node_row = {}
        self.node_col = {}
        self.node_offset = {}
        self.current_leaf_row = 0
        self.line_segments = []
        self.internal_nodes = []
        self.internal_node_positions = []
        self.leaf_nodes = []
        self.leaf_node_positions = []

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

    def compose(self, tree):
        self.reset()
        if self.display_width is None:
            display_width = terminal.terminal_width() - 1
        else:
            display_width = self.display_width
        
        # Calculate effective display width (simplified from original)
        effective_display_width = display_width - 10  # Leave some margin
        
        self._calc_node_offsets(tree)
        widths = [
            self.node_offset[i]
            for i in tree.leaf_node_iter()
            if self.node_offset[i] is not None
        ]
        max_width = float(max(widths))
        if max_width == 0:
            raise TreeDrawComposer.NullEdgeLengthError(
                "Tree cannot be plotted under metric '%s' due to zero or null edge"
                " lengths: '%s'" % (self.plot_metric, tree._as_newick_string())
            )
        edge_scale_factor = float(effective_display_width) / max_width
        self.calc_plot(tree.seed_node, edge_scale_factor=edge_scale_factor)
        
        # Collect coordinates instead of drawing
        self.collect_coordinates(tree.seed_node)
        
        # Calculate bounds
        all_x_coords = [pos[0] for pos in self.internal_node_positions + self.leaf_node_positions]
        all_y_coords = [pos[1] for pos in self.internal_node_positions + self.leaf_node_positions]
        
        # Also consider line segment coordinates for bounds
        for segment in self.line_segments:
            all_x_coords.extend([segment[0][0], segment[1][0]])
            all_y_coords.extend([segment[0][1], segment[1][1]])
        
        bounds_x_max = max(all_x_coords) if all_x_coords else 0
        bounds_y_max = max(all_y_coords) if all_y_coords else 0
        
        return {
            'line_segments': self.line_segments,
            'internal_nodes': self.internal_nodes,
            'internal_node_positions': self.internal_node_positions,
            'leaf_nodes': self.leaf_nodes,
            'leaf_node_positions': self.leaf_node_positions,
            'bounds_x_max': bounds_x_max,
            'bounds_y_max': bounds_y_max
        }

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
        
        self.node_col[node] = int(float(self.node_offset[node]) * edge_scale_factor)

    def collect_coordinates(self, node):
        """
        Second pass through tree, collecting coordinates instead of drawing.
        """
        child_nodes = node.child_nodes()
        
        if child_nodes:
            # This is an internal node
            self.internal_nodes.append(node)
            self.internal_node_positions.append((self.node_col[node], self.node_row[node]))
            
            for i, child_node in enumerate(child_nodes):
                # Recursively process child
                self.collect_coordinates(child_node)
                
                # Add horizontal line from parent to child
                self.line_segments.append([
                    (self.node_col[node], self.node_row[child_node]),
                    (self.node_col[child_node], self.node_row[child_node])
                ])
            
            # Add vertical line connecting all children (if more than one child)
            if len(child_nodes) > 1:
                child_rows = [self.node_row[child] for child in child_nodes]
                min_row = min(child_rows)
                max_row = max(child_rows)
                self.line_segments.append([
                    (self.node_col[node], min_row),
                    (self.node_col[node], max_row)
                ])
        else:
            # This is a leaf node
            self.leaf_nodes.append(node)
            self.leaf_node_positions.append((self.node_col[node], self.node_row[node]))