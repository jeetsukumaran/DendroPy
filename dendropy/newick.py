#! /usr/bin/env python

############################################################################
##  newick.py
##
##  Part of the DendroPy phylogenetic computation library.
##
##  Copyright 2008 Jeet Sukumaran and Mark T. Holder.
##
##  This program is free software; you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation; either version 3 of the License, or
##  (at your option) any later version.
##
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License along
##  with this programm. If not, see <http://www.gnu.org/licenses/>.
##
############################################################################

"""
This module handles the reading and writing of trees in NEWICK format.
"""

import re
from dendropy import datasets
from dendropy import trees


def from_string(trees_string):
    newick_reader = NewickTreeReader()
    trees = newick_reader.read_trees(text=trees_string)
    if len(trees) == 0:
        return None
    elif len(trees) == 1:
        return trees[0]
    else:
        return trees

def to_string(tree):
    results = []
    newick_writer = NewickTreeWriter()
    if isinstance(tree, list) or isinstance(tree, trees.TreesBlock):
        for member_tree in tree:
            results.append(newick_writer.compose_tree(member_tree))
    else:
        results.append((newick_writer.compose_tree(tree)))
    return ';\n'.join(results)
            
class NewickTreeReader(datasets.Reader):
    """
    Implementation of TreeReader for NEWICK files and strings.
    """
    
    def __init__(self):
        """
        This has changed so much so many times that any documentation
        I put in here will probably be obselete in a matter of hours
        so this comment is all you are getting.
        """
        datasets.Reader.__init__(self)

    def read_trees(self, fileobj=None, text=None, trees_block=None):
        """
        Instantiates and returns a TreesBlock object based
        on the Newick-formatted contents read from the file
        descriptor object `fileobj`.
        """
        if trees_block is None:
            trees_block = self.trees_block_factory()
        if fileobj:
            return self.parse_trees(fileobj.read(), trees_block)
        else:
            return self.parse_trees(text, trees_block)
            
    def tree_iter(self, filepath=None, fileobj=None, text=None, trees_block=None):
        """
        Instantiates and returns a TreesBlock object based
        on the Newick-formatted contents read from the file
        descriptor object `fileobj`.
        """
        if trees_block is None:
            trees_block = self.trees_block_factory()
        filehandle = datasets.Reader.get_file_handle(filepath=filepath, fileobj=fileobj, text=text)
        statement_block = filehandle.read()           
        statement_block = statement_block.replace('\n','').replace('\r','')
        for statement in statement_block.split(';'):
            statement = statement.strip() + ';'
            newick_parser = NewickTreeParser()
            tree = newick_parser.parse_tree_statement(statement, trees_block)
            trees_block.pop()
            yield tree

    ## Following methods are class-specific ##

    def parse_trees(self, statement_block, trees_block, translate_dict=None):
        """
        Given a string block which defines trees in Newick format,
        this parses the Newick strings and adds the trees found to
        dataset.
        """
        statement_block = statement_block.replace('\n','').replace('\r','')
        tree_statements = []
        for statement in statement_block.split(';'):
            statement = statement.strip()
            if statement:
                tree_statements.append(statement + ';')
        newick_parser = NewickTreeParser()
        trees = []
        for tree_statement in tree_statements:
            newick_parser.parse_tree_statement(tree_statement, trees_block, translate_dict)
        return trees_block

class NewickTreeWriter(datasets.Writer):
    """
    Handles representation and serialization of a DendroPy Tree object
    in NEWICK format.
    """

    def __init__(self, **kwargs):
        """
        Instantiates the object, setting default for various
        formatting/representation options.
        """
        self.edge_lengths = True
        self.internal_labels = True
        self.support_as_edge_lengths = False
        self.support_as_labels = False
        self.support_as_percentages = False
        self.support_decimals = None

    def write_dataset(self, dataset, dest):
        """
        Writes a DataSet object to a full document-level
        representation of the format being implemented by the
        deriving class. 
        """
        for trees_block in dataset.trees_blocks:
            for tree in trees_block:
                dest.write(self.compose_node(tree.seed_node) + ';\n')                                

    ### Derived-class specific methods ###
    
    def compose_taxlabel(self, label):
        if re.search('[' + NewickTreeParser.whitespace + NewickTreeParser.punctuation + ']', label) != None:
            return "'" + label + "'"
        else:
            return label    
    
    def compose_tree(self, tree):
        """
        Convienience method.        
        """
        return self.compose_node(tree.seed_node)

    def choose_display_tag(self, node):
        """
        Based on current settings, the attributes of a node, and
        whether or not the node is a leaf, returns an appropriate tag.
        """
        if hasattr(node, 'taxon') and node.taxon:
            return self.compose_taxlabel(node.taxon.label)
        elif hasattr(node, 'label') and node.label:
            return self.compose_taxlabel(node.label)
        elif len(node.children()) == 0:
            # force label if a leaf node
            return self.compose_taxlabel(node.elem_id)
        else:
            return ""
        
    def compose_node(self, node):
        """
        Given a DendroPy Node, this returns the Node as a NEWICK
        statement according to the class-defined formatting rules.
        """
        children = node.children()
        if children:
            subnodes = [self.compose_node(child) for child in children]
            statement = '(' + ','.join(subnodes) + ')'
            if self.internal_labels:
                statement = statement + self.choose_display_tag(node)
            if node.edge.length != None and self.edge_lengths:
                try:
                    statement =  "%s:%f" \
                                % (statement, float(node.edge.length))
                except ValueError:
                    statement =  "%s:%s" \
                                % (statement, node.edge.length)
            return statement
        else:
            if self.internal_labels:
                statement = self.choose_display_tag(node)
            if node.edge.length != None and self.edge_lengths:
                try:
                    statement =  "%s:%0.10f" \
                                % (statement, float(node.edge.length))
                except ValueError:
                    statement =  "%s:%s" \
                                % (statement, node.edge.length)
            return statement

class NewickTreeParser(object):
    """
    Encapsulates process of generating a (single) DendroPy Tree object
    based on a (single) Newick tree statement string. Slow, but
    (fairly) robust.
    """

    punctuation = '\(\)\[\]\{\}\\\/\,\;\:\=\*\'\"\`\+\-\<\>'
    whitespace = ' \0\t\n\r'

    def __init__(self, trees_block_factory=None, tree_factory=None, node_factory=None, edge_factory=None):
        """
        Must be given tree factory to create trees.
        """
        self.statement = ''
        self.curr_pos = 0
        self.current_token = None
        if trees_block_factory is None:
            self.trees_block_factory = trees.TreesBlock
        else:
            self.trees_block_factory = trees_block_factory
        if tree_factory is None:
            self.tree_factory = trees.Tree
        else:
            self.tree_factory = tree_factory
        if node_factory is None:
            self.node_factory = trees.Node
        else:
            self.node_factory = node_factory
        if edge_factory is None:
            self.edge_factory = trees.TreesBlock
        else:
            self.edge_factory = edge_factory

    def parse_tree_statement(self, tree_statement, trees_block, translate_dict=None):
        """
        Processes a TREE command. Assumes that the input stream is
        located at the beginning of the statement (i.e., the first
        parenthesis that defines the tree).
        """      
        self.statement = tree_statement
        self.curr_pos = 0
        self.current_token = None
        child_nodes = []
        tree = self.tree_factory()
        token = self.read_next_token()
        while token and token != ';' and token != ':':
            # process nodes until no more tokens, end of tree
            # statement, or ':' is encountered, presumably outside
            # main tree parenthetical statement (i.e., length of root
            node = self.parse_tree_node(tree)
            if node:
                child_nodes.append(node)
            token = self.current_token
            if self.current_token == ')':
                # OK, I'll be the first to admit that this is rather
                # hacky but it works.
                # If an end-parenthesis is encountered ...
                token = self.read_next_token()
                if token and not token in NewickTreeParser.punctuation:
                    break
        for node in child_nodes:
            tree.seed_node.add_child(node)
        if token and not token in NewickTreeParser.punctuation:
            tree.seed_node.label = token
            token = self.read_next_token()
        if token and token == ':':
            length = self.read_next_token(ignore_punct='-')
            tree.seed_node.edge.length = length
            
        # convert labels at terminal nodes to taxa
        for node in tree.leaves():
            if node.label:
                if translate_dict and node.label in translate_dict:
                    label = translate_dict[node.label]
                else:
                    label = node.label                      
                node.taxon = trees_block.taxa_block.find_taxon(label=label, update=True)            
        trees_block.append(tree)
        return tree

    def parse_tree_node(self, tree):
        """
        Processes a TREE statement. Assumes that the file reader is
        positioned right after the '(' token in a TREE statement or
        right after a comma following a node inside a tree statement.
        """
        node = self.node_factory()
        token = self.read_next_token()
        while token and token != ')' and token != ',':            
            if token=="." or token not in NewickTreeParser.punctuation:
#                 if translate_dict and token in translate_dict:
#                     label = translate_dict[token]
#                 else:
#                     label = token
#                 node.taxon = taxa_block.find_taxon(label=label, update=True)
                if node.label is None:
                    node.label = token
                else:                    
                    node.label = node.label + token
            if token == ':':
                edge_length_str = self.read_next_token(ignore_punct='-')
                try:
                    node.edge.length = float(edge_length_str)
                except ValueError:
                    node.edge.length = edge_length_str
            if token == '(':
                while token and token != ')':
                    child_node = self.parse_tree_node(tree)
                    if child_node:
                        node.add_child(child_node)
                    token = self.current_token                  
            token = self.read_next_token()             
        return node

    def is_eos(self):
        """
        Returns True if currently at end of the string stream.
        """
        return self.curr_pos >= len(self.statement)

    def read_next_char(self):
        """
        Advances the stream cursor to the next character and returns
        it.
        """
        self.curr_pos = self.curr_pos + 1
        if self.curr_pos < len(self.statement):
            return self.statement[self.curr_pos]
        else:
            return ''

    def skip_to_significant_character(self):
        """
        Advances to the first non-whitespace character.
        """
        while (self.statement[self.curr_pos] in NewickTreeParser.whitespace) \
                  and not self.is_eos():
            self.read_next_char()

    def read_next_token(self, ignore_punct=None):
        """
        Reads the next token in the file stream. A token in this
        context is any word or punctuation character outside of a
        comment block.
        """
        if ignore_punct == None:
            ignore_punct = []
        if not self.is_eos():
            token = ''
            self.skip_to_significant_character()
            if not self.is_eos():
                if self.statement[self.curr_pos] == "'":
                    self.read_next_char()
                    end_quote = False
                    while not end_quote and not self.is_eos():
                        if self.statement[self.curr_pos] == "'":
                            self.read_next_char()
                            if self.statement[self.curr_pos] == "'":
                                token = token + "'"
                                self.read_next_char()
                            else:
                                end_quote = True
                        else:
                            token = token + self.statement[self.curr_pos]
                            self.read_next_char()
                else:
                    # it gets pretty hairy here ...

                    if (self.statement[self.curr_pos] \
                        in NewickTreeParser.punctuation) \
                           and (self.statement[self.curr_pos] \
                                not in ignore_punct):
                        token = self.statement[self.curr_pos]
                        self.read_next_char()
                    else:
                        while not self.is_eos() \
                                  and not ((self.statement[self.curr_pos] \
                                            in NewickTreeParser.whitespace) \
                              or (self.statement[self.curr_pos] in \
                                  NewickTreeParser.punctuation \
                                  and self.statement[self.curr_pos] \
                                  not in ignore_punct)):
                            token = token + self.statement[self.curr_pos]
                            self.read_next_char()
                self.current_token = token
            else:
                self.current_token = None
        else:
            self.current_token = None
        return self.current_token


if __name__ == "__main__":
    source = "tests/files/newick_trees.tre"
    sourcef = open(source, 'r')
    newick = NewickTreeReader()
    trees_block = newick.read_trees(sourcef)
    for tree_idx, tree in enumerate(trees_block):
        print
        print "%d.  Tree %s (%s)" % (tree_idx, tree.elem_id, tree.label)
        for node in tree.postorder_node_iter():
            print "%s (%s): Taxon=%s  Edge=%s" % (node.elem_id, node.label, node.taxon, node.edge.length)
    
    
