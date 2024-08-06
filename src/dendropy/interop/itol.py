#! /usr/bin/env python
# -*- coding: utf-8 -*-

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010-2015 Jeet Sukumaran and Mark T. Holder.
##  All rights reserved.
##
##  See "LICENSE.rst" for terms and conditions of usage.
##
##  If you use this work or any portion thereof in published work,
##  please cite it as:
##
##     Sukumaran, J. and M. T. Holder. 2010. DendroPy: a Python library
##     for phylogenetic computing. Bioinformatics 26: 1569-1571.
##
##############################################################################

"""
Wrappers for interacting with ITOL.
"""

import os
import tempfile
from urllib.parse import urlencode
from urllib.request import urlopen
from urllib import request
import zipfile
try:
    import zlib
    COMPRESSION_TYPE = zipfile.ZIP_DEFLATED
except:
    COMPRESSION_TYPE = zipfile.ZIP_STORED
from dendropy.utility import urlio
from dendropy.utility import error

class ItolService(object):

    DISPLAY_MODE_NORMAL = 1
    DISPLAY_MODE_CIRCULAR = 2
    DISPLAY_MODE_UNROOTED = 3

    RENDER_PARAM_NAMES = set([
        "newick_format",
        "display_mode",
        "tree_x",
        "tree_y",
        "vertical_shift_factor",
        "horizontal_scale_factor",
        "current_font_size",
        "current_font_name",
        "current_font_style",
        "leaf_sorting",
        "label_display",
        "align_labels",
        "label_shift",
        "dashed_lines",
        "ignore_branch_length",
        "slanted_normal",
        "arc",
        "rotation",
        "normal_rotation",
        "unrooted_rotation",
        "line_width",
        "default_branch_color",
        "default_label_color",
        "inverted",
        "circle_size_inverted",
        "range_mode",
        "include_ranges_legend",
        "ranges_legend_title",
        "internal_marks",
        "internal_scale",
        "internalScale1",
        "internalScale2",
        "internalScale1Color",
        "internalScale2Color",
        "branchlength_display",
        "branchlength_label_size",
        "branchlength_label_position",
        "branchlength_label_sci",
        "branchlength_label_rounding",
        "branchlength_label_age",
        "metadata_source",
        "bootstrap_display",
        "bootstrap_type",
        "bootstrap_symbol",
        "bootstrap_symbol_min",
        "bootstrap_symbol_max",
        "bootstrap_symbol_position",
        "bootstrap_symbol_color",
        "bootstrap_slider_min",
        "bootstrap_slider_max",
        "bootstrap_label_position",
        "bootstrap_label_size",
        "bootstrap_label_sci",
        "bootstrap_label_rounding",
        "bootstrap_width_min",
        "bootstrap_width_max",
        "bootstrap_min_color",
        "bootstrap_use_mid_color",
        "bootstrap_mid_color",
        "bootstrap_max_color",
        "datasets_visible",
        ])

    def __init__(self, **kwargs):
        """

        Keyword Arguments
        -----------------
        newick_format: int
            Possible value: ID. Only used when export format is set to newick. If specified, internal node IDs will be included in exported Newick trees.
        display_mode: int
            Possible values: 1,2 or 3 (1=normal, 2=circular, 3=unrooted)
        tree_x: int
            Tree will be shifted horizontally by this value, positive or negative (value in pixels)
        tree_y: int
            Tree will be shifted vertically by this value, positive or negative (value in pixels)
        vertical_shift_factor: int
            Vertical tree scaling factor (positive number)
        horizontal_scale_factor: int
            Horizontal tree scaling factor (positive number)
        current_font_size: int
            Main label font size in pixels (integer >= 9)
        current_font_name: int
            Possible values: 'Arial', 'Courier', 'Courier New', 'Verdana', 'Impact', 'Georgia', 'Times New Roman' or 'Monotype Corsiva'
        current_font_style: int
            Possible values: 'normal', 'bold' or 'bold italic'
        leaf_sorting: int
            Possible values: 1 or 2 (1=normal sorting, 2=no sorting)
        label_display: int
            Possible values: 0 or 1 (0=hide labels, 1=show labels)
        align_labels: int
            Possible values: 0 or 1 (0=labels not aligned, 1=labels aligned)
        label_shift: int
            Amount to shift the leaf labels (positive or negative, in pixels)
        dashed_lines: int
            Display dashed lines connecting branches to leaf labels. Possible values: 0 or 1 (0=hide lines, 1=display lines)
        ignore_branch_length: int
            Possible values: 0 or 1
        slanted_normal: int
            Display tree in slanted (triangular) mode. Note that display_mode
            must be 1. Possible values: 0 or 1
        arc: int
            Angle of the display arc for the circular tree (in degrees, a number between 0 and 360)
        rotation: int
            Rotation angle for circular tree (in degrees, a number between 0 and 360)
        normal_rotation: int
            Rotation angle for normal tree (in degrees, a number between 0 and 360)
        unrooted_rotation: int
            Rotation angle for unrooted tree (in degrees, a number between 0 and 360).
        line_width: int
            Width of tree lines (in pixels).
        default_branch_color: str
            Default color of tree's branches (hex, RGB or RGBA notation).
        default_label_color: str
            Default color of tree's labels (hex, RGB or RGBA notation).
        inverted: int
            Inverted display (only for circular and normal mode). Possible values: 0 or 1.
        circle_size_inverted: int
            For inverted circular display only. Internal radius will be
            increased by this amount (value in pixels)
        range_mode: int
            Colored ranges display style. Possible values: 0,1 or 2 (0=off,
            1=cover labels only, 2=cover full clades)<
        include_ranges_legend: int
            Include colored ranges legend. Possible values: 0 or 1.
        ranges_legend_title: str
            Title of the colored ranges legend.
        internal_marks: int
            Draw circles to mark the location of internal nodes. Possible
            values: 0, 1 or 2 (0=Do not display, 1=Display on nodes with one
            child, 2=Display always).
        internal_scale: int
            Draw internal tree scale. Possible values: 0 or 1.
        internalScale1: int
            Repeat value for the first set of internal scale lines (branch
            length value)
        internalScale2: int
            Repeat value for the second set of internal scale lines (branch
            length value).
        internalScale1Color: str
            Color for the first set of internal scale lines (hex, RGB or RGBA).
        internalScale2Color: str
            Color for the second set of internal scale lines (hex, RGB or RGBA).
        branchlength_display: int
            Display branch length text labels on branches (possible values: 0 or 1).
        branchlength_label_size: int
            Font size for the branch length labels (in pixels, integer >= 0)
        branchlength_label_position: int
            Position of the branch length labels along the branch (in percent,
            default is 50%).
        branchlength_label_sci: int
            Display branch length labels in scientific notation (if set to 1).
        branchlength_label_rounding
            Round branch length labels to this number of decimals (integer >= 0).
        branchlength_label_age
            Display node age instead of raw branch length values (if set to 1).
        metadata_source: str
            Which metadata source to use for bootstrap display options
            (default: 'bootstrap', other possible options depend on the data in
            the tree).
        bootstrap_display: int
            Display metadata values (possible values: 0 or 1).
        bootstrap_type: int
            Type of metadata display. Possible values: 1, 2, 3 or 4 (1=Symbol,
            2=Text label, 3=Branch color and 4=Branch width.
        bootstrap_symbol: int
            Symbol used to display metadata values. Possible values: 1, 2, 3 or 4 (1=Circle, 2=Triangle, 3=Square and 4=Star).
        bootstrap_symbol_min: int
            Minimum size for the metadata symbol (in pixels).
        bootstrap_symbol_max: int
            Maximum size for the metadata symbol (in pixels).
        bootstrap_symbol_position: int
            Position of the metadata symbol along the branch (in percent,
            default is 50%).
        bootstrap_symbol_color: str
            Bootstrap symbol color (hex, RGB or RGBA).
        bootstrap_slider_min: int
            Minimum metadata value to display.
        bootstrap_slider_max: int
            Maximum metadata value to display.
        bootstrap_label_position: int
            Position of the metadata text label along the branch (in percent,
            default is 50%).
        bootstrap_label_size: int
            Font size for the metadata text labels (in pixels, integer >= 9).
        bootstrap_label_sci: int
            Display metadata labels in scientific notation (if set to 1).
        bootstrap_label_rounding: int
            Round metadata labels to this number of decimals (integer >= 0).
        bootstrap_width_min: int
            Branch width for minimum metadata value (in pixels), when
            bootstrap_type=4.
        bootstrap_width_max: int
            Branch width for maximum metadata value (in pixels), when
            bootstrap_type=4.
        bootstrap_min_color: int
            Branch color for minimum metadata value (hex, RGB or RGBA), when
            bootstrap_type=3.
        bootstrap_use_mid_color: int
            Use a three color gradient for metadata branch colors (possible
                values: 0 or 1).
        bootstrap_mid_color: str
            Branch color for midpoint metadata value (hex, RGB or RGBA), when
            bootstrap_type=3 and bootstrap_use_mid_color=1.
        bootstrap_max_color: str
            Branch color for maximum metadata value (hex, RGB or RGBA), when
            bootstrap_type=3.
        datasets_visible: int
            Comma delimited list of datasets to display (starting with 0, e.g.
            datasets_visible=0,2,5).
        """
        self._base_upload_url = "https://itol.embl.de/batch_uploader.cgi"
        self._base_download_url = "https://itol.embl.de/batch_downloader.cgi"
        self._headers = {"Content-type": "application/x-www-form-urlencoded",
           "Accept": "text/plain"}
        self._render_params = {}
        for k, v in kwargs:
            setattr(self, k, v)

    def __setattr__(self, name, value):
        if name in ItolService.RENDER_PARAM_NAMES:
            if value is None:
                try:
                    del self._render_params[name]
                except KeyError:
                    pass
            else:
                self._render_params[name] = value
        object.__setattr__(self, name, value)

    def read_remote(self,
            tree_id,
            result_format):
        """
        Executes query and returns binary result.
        """
        params_d = {
                "tree": tree_id,
                "format": result_format,
        }
        params_d.update(self._render_params)
        params = urlencode(params_d)
        req = request.Request(
                self._base_download_url,
                params.encode("ascii"),
                self._headers)
        response = urlopen(req)
        result = response.read()
        return result

    def download(self,
            dest,
            tree_id,
            result_format):
        """
        Executes query and writes binary result to 'dest'.

        Examples
        --------

        ::

            from dendropy.interop import itol
            itol_service = itol.ItolService()
            itol_service.display_mode = 1
            itol_service.label_display = 0
            itol_service.ignore_branch_length = 1
            itol_service.normal_rotation = 90
            itol_service.slanted_normal = 1
            itol_service.line_width = 13
            itol_service.internal_scale = 0
            result_format = "pdf" # "png", "svg", etc.
            itol_service.download(
                    dest="x1.{}".format(result_format),
                    tree_id=709552230231221546743590,
                    result_format=result_format,)

        """
        result = self.read_remote(
                tree_id=tree_id,
                result_format=result_format,)
        if isinstance(dest, str):
            dest = open(dest, "wb")
        with dest:
            dest.write(result)

    def upload(self, tree):
        """
        Uploads a tree to iTOL.

        Examples
        --------

        ::

            import dendropy
            from dendropy.utility import itol
            tree = dendropy.Tree.get(data="[&R] (A1,(B2,(C3,(D4,E5))));",
                    schema="newick")
            itol_service = ItolService()
            tree_id = itol_service.upload(tree)
            print(tree_id)

        """
        import requests
        with tempfile.TemporaryDirectory() as temp_dirname:
            zf_path = os.path.join(temp_dirname, "data.zip")
            tree_str = tree.as_string(schema="newick")
            zf = zipfile.ZipFile(zf_path, mode="w", compression=COMPRESSION_TYPE)
            # zf = zipfile.ZipFile(zf_path, mode="w")
            with zf:
                zf.writestr("tree.tree", tree_str,)
            with open(zf_path, "rb") as zf:
                files = {"zipFile": open(zf_path, "rb")}
            response = urlio.post_request(
                    url=self._base_upload_url,
                    files=files)
            response_text = response.text
            if "SUCCESS" not in response_text:
                raise error.ExternalServiceError(
                        service_name="iTOL",
                        invocation_command=self._base_upload_url,
                        service_input="",
                        returncode=-1,
                        stdout=response_text,
                        stderr="")
            response_lines = [row for row in response_text.split("\n") if row]
            warnings = response_lines[0:-1]
            tree_id = response_lines[-1].split()[1]
            return tree_id

    def render_tree(self,
            tree,
            dest,
            result_format,):
        """
        Renders a tree using the iTOL web service.

        Examples
        --------

        ::

            itol_service = ItolService()
            result_format = "pdf"
            tree = dendropy.Tree.get(data="[&R] (A1,(B2,(C3,(D4,E5))));",
                    schema="newick")
            itol_service.render_tree(
                    tree=tree,
                    dest="sample.{}".format(result_format),
                    result_format=result_format,)
        """
        tree_id = self.upload(tree)
        self.download(
                dest=dest,
                tree_id=tree_id,
                result_format=result_format,)
        return tree_id

