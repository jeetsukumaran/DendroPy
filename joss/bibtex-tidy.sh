#!/bin/bash

set -e # exit with error if any of this fails

script_dir="$(dirname "$(readlink -f "$0")")"

cd "${script_dir}"

bibtex-tidy --omit=abstract,keywords --curly --numeric --space=2 --align=0 --sort=key --duplicates=key,doi --merge=combine --strip-enclosing-braces --drop-all-caps --no-escape --sort-fields=title,shorttitle,author,year,month,day,journal,booktitle,location,on,publisher,address,series,volume,number,pages,doi,isbn,issn,url,urldate,copyright,category,note,metadata --trailing-commas --remove-empty-fields --no-backup paper.bib
