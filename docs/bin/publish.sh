#! /bin/bash

set -e -o pipefail

netlify deploy --prod --site "b1510f2f-5f7e-400a-82e5-0c78ffedcc8b" -d build/html
