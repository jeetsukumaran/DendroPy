#!/usr/bin/env bash

python3 -m uv pip compile requirements.in --python 3.8 > requirements.txt
