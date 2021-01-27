#!/usr/bin/env python3

import sys
import os

# append the parent directory to path (to get algs module in path)
sys.path.append(
    os.path.abspath(
        os.path.join(
            os.path.curdir, os.path.pardir
        )
    )
)

from align import algs

algs.main()
