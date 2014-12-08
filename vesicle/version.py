# -*- coding: utf-8 -*-

import os.path
import sys

version = "1.0pre1"
date = ("November", "11", "2014")
title = "Vesicle"
author = "Max Larsson"
email = "max.larsson@liu.se"
homepage = "www.hu.liu.se/forskning/larsson-max/software"
if hasattr(sys, 'frozen'):
    if '_MEIPASS2' in os.environ:
        path = os.environ['_MEIPASS2']
    else:
        path = sys.argv[0]
else:
    path = __file__
app_path = os.path.dirname(path)
icon = os.path.join(app_path, "ves.ico")
