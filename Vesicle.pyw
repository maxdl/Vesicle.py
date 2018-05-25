#!/usr/bin/env python3

import wx
from vesicle import frame


def main():
    app = wx.App()
    mainframe = frame.Frame(None)
    mainframe.Show(True)
    app.MainLoop()

if __name__ == '__main__':
    main()