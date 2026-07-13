"""Minimal tkinter compatibility shim for CPA on Python builds without _tkinter.

CPA imports only `N` from tkinter in its model module. The cluster Python used
for the CPA environment has no `_tkinter` extension, so this shim supplies that
constant without pulling in the GUI toolkit.
"""

N = "n"

