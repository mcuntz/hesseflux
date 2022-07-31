#!/usr/bin/env python
"""
fgui : GUI dialogs to choose files and directories using Tkinter.

This module was written by Matthias Cuntz while at Department of
Computational Hydrosystems, Helmholtz Centre for Environmental
Research - UFZ, Leipzig, Germany, and continued while at Institut
National de Recherche pour l'Agriculture, l'Alimentation et
l'Environnement (INRAE), Nancy, France.

Copyright (c) 2015-2020 Matthias Cuntz - mc (at) macu (dot) de
Released under the MIT License; see LICENSE file for details.

* Written Jun 2014 by Matthias Cuntz (mc (at) macu (dot) de)
* Added directories_from_gui, Oct 2015, Matthias Cuntz
* Using numpy docstring format, May 2020, Matthias Cuntz

.. moduleauthor:: Matthias Cuntz

The following functions are provided

.. autosummary::
   directory_from_gui
   directories_from_gui
   file_from_gui
   files_from_gui
"""
from __future__ import division, absolute_import, print_function


__all__ = ['directory_from_gui', 'directories_from_gui', 'file_from_gui', 'files_from_gui']


# -------------------------------------------------------------------------
# Choose directories in GUI
#

def directories_from_gui(initialdir='.', title='Choose one or several directories.'):
    """
    Open dialog to select several directories.

    Parameters
    ----------
    initialdir : str, optional
        Initial directory, in which opens GUI (default: '.')
    title : str, optional
        Title of GUI (default: 'Choose one or several directories.'

    Returns
    ------
    list
        Selected directories.

    Examples
    --------
    .. code-block:: python

       if not dirs:
           dirs = directories_from_gui()
           if not dirs:
               raise ValueError('Error: no directories given.')

    History
    -------
    Written,  Matthias Cuntz, Oct 2015
    Modified, Matthias Cuntz, May 2020 - numpy docstring format
    """
    try:                # Python 3
        import tkinter as Tkinter
        import tkinter.filedialog as tkFileDialog
    except ImportError: # Python 2
        import Tkinter, tkFileDialog
    root = Tkinter.Tk()
    root.withdraw()                                      # hide root window, i.e. white square

    # always on top
    root.tk.call('wm', 'attributes', '.', '-topmost', 1) # focus on (hidden) window so that child is on top

    # # Make it almost invisible - no decorations, 0 size, top left corner.
    # # Then show window again and lift it to top so it can get focus
    # root.overrideredirect(True)
    # root.geometry('0x0+0+0')
    # root.deiconify()
    # root.lift()
    # root.focus_force()

    idir = initialdir
    alldirs = []
    while True:
        dirs = tkFileDialog.askdirectory(parent=root,
                                         title=title,
                                         initialdir=idir)
        if not dirs: break
        alldirs.append(dirs)
        idir = dirs

    root.destroy()

    return alldirs


# -------------------------------------------------------------------------
# Choose one directory in GUI
#

def directory_from_gui(initialdir='.', title='Choose directory.'):
    """
    Opens dialog to select directory.

    Parameters
    ----------
    initialdir : str, optional
        Initial directory, in which opens GUI (default: '.')
    title : str, optional
        Title of GUI (default: 'Choose directory.')

    Returns
    ------
    str
        Selected directory.

    Examples
    --------
    .. code-block:: python

       if not idir:
           idir = directory_from_gui()
           if not idir:
               raise ValueError('Error: no directory given.')

    History
    -------
    Written,  Matthias Cuntz, Jun 2014
    Modified, Matthias Cuntz, May 2020 - numpy docstring format
    """
    try:                # Python 3
        import tkinter as Tkinter
        import tkinter.filedialog as tkFileDialog
    except ImportError: # Python 2
        import Tkinter, tkFileDialog
    root = Tkinter.Tk()
    root.withdraw()                                      # hide root window, i.e. white square

    # always on top
    root.tk.call('wm', 'attributes', '.', '-topmost', 1) # focus on (hidden) window so that child is on top

    # # Make it almost invisible - no decorations, 0 size, top left corner.
    # # Then show window again and lift it to top so it can get focus
    # root.overrideredirect(True)
    # root.geometry('0x0+0+0')
    # root.deiconify()
    # root.lift()
    # root.focus_force()

    dirs = tkFileDialog.askdirectory(parent=root,
                                     title=title,
                                     initialdir=initialdir)

    root.destroy()

    return dirs



# -------------------------------------------------------------------------
# Choose one file in GUI
#

def file_from_gui(initialdir='.', title='Choose file', multiple=False):
    """
    Wrapper for :func:`files_from_gui` with multiple=False, i.e.
    open dialog to select one file.

    Examples
    --------
    .. code-block:: python

       if not file:
           file = file_from_gui()
           if not file:
               raise ValueError('Error: no input file given.')

    History
    -------
    Written,  Matthias Cuntz, Jun 2014
    Modified, Matthias Cuntz, May 2020 - numpy docstring format
    """
    return files_from_gui(initialdir=initialdir, title=title, multiple=multiple)


# -------------------------------------------------------------------------
# Choose files in GUI
#

def files_from_gui(initialdir='.', title='Choose file(s).', multiple=True):
    """
    Open dialog to select one or several files.

    Parameters
    ----------
    initialdir : str, optional
        Initial directory, in which opens GUI (default: '.')
    title : str, optional
        Title of GUI (default: 'Choose file(s).')
    multiple : bool, optional
        True:  allow selection of multiple files (default).
        False: only one single file possible to select.

    Returns
    ------
    list
        Selected files.

    Note
    ----
    It always returns a list even with `multiple=False`.

    Examples
    --------
    .. code-block:: python

       if not files:
           files = files_from_gui()
           if not files:
               raise ValueError('Error: no input file(s) given.')

    History
    -------
    Written,  Matthias Cuntz, Jun 2014
    Modified, Matthias Cuntz, May 2020 - numpy docstring format
    """
    try:                # Python 3
        import tkinter as Tkinter
        import tkinter.filedialog as tkFileDialog
    except ImportError: # Python 2
        import Tkinter, tkFileDialog
    root = Tkinter.Tk()
    root.withdraw()                                      # hide root window, i.e. white square

    # always on top
    root.tk.call('wm', 'attributes', '.', '-topmost', 1) # focus on (hidden) window so that child is on top

    # # Make it almost invisible - no decorations, 0 size, top left corner.
    # # Then show window again and lift it to top so it can get focus
    # root.overrideredirect(True)
    # root.geometry('0x0+0+0')
    # root.deiconify()
    # root.lift()
    # root.focus_force()

    files = tkFileDialog.askopenfilename(parent=root,
                                         title=title,
                                         multiple=multiple,
                                         initialdir=initialdir)
    files = list(root.tk.splitlist(files))
    
    root.destroy()

    return files


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
