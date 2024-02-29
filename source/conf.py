# -*- coding: utf-8 -*-
#
# Configuration file for the Sphinx documentation builder.
#
# This file does only contain a selection of the most common options. For a
# full list see the documentation:
# http://www.sphinx-doc.org/en/master/config

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#

import os
import sys
from importlib.metadata import version as get_version
# sys.path.insert(0, os.path.abspath('./seispy/seispy'))
# sys.path.insert(0, os.path.abspath('/Users/xumj/Codes/seispy/seispy'))


# -- Project information -----------------------------------------------------

project = 'Seispy'
copyright = '2020, Mijian Xu'
author = 'Mijian Xu'

# The short X.Y version
version = f"v{get_version('python-seispy')}"
# The full version, including alpha/beta/rc tags
release = version


# -- General configuration ---------------------------------------------------

# If your documentation needs a minimal Sphinx version, state it here.
#
# needs_sphinx = '1.0'

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.todo',
    'sphinx.ext.coverage',
    'sphinx.ext.viewcode',
    'sphinx.ext.mathjax',
    "myst_nb",
    # "myst_parser",
    'sphinx.ext.githubpages',
    "sphinx.ext.intersphinx",
    "sphinx_cjkspace.cjkspace",
    "sphinx_copybutton",
    "sphinx_design",
    # "sphinx_panels",
    "numpydoc",
]

panels_add_bootstrap_css = False

myst_enable_extensions = [
  "colon_fence",
  "dollarmath",
  "amsmath",
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
#
source_suffix = ['.rst','.md']
# source_suffix = '.rst'

# The master toctree document.
master_doc = 'index'

# The language for content autogenerated by Sphinx. Refer to documentation
# for a list of supported languages.
#
# This is also used if you do content translation via gettext catalogs.
# Usually you set "language" from the command line for these cases.
language = 'en'

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path .
exclude_patterns = []

# The name of the Pygments (syntax highlighting) style to use.
# pygments_style = 'sphinx'

myst_heading_anchors = 3

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#

html_baseurl = "seispy.xumijian.me"
html_theme = 'furo'
# html_theme = 'pydata_sphinx_theme'
# html_theme_path = [lsst_dd_rtd_theme.get_html_theme_path()]

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
#
# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# Custom sidebar templates, must be a dictionary that maps document names
# to template names.
templates_path = ['_templates']

#
# The default sidebars (for documents that don't match any pattern) are
# defined by theme itself.  Builtin themes are using these templates by
# default: ``['localtoc.html', 'relations.html', 'sourcelink.html',
# 'searchbox.html']``.
html_theme_options = {
#     "icon_links": [
#         {
#             # Label for this link
#             "name": "GitHub",
#             # URL where the link will redirect
#             "url": "https://github.com/xumi1993/seispy",  # required
#             # Icon class (if "type": "fontawesome"), or path to local image (if "type": "local")
#             "icon": "fab fa-github-square",
#             # Whether icon should be a FontAwesome class, or a local file
#             # "type": "fontawesome OR local",  # Default is fontawesome
#         }
#    ],
    "light_css_variables": {
        "color-brand-primary": "#0080C4",
        "color-brand-content": "#0080C4",
    },
}

#
# html_sidebars = {}
html_favicon = os.path.abspath(os.path.join('.', '_static', 'seispy_100.png'))
html_logo = os.path.abspath(os.path.join('.', '_static', 'seispy-long.png'))

# -- Options for HTMLHelp output ---------------------------------------------

# Output file base name for HTML help builder.
htmlhelp_basename = 'seispydoc'


# -- Options for LaTeX output ------------------------------------------------

latex_elements = {
    # The paper size ('letterpaper' or 'a4paper').
    #
    'papersize': 'a4paper',

    # The font size ('10pt', '11pt' or '12pt').
    #
    'pointsize': '10pt',

    # Additional stuff for the LaTeX preamble.
    #
    'preamble': '',

    # Latex figure (float) alignment
    #
    'figure_align': 'htbp',
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_documents = [
    (master_doc, 'seispy.tex', 'seispy Documentation',
     'Mijian Xu', 'manual'),
]


# -- Options for manual page output ------------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [
    (master_doc, 'seispy', 'seispy Documentation',
     [author], 1)
]


# -- Options for Texinfo output ----------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
    (master_doc, 'seispy', 'seispy Documentation',
     author, 'seispy', 'One line description of project.',
     'Miscellaneous'),
]


# -- Extension configuration -------------------------------------------------

# -- Options for todo extension ----------------------------------------------

# If true, `todo` and `todoList` produce output, else they produce nothing.
todo_include_todos = True