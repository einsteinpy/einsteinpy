import os
from datetime import datetime

import alabaster

import einsteinpy

project = "EinsteinPy"
year = datetime.now().year
copyright = "%d EinsteinPy Development Team" % year

version = einsteinpy.__version__

release = version
highlight_language = "python"
pygments_style = "sphinx"
autoclass_content = "both"

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "astropy": ("http://docs.astropy.org/en/stable/", None),
    "numpy": ("https://docs.scipy.org/doc/numpy/", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/reference", None),
    "matplotlib": ("https://matplotlib.org", None),
    "sympy": ("https://docs.sympy.org/latest", None),
}

autodoc_member_order = "bysource"

html_theme_options = {
    "logo": "logo_small.png",
    "logo_name": True,
    "logo_text_align": "center",
    "travis_button": True,
    "codecov_button": True,
    "description": "General Relativity in Python",
    "body_text_align": "left",
    "github_user": "einsteinpy",
    "github_repo": "einsteinpy",
    "show_relbars": True,
    "show_powered_by": False,
    "page_width": "80%",
    "github_banner": True,
    "extra_nav_links": {"Blog": "https://docs.einsteinpy.org/"},
}

add_function_parentheses = True

add_module_names = True

needs_sphinx = "1.3"
extensions = [
    "alabaster",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.todo",
    "sphinx.ext.mathjax",
    "sphinx.ext.intersphinx",
    "nbsphinx",
    "IPython.sphinxext.ipython_console_highlighting",
    "sphinx.ext.mathjax",  # New module for matrix visualization
    "sphinx.ext.graphviz",  # For creating the diagrams
]

#Nbsphinx configuration
# See https://github.com/jupyter/nbconvert/issues/878#issuecomment-419655951
# Should not be needed after nbconvert 5.5 is out
nbsphinx_kernel_name = "python3"

def setup(app):
    # https://docs.readthedocs.io/en/latest/guides/adding-custom-css.html
    # https://www.sphinx-doc.org/en/master/extdev/appapi.html#sphinx.application.Sphinx.add_js_file
    app.add_js_file('https://cdnjs.cloudflare.com/ajax/libs/require.js/2.3.6/require.min.js')
    app.add_js_file("https://unpkg.com/@jupyter-widgets/html-manager@^0.14.0/dist/embed-amd.js")

if os.environ.get('READTHEDOCS') == 'True':
    nbsphinx_execute = 'never'
else:
    nbsphinx_execute = 'always'

    # Controls when a cell will time out (defaults to 30; use -1 for no timeout):
    nbsphinx_timeout = 60

templates_path = ["_templates"]

source_suffix = ".rst"

master_doc = "index"

html_theme = "alabaster"

html_theme_path = [alabaster.get_path()]

html_title = "EinsteinPy"

html_static_path = ["_static"]

htmlhelp_basename = "einsteinpydoc"

html_sidebars = {
    "**": [
        "about.html",
        "navigation.html",
        "relations.html",
        "searchbox.html",
        "donate.html",
    ]
}
