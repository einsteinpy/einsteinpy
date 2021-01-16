from datetime import datetime

import alabaster

project = "EinsteinPy"
year = datetime.now().year
copyright = "%d EinsteinPy Development Team" % year

version = "0.3"
release = "0.3.1"
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


def setup(app):
    app.add_js_file('https://cdnjs.cloudflare.com/ajax/libs/require.js/2.3.6/require.min.js')

autodoc_member_order = "bysource"

html_theme_options = {
    "logo": "logo_small.png",
    "logo_name": True,
    "logo_text_align": "center",
    "codecov_button": True,
    "description": "General Relativity in Python",
    "body_text_align": "left",
    "github_user": "einsteinpy",
    "github_repo": "einsteinpy",
    "show_relbars": True,
    "show_powered_by": False,
    "page_width": "80%",
    "github_banner": True,
    "donate_url": "https://opencollective.com/einsteinpy",
    "github_button": True,
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
    "sphinx.ext.viewcode",  # View Source button
]
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
