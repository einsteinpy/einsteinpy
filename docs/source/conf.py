import alabaster
from datetime import datetime
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
    "matplotlib": ("http://matplotlib.org", None),
}

autodoc_member_order = 'bysource'

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
    "extra_nav_links": {"Blog": "https://einsteinpy.github.io/"},
}

add_function_parentheses = True

add_module_names = True

needs_sphinx = "1.3"
extensions = [
    'alabaster',
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.todo',
    'sphinx.ext.mathjax',
    'sphinx.ext.intersphinx',
    'nbsphinx',
    'IPython.sphinxext.ipython_console_highlighting',
    'sphinx.ext.mathjax', #New module for matrix visualization
    'sphinx.ext.graphviz', # For creating the diagrams
]

def setup(app):
    # https://docs.readthedocs.io/en/latest/guides/adding-custom-css.html
    # https://www.sphinx-doc.org/en/master/extdev/appapi.html#sphinx.application.Sphinx.add_js_file
    app.add_js_file('https://cdnjs.cloudflare.com/ajax/libs/require.js/2.3.6/require.min.js')
    app.add_js_file("https://unpkg.com/@jupyter-widgets/html-manager@*/dist/embed-amd.js")

templates_path = ["_templates"]

source_suffix = ".rst"

master_doc = "index"

html_theme = "alabaster"

html_theme_path = [alabaster.get_path()]

html_title = "EinsteinPy"

htmlhelp_basename = 'einsteinpydoc'
html_sidebars = {
    "**": [
        "about.html",
        "navigation.html",
        "relations.html",
        "searchbox.html",
        "donate.html",
    ]
}
