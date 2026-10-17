from datetime import datetime

import alabaster

project = "EinsteinPy"
year = datetime.now().year
copyright = "%d EinsteinPy Development Team" % year

version = "0.4"
release = "0.4.dev0"
highlight_language = "python"
pygments_style = "sphinx"
autoclass_content = "both"

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "astropy": ("http://docs.astropy.org/en/stable/", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
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
    "badge_branch": "main",
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
    "extra_nav_links": {"Blog": "https://blog.einsteinpy.org/"},
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

# pdflatex (used for the Read the Docs PDF) cannot print these characters
# on its own. They appear in the bodies module docstring, a notebook title and
# the developer guide.
latex_elements = {
    "preamble": r"""
\usepackage[nointegrals]{wasysym}
\DeclareUnicodeCharacter{039B}{\ensuremath{\Lambda}}
\DeclareUnicodeCharacter{2609}{\astrosun}
\DeclareUnicodeCharacter{2641}{\earth}
\DeclareUnicodeCharacter{263E}{\leftmoon}
\DeclareUnicodeCharacter{263F}{\mercury}
\DeclareUnicodeCharacter{2640}{\venus}
\DeclareUnicodeCharacter{2642}{\mars}
\DeclareUnicodeCharacter{2643}{\jupiter}
\DeclareUnicodeCharacter{2644}{\saturn}
\DeclareUnicodeCharacter{26E2}{\uranus}
\DeclareUnicodeCharacter{2646}{\neptune}
\DeclareUnicodeCharacter{2647}{\pluto}
\DeclareUnicodeCharacter{1F38A}{}
\DeclareUnicodeCharacter{1F600}{:)}
""",
}

html_sidebars = {
    "**": [
        "about.html",
        "navigation.html",
        "relations.html",
        "searchbox.html",
        "donate.html",
    ]
}
