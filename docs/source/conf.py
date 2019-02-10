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

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "astropy": ("http://docs.astropy.org/en/stable/", None),
    "numpy": ("https://docs.scipy.org/doc/numpy/", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/reference", None),
    "matplotlib": ("http://matplotlib.org", None),
}

html_theme_options = {
    "logo": "logo_HD.png",
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
extensions = ["alabaster"]

templates_path = ["_templates"]

source_suffix = ".rst"

master_doc = "index"

html_theme = "alabaster"

html_theme_path = [alabaster.get_path()]

html_title = "EinsteinPy"

html_sidebars = {
    "**": [
        "about.html",
        "navigation.html",
        "relations.html",
        "searchbox.html",
        "donate.html",
    ]
}
